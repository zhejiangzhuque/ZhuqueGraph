import datetime
import os
import time

import torch
import torch.nn.functional as F
import torchmetrics.functional as MF
from torch.utils.data import DataLoader
from torch.utils.data.distributed import DistributedSampler
from wg_torch import comm as comm
from wg_torch import embedding_ops as embedding_ops
from wg_torch import graph_ops as graph_ops
from wg_torch.wm_tensor import *

from wholegraph.torch import wholegraph_pytorch as wg


def parse_max_neighbors(num_layer, neighbor_str):
    neighbor_str_vec = neighbor_str.split(",")
    max_neighbors = []
    for ns in neighbor_str_vec:
        max_neighbors.append(int(ns))
    assert len(max_neighbors) == 1 or len(max_neighbors) == num_layer
    if len(max_neighbors) != num_layer:
        for i in range(1, num_layer):
            max_neighbors.append(max_neighbors[0])
    # max_neighbors.reverse()
    return max_neighbors

def create_gnn_layers(in_feat_dim, hidden_feat_dim, num_layer, num_head, args):
    if args.framework == "dgl":
        import dgl
        from dgl.nn.pytorch.conv import SAGEConv, GATConv
    elif args.framework == "pyg":
        from torch_sparse import SparseTensor
        from torch_geometric.nn import SAGEConv, GATConv
    elif args.framework == "wg":
        from wg_torch.gnn.SAGEConv import SAGEConv
        from wg_torch.gnn.GATConv import GATConv
    gnn_layers = torch.nn.ModuleList()
    for i in range(num_layer):
        layer_output_dim = hidden_feat_dim // num_head if i != num_layer - 1 else args.classnum
        layer_input_dim = in_feat_dim if i == 0 else hidden_feat_dim
        mean_output = True if i == num_layer - 1 else False
        if args.framework == "pyg":
            if args.model == "sage":
                gnn_layers.append(SAGEConv(layer_input_dim, layer_output_dim))
            elif args.model == "gat":
                concat = not mean_output
                gnn_layers.append(
                    GATConv(
                        layer_input_dim, layer_output_dim, heads=num_head, concat=concat
                    )
                )
            else:
                assert args.model == "gcn"
                gnn_layers.append(
                    SAGEConv(layer_input_dim, layer_output_dim, root_weight=False)
                )
        elif args.framework == "dgl":
            if args.model == "sage":
                gnn_layers.append(SAGEConv(layer_input_dim, layer_output_dim, "mean"))
            elif args.model == "gat":
                gnn_layers.append(
                    GATConv(
                        layer_input_dim,
                        layer_output_dim,
                        num_heads=num_head,
                        allow_zero_in_degree=True,
                    )
                )
            else:
                assert args.model == "gcn"
                gnn_layers.append(SAGEConv(layer_input_dim, layer_output_dim, "gcn"))
        elif args.framework == "wg":
            if args.model == "sage":
                gnn_layers.append(SAGEConv(layer_input_dim, layer_output_dim))
            elif args.model == "gat":
                gnn_layers.append(
                    GATConv(
                        layer_input_dim,
                        layer_output_dim,
                        num_heads=num_head,
                        mean_output=mean_output,
                    )
                )
            else:
                assert args.model == "gcn"
                gnn_layers.append(
                    SAGEConv(layer_input_dim, layer_output_dim, aggregator="gcn")
                )
    return gnn_layers

def create_sub_graph(
    target_gid,
    target_gid_1,
    edge_data,
    csr_row_ptr,
    csr_col_ind,
    sample_dup_count,
    add_self_loop: bool,
    args
):
    if args.framework == "pyg":
        neighboor_dst_unique_ids = csr_col_ind
        neighboor_src_unique_ids = edge_data[1]
        target_neighbor_count = target_gid.size()[0]
        if add_self_loop:
            self_loop_ids = torch.arange(
                0,
                target_gid_1.size()[0],
                dtype=neighboor_dst_unique_ids.dtype,
                device=target_gid.device,
            )
            edge_index = SparseTensor(
                row=torch.cat([neighboor_src_unique_ids, self_loop_ids]).long(),
                col=torch.cat([neighboor_dst_unique_ids, self_loop_ids]).long(),
                sparse_sizes=(target_gid_1.size()[0], target_neighbor_count),
            )
        else:
            edge_index = SparseTensor(
                row=neighboor_src_unique_ids.long(),
                col=neighboor_dst_unique_ids.long(),
                sparse_sizes=(target_gid_1.size()[0], target_neighbor_count),
            )
        return edge_index
    elif args.framework == "dgl":
        if add_self_loop:
            self_loop_ids = torch.arange(
                0,
                target_gid_1.numel(),
                dtype=edge_data[0].dtype,
                device=target_gid.device,
            )
            block = dgl.create_block(
                (
                    torch.cat([edge_data[0], self_loop_ids]),
                    torch.cat([edge_data[1], self_loop_ids]),
                ),
                num_src_nodes=target_gid.size(0),
                num_dst_nodes=target_gid_1.size(0),
            )
        else:
            block = dgl.create_block(
                (edge_data[0], edge_data[1]),
                num_src_nodes=target_gid.size(0),
                num_dst_nodes=target_gid_1.size(0),
            )
        return block
    else:
        assert args.framework == "wg"
        return [csr_row_ptr, csr_col_ind, sample_dup_count]
    return None

def layer_forward(layer, x_feat, x_target_feat, sub_graph, args):
    if args.framework == "pyg":
        x_feat = layer((x_feat, x_target_feat), sub_graph)
    elif args.framework == "dgl":
        x_feat = layer(sub_graph, (x_feat, x_target_feat))
    elif args.framework == "wg":
        x_feat = layer(sub_graph[0], sub_graph[1], sub_graph[2], x_feat, x_target_feat)
    return x_feat

class NodeClassificationWgGNNModel(torch.nn.Module):
    def __init__(
        self,
        graph: graph_ops.HomoGraph,
        num_layer,
        hidden_feat_dim,
        max_neighbors: str,
        args
    ):
        super().__init__()
        self.args = args
        self.graph = graph
        self.num_layer = num_layer
        self.hidden_feat_dim = hidden_feat_dim
        self.max_neighbors = parse_max_neighbors(num_layer, max_neighbors)
        num_head = args.heads if (args.model == "gat") else 1
        assert hidden_feat_dim % num_head == 0
        in_feat_dim = self.graph.node_feat_shape()[1]
        self.gnn_layers = create_gnn_layers(
            in_feat_dim, hidden_feat_dim, num_layer, num_head, args
        )
        self.mean_output = True if args.model == "gat" else False
        self.add_self_loop = True if args.model == "gat" else False
        self.gather_fn = embedding_ops.EmbeddingLookUpModule(need_backward=False)

    def forward(self, ids):
        ids = ids.to(self.graph.id_type()).cuda()
        (
            target_gids,
            edge_indice,
            csr_row_ptrs,
            csr_col_inds,
            sample_dup_counts,
        ) = self.graph.unweighted_sample_without_replacement(ids, self.max_neighbors)
        x_feat = self.gather_fn(target_gids[0], self.graph.node_feat)
        # x_feat = self.graph.gather(target_gids[0])
        for i in range(self.num_layer):
            x_target_feat = x_feat[: target_gids[i + 1].numel()]
            sub_graph = create_sub_graph(
                target_gids[i],
                target_gids[i + 1],
                edge_indice[i],
                csr_row_ptrs[i],
                csr_col_inds[i],
                sample_dup_counts[i],
                self.add_self_loop,
                self.args
            )
            x_feat = layer_forward(self.gnn_layers[i], x_feat, x_target_feat, sub_graph, self.args)
            if i != self.num_layer - 1:
                if self.args.framework == "dgl":
                    x_feat = x_feat.flatten(1)
                x_feat = F.relu(x_feat)
                x_feat = F.dropout(x_feat, self.args.dropout, training=self.training)
        if self.args.framework == "dgl" and self.mean_output:
            out_feat = x_feat.mean(1)
        else:
            out_feat = x_feat
        return out_feat
