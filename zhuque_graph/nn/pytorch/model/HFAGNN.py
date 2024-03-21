import torch.nn as nn
from torch_geometric.nn.dense import Linear
from torch_geometric.nn import global_add_pool
import torch.nn.functional as F
import torch
from ogb.graphproppred.mol_encoder import AtomEncoder, BondEncoder
from zhuque_graph.nn.pytorch.layer.pos_feature_encoder import angle_emb, torsion_emb
from zhuque_graph.utils.dist_calc import dist_calc
from torch_cluster import radius_graph
from ..layer.hybrid_block import HybridBlock, MLP

class  HFAGNN(nn.Module):
    def __init__(self, cutoff=8.0, num_layers=4, hidden_channels=256, middle_channels=128, residual=True, 
    dropout=0., num_radial=3, num_spherical=2, out_channels=1, norm='layer'):
        super().__init__()
        self.num_layers = num_layers
        self.cutoff = cutoff
        self.residual = residual
        self.norms_list = nn.ModuleList()
        self.act = nn.ReLU()
        self.dropout = dropout

        self.x_emb = AtomEncoder(hidden_channels)
        self.edge_emb = BondEncoder(hidden_channels)

        self.feat1 = torsion_emb(num_radial=num_radial, num_spherical=num_spherical, cutoff=cutoff)
        self.feat2 = angle_emb(num_radial=num_radial, num_spherical=num_spherical, cutoff=cutoff)
    
        self.convs = nn.ModuleList()
        self.norms_list = nn.ModuleList() ### List of norms
        ### List of MLPs to transform virtual node at every layer
        self.mlp_virtualnode_list = nn.ModuleList()
        
        for layer in range(num_layers):
            self.convs.append(HybridBlock(hidden_channels, middle_channels,
                    dropout, num_radial, num_spherical, norm=norm))
            if norm and norm == 'batch':
                self.norms_list.append(nn.BatchNorm1d(hidden_channels))
            elif norm and norm == 'layer':
                self.norms_list.append(nn.LayerNorm(hidden_channels))
            elif norm and norm == 'instance':
                self.norms_list.append(nn.InstanceNorm1d(hidden_channels))
            elif norm:
                raise NotImplementedError(
                    f'Normalization layer "{norm}" not supported.')

            if layer < num_layers - 1:
                self.mlp_virtualnode_list.append(MLP([hidden_channels, middle_channels, hidden_channels],
                norm=norm))

        self.pool = global_add_pool
        self.graph_pred_linear = Linear(hidden_channels, out_channels)

        self.reset_parameters()

    def reset_parameters(self):
        for layer in range(self.num_layers):
            self.convs[layer].reset_parameters()
            if layer < self.num_layers - 1:
                self.mlp_virtualnode_list[layer].reset_parameters()
        self.graph_pred_linear.reset_parameters()

    def forward(self, batched_data):
        x, edge_index, edge_attr, batch, pos = batched_data.x, batched_data.edge_index, \
            batched_data.edge_attr, batched_data.batch, batched_data.pos

        num_nodes = x.size(0)

        x = self.x_emb(x)
        edge_attr = self.edge_emb(edge_attr)

        pos_edge_index = radius_graph(pos, r=self.cutoff, batch=batch)
        i, j = pos_edge_index[0], pos_edge_index[1]
        posi, posj = checkeq(pos[i], pos[j])
        vecs = posj - posi
        dist = vecs.norm(dim=-1)
        theta, phi, tau = dist_calc(self.cutoff, vecs, dist, i, j, num_nodes)
        feat1 = self.feat1(dist, theta, phi)
        feat2 = self.feat2(dist, tau)

        # first conv, no normalization
        h_in = self.convs[0](x, feat1, feat2, pos_edge_index, edge_index, edge_attr)
        h_virt = self.pool(h_in, batch)
        h_virt = self.mlp_virtualnode_list[0](h_virt)
        h = 0

        # norm -> act -> dropout -> conv -> res
        for layer in range(1, self.num_layers):
            # add message from virtual node to graph nodes.
            h_in = h_in + h_virt[batch]

            h = self.norms_list[layer](h_in)
            h = self.act(h)
            h = F.dropout(h, p=self.dropout, training=self.training)

            # virtual node conv
            h = self.convs[layer](h, feat1, feat2, pos_edge_index, edge_index, edge_attr)

            # add message from graph nodes to virtual node
            if layer < self.num_layers - 1:          
                h_virt_temp = h_virt + global_add_pool(h, batch)
                if self.residual:
                    h_virt = h_virt + self.mlp_virtualnode_list[layer](h_virt_temp)
                else:
                    h_virt = self.mlp_virtualnode_list[layer](h_virt_temp)
            
            # residual
            if self.residual:
                h = h + h_in
            h_in = h

        h = self.act(self.norms_list[0](h))
        h = F.dropout(h, p=self.dropout, training=self.training)

        h_graph = self.pool(h, batch)
        output = self.graph_pred_linear(h_graph)

        if self.training:
            return output
        else:
            # At inference time, we clamp the value between 0 and 20
            return torch.clamp(output, min=0, max=20)

# This function for sloving the problem of two atom have the same 3D position
# generated from rdkit (except case) 
def checkeq(posi, posj):
    eq_calcu = torch.eq(posi, posj).sum(-1)
    eq_indcies = torch.argwhere(eq_calcu == 3).squeeze(-1)
    if not torch.any(torch.isnan(eq_indcies)): 
        add_random = torch.randn(eq_indcies.shape[0], dtype=posj.dtype).uniform_(0.1000, 0.3000).to(posi.device)
        posi[eq_indcies,-1] += add_random
        posj[eq_indcies,-1] += torch.flip(add_random, dims=[0]).to(posi.device)
    return posi, posj

if __name__ == '__main__':
    model = HFAGNN()
    print(model)
