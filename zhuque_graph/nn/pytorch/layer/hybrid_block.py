from typing import List
import torch.nn as nn
from torch_geometric.nn.dense import Linear
from torch_geometric.nn import MessagePassing, GraphConv
import torch


class MLP(nn.Module):
    def __init__(self, channels: List[int], norm = None, bias = True, dropout = 0.):
        super().__init__()
        self.mlp = nn.ModuleList()
        for i in range(1, len(channels)):
            self.mlp.append(Linear(channels[i - 1], channels[i], bias=bias))
            # linear -> norm -> act -> dropout
            if i < len(channels) - 1:
                if norm and norm == 'batch':
                    self.mlp.append(nn.BatchNorm1d(channels[i]))
                elif norm and norm == 'layer':
                    self.mlp.append(nn.LayerNorm(channels[i]))
                elif norm and norm == 'instance':
                    self.mlp.append(nn.InstanceNorm1d(channels[i]))
                elif norm:
                    raise NotImplementedError(
                        f'Normalization layer "{norm}" not supported.')
                self.mlp.append(nn.ReLU())
                self.mlp.append(nn.Dropout(dropout))

    def reset_parameters(self):
        for moudle in self.mlp:
            if hasattr(moudle, "reset_paramters"):
                moudle.reset_parameters()

    def forward(self, x):
        for moudle in self.mlp:
            x = moudle(x)
        return x

class GINEConv(MessagePassing):
    """
    expansion: expansion factor of hidden channels in MLP layers.
    num_layers: num of mlp layers.
    """
    def __init__(self, emb_dim, num_layers = 2, norm = 'batch', bias = True, expansion = 1,
    dropout = 0.):
        super().__init__(aggr="add")
        channels = [emb_dim]
        for _ in range(num_layers - 1):
            channels.append(int(emb_dim * expansion))
        channels.append(emb_dim)
        self.mlp = MLP(channels, norm=norm, bias=bias, dropout=dropout)
        self.eps = nn.Parameter(torch.Tensor([0]))

    def reset_parameters(self):
        self.mlp.reset_parameters()

    def forward(self, x, edge_index, edge_attr):
        out = self.mlp((1 + self.eps) * x + self.propagate(edge_index, x=x, edge_attr=edge_attr))

        return out

    def message(self, x_j, edge_attr):
        return x_j if edge_attr is None else x_j + edge_attr

    def update(self, aggr_out):
        return aggr_out

class EdgeGraphConv(GraphConv):
    def message(self, x_j, edge_weight):
        return x_j if edge_weight is None else edge_weight * x_j

class HybridBlock(nn.Module):
    def __init__(self, hidden_channels, middle_channels, 
    dropout, num_radial, num_spherical, norm):
        super(HybridBlock, self).__init__()
        self.act = nn.ReLU()

        self.conv1 = EdgeGraphConv(hidden_channels, hidden_channels)
        self.conv2 = EdgeGraphConv(hidden_channels, hidden_channels)
        self.conv3 = GINEConv(hidden_channels, dropout=dropout, norm=norm)

        self.lin_feat1 = MLP([num_radial * num_spherical ** 2, middle_channels, hidden_channels],
        dropout=dropout)
        self.lin_feat2 = MLP([num_radial * num_spherical, middle_channels, hidden_channels],
        dropout=dropout)


        self.linear_cat = Linear(2 * hidden_channels, hidden_channels)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.conv3.reset_parameters()
        self.lin_feat1.reset_parameters()
        self.lin_feat2.reset_parameters()
        self.linear_cat.reset_parameters()

    def forward(self, x, feature1, feature2, pos_edge_index, edge_index, edge_attr):
        feature1 = self.lin_feat1(feature1)
        h1 = self.conv1(x, pos_edge_index, feature1)
        h1 = self.act(h1)

        feature2 = self.lin_feat2(feature2)
        h2 = self.conv2(x, pos_edge_index, feature2)
        h2 = self.act(h2)

        h = self.linear_cat(torch.cat([h1, h2], 1))
        h = self.act(h)

        h3 = self.conv3(x, edge_index, edge_attr)

        h = h + h3
        return h