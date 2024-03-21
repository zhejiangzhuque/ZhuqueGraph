import math
import torch
from torch import nn
import torch.nn.functional as F
from torch.nn import Parameter
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.utils import remove_self_loops, add_self_loops, degree
from torch_scatter import scatter


def uniform(size, tensor):
    ''' sampled from the continuous uniform distribution '''
    bound = 1.0 / math.sqrt(size)
    tensor.data.uniform_(-bound, bound)


class SAGEConv(MessagePassing):
    def __init__(self, in_channels, out_channels, alphas=[0,1],
                 activate=False, bias=True, normalize=False, shared_weight=False,
                 aggr='mean', node_dim=-3,
                 **kwargs):
        super(SAGEConv, self).__init__(aggr=aggr, node_dim=node_dim, **kwargs)
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.alphas = alphas #[self_alpha, pro_alpha]
        self.activate = activate
        self.weight = Parameter(torch.Tensor(self.in_channels, out_channels))
        if shared_weight:
            self.self_weight = self.weight
        else:
            self.self_weight = Parameter(torch.Tensor(self.in_channels, out_channels))
        self.normalize = normalize
        if bias:
            self.bias = Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        uniform(self.in_channels, self.weight)
        uniform(self.in_channels, self.bias)
        uniform(self.in_channels,self.self_weight)

    def forward(self, x, edge_index, edge_weight=None, size=None):
        out  =  torch.matmul(x, self.self_weight)
        out2 = self.propagate(x, edge_index, edge_weight=edge_weight)
        # torch.Size([17499, 64, 1])
        # torch.Size([2, 102762])
        # torch.Size([102762])
        return self.alphas[0]*out + self.alphas[1]*out2

    def propagate(self, x, edge_index, edge_weight):
        out = self.message(x, edge_index, edge_weight)
        out = self.aggregate(out, edge_index)
        out = self.update(out)
        return out

    def message(self, x, edge_index, edge_weight):
        # x [node(genes), sample, in_channels]
        # edge_index [2, E]
        row, col = edge_index
        x_j = x[row] #source
        if edge_weight is None:
            return  x_j #torch.Size([102762, 64, 1])
        else:
            return edge_weight.view(-1,1,1) * x_j #torch.Size([102762, 64, 1])

    def aggregate(self, x_j, edge_index):
        row, col = edge_index
        aggr_out = scatter(x_j, col, dim=self.node_dim, reduce=self.aggr)
        return aggr_out


    def update(self, aggr_out):
        if self.activate:
            aggr_out = F.relu(aggr_out)
        aggr_out = torch.matmul(aggr_out, self.weight)
        if self.bias is not None:
            aggr_out += self.bias
        if self.normalize:
            aggr_out = F.normalize(aggr_out, p=2, dim=-1)
        return aggr_out


class scGraph(nn.Module):
    def __init__(self, in_channel, mid_channel, out_channel, num_nodes, num_edge,
                 global_conv1_dim, global_conv2_dim, FC1_dim, FC2_dim,
                 train_edges=True, dropout_ratio=None,
                 **args):
        super(scGraph,self).__init__()
        self.dropout_ratio = dropout_ratio
        if train_edges:
            self.edge_weight = nn.Parameter(torch.ones(num_edge).float() * 0.01)
        else:
            self.edge_weight = None

        # 1. Graph representation module
        self.graph_conv = SAGEConv(in_channel, mid_channel)
        self.bn1 = torch.nn.LayerNorm((num_nodes, mid_channel))
        self.act = nn.ReLU()

        # 2. Feature extraction module
        self.global_conv1 = torch.nn.Conv2d(mid_channel, global_conv1_dim,1)
        self.global_bn1 = torch.nn.BatchNorm2d(global_conv1_dim)

        self.global_conv2 = torch.nn.Conv2d(global_conv1_dim, global_conv2_dim, [1,1])
        self.global_bn2 = torch.nn.BatchNorm2d(global_conv2_dim)

        self.nn = []
        channel_list = [global_conv2_dim * num_nodes, FC1_dim, FC2_dim]
        for idx, num in enumerate(channel_list[:-1]):
            self.nn.append(nn.Linear(channel_list[idx], channel_list[idx + 1]))
            self.nn.append(nn.BatchNorm1d(channel_list[idx + 1]))
            if self.dropout_ratio > 0:
                self.nn.append(nn.Dropout(0.3))
            self.nn.append(nn.ReLU())
        self.global_fc_nn = nn.Sequential(*self.nn)

        # 3. Classification module
        self.fc1 = nn.Linear(FC2_dim,out_channel)
        self.reset_parameters()

    def init_weights(self, m):
        if type(m) == nn.Linear:
            nn.init.xavier_uniform_(m.weight)
            m.bias.data.fill_(0.01)
        if type(m) == torch.nn.Conv2d:
            nn.init.kaiming_normal_(m.weight, mode='fan_out')
            uniform(m.weight.shape[1], m.bias)

    def reset_parameters(self,):
        self.graph_conv.apply(self.init_weights)
        self.global_conv1.apply(self.init_weights)
        self.global_conv2.apply(self.init_weights)
        self.global_fc_nn.apply(self.init_weights)
        self.fc1.apply(self.init_weights)
        pass

    # def gcn_weight_penalty(self, mode='L2'):
    #     if mode == 'L1':
    #         func = lambda x:  torch.sum(torch.abs(x))
    #     elif mode =='L2':
    #         func  = lambda x: torch.sqrt(torch.sum(x**2))
    #     tmp = getattr(self.graph_conv, 'weight',None)
    #     if tmp is not None:
    #         loss = func(tmp)
    #     return loss

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
            # x: torch.Size([23459, 64, 1])
            # edge_index: torch.Size([2, 109914])
        # whether the edge is trainable
        if self.edge_weight is not None:
            edge_weight = torch.sigmoid(self.edge_weight) # edge_weight: torch.Size([109914])
        else:
            edge_weight = None

        x = self.graph_conv(x, edge_index, edge_weight=edge_weight)
        x = self.act(x)
        # LayerNorm
        x = x.permute(1, 0, 2)  # #samples x #nodes x #features
        x = self.bn1(x)
        x = x.permute(1, 0, 2)  # #nodes x #samples x #features

        if self.dropout_ratio is not None:
            x = F.dropout(x, p=0.1, training=self.training)

        x = x.permute(1,2,0)  # #samples x #features x #nodes
        x = x.unsqueeze(dim=-1) # #samples x #features x #nodes x 1

        x = self.global_conv1(x)  # #samples x #features x #nodes x 1
        x = self.act(x)
        x = self.global_bn1(x)
        if self.dropout_ratio >0: x = F.dropout(x, p=0.3, training=self.training)

        x = self.global_conv2(x)
        x = self.act(x)
        x = self.global_bn2(x)
        if self.dropout_ratio >0: x = F.dropout(x, p=0.3, training=self.training)

        x = x.squeeze(dim=-1)  # #samples  x #features  x #nodes
        num_samples = x.shape[0]

        # flatten
        x = x.view(num_samples,-1)

        x = self.global_fc_nn(x)
        x = self.fc1(x)

        return F.softmax(x, dim=-1), x #prob, embedding



