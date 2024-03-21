# import os.path as osp
import argparse
import torch
import numpy as np
import random
from launch.zhuqueSages.planetoid import Planetoid
from zhuque_graph.nn.pytorch.model.sages_nn import SAGES
from torch_geometric.nn import GCNConv, GAE
from sklearn.cluster import SpectralClustering
from zhuque_graph.utils.sages_classifier import Classifier


def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True

# 设置随机数种子
setup_seed(1234)

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='SAGES')
parser.add_argument('--data_name', type=str, default='Cora')
parser.add_argument('--channels', type=int, default=512)
parser.add_argument('--num_layers', type=int, default=5)

args = parser.parse_args()
assert args.model in ['GAE', 'VGAE', 'SAGES']
assert args.data_name in ['cora', 'citeseer', 'pubmed']
kwargs = {'GAE': GAE, 'SAGES':SAGES}

root_file = '/INPUT/datasets/' + args.data_name
dataset = Planetoid(root=root_file, name=args.data_name)
data = dataset[0]

idx_train, idx_val, idx_test = [], [], []
for i in range(data.x.size(0)):
    if data.train_mask[i]:
        idx_train.append(i)
    elif data.val_mask[i]:
        idx_val.append(i)
    elif data.test_mask[i]:
        idx_test.append(i)

# 设置谱聚类，进行测试
if args.data_name == 'cora':
    n_clusters = 7
    Cluster = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=0)
elif args.data_name == 'citeseer':
    n_clusters = 6
    Cluster = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=0)
elif args.data_name == 'pubmed':
    n_clusters = 3
    Cluster = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=0)

class Encoder(torch.nn.Module):
    # layer_dims包括初始feature维度，如cora[1433, 512, 512]等
    def __init__(self, layer_dims):
        super(Encoder, self).__init__()

        self.num_layers = len(layer_dims) - 1
        self.convs = torch.nn.ModuleList()

        for i in range(self.num_layers):
            self.convs.append(GCNConv(layer_dims[i], layer_dims[i+1], cached=True))


    # 调试：对隐层Z加入归一化，以保证其可以进行谱聚类
    def scale(self, z):
        zmax = z.max(dim=1, keepdim=True)[0]
        zmin = z.min(dim=1, keepdim=True)[0]
        z_std = (z - zmin) / (zmax - zmin)
        z_scaled = z_std
        return z_scaled

    def forward(self, x, edge_index):
        for i in range(self.num_layers):
            x = self.convs[i](x, edge_index)

        out = x
        # out = self.scale(out)
        # out = F.normalize(out)
        return out

# 这里是DGI的腐蚀函数
def corruption(x, edge_index):
    return x[torch.randperm(x.size(0))], edge_index



channels = args.channels
layer_dims = [dataset.num_features, channels, channels]
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


model = SAGES(feat_channels=dataset.num_features, hidden_channels=channels, encoder=Encoder(layer_dims),
              summary=lambda z, *args, **kwargs: torch.sigmoid(z.mean(dim=0)),
              corruption=corruption).to(device)

x, train_pos_edge_index = data.x.to(device), data.edge_index.to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)


def train():
    model.train()
    optimizer.zero_grad()
    pos_z, neg_z, summary = model(x, train_pos_edge_index)
    loss = model.loss(pos_z, neg_z, summary, train_pos_edge_index, x, 1, 0.01, 0.01)

    loss.backward()
    optimizer.step()



def test():
    model.eval()
    with torch.no_grad():
        z, _, _ = model(data.x.to(device), data.edge_index.to(device))

    node_embed = z.cpu().numpy()
    classifier = Classifier(vectors=node_embed)
    acc, mic_f1s, macro_f1s = classifier(idx_train, idx_test, idx_val, data.y.numpy(), seed=0)
    # db, cluster_acc, cluster_nmi, cluster_adj = clustering(Cluster, node_embed, data.y.numpy())
    # return acc, mic_f1s, macro_f1s, db, cluster_acc, cluster_nmi, cluster_adj, node_embed
    return acc, mic_f1s, macro_f1s


best_acc, best_c_acc, best_c_nmi, best_c_adj = 0, 0, 0, 0
for epoch in range(1, 200):
    train()
    # acc, mic_f1s, macro_f1s, db, cluster_acc, cluster_nmi, cluster_adj, node_embed = test()
    # print('Epoch: {:03d}, acc: {:.4f}, mic_f1: {:.4f}， c_acc: {:.4f}, c_nmi: {:.4f}, c_adj: {:.4f}'.format(epoch, acc, mic_f1s, cluster_acc, cluster_nmi, cluster_adj))
    # best_acc, best_c_acc, best_c_nmi, best_c_adj = max(acc, best_acc), max(cluster_acc,best_c_acc), max(cluster_nmi, best_c_nmi), max(best_c_adj, cluster_adj)
    # 只进行节点分类
    acc, mic_f1s, macro_f1s = test()
    print('Epoch: {:03d}, acc: {:.4f}, mic_f1: {:.4f}'.format(epoch, acc, mic_f1s))

    best_acc = max(acc, best_acc)


# print("Best ACC: {:.4f}, Cluster_acc: {:.4f}, NMI: {:.4f}, Adj: {:.4f}".format(best_acc, best_c_acc, best_c_nmi, best_c_adj))
print("Best ACC: {:.4f}".format(best_acc))