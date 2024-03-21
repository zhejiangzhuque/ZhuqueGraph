import pickle
import argparse
import numpy as np
from sklearn.manifold import TSNE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', type=str, default="cora", choices=['cora', 'citeseer', 'pubmed'])
args = parser.parse_args()


with open('Cora_dgi_node_embed512_epoch_100_label.pkl','rb') as file:
    graph_dic = pickle.load(file)

node_embed, node_label = graph_dic['embed'], graph_dic['label']


tsne = TSNE(n_components=2, random_state=0)
X_2d = tsne.fit_transform(node_embed)
target_ids = range(max(node_label) + 1)

plt.figure(figsize=(6, 6))

if args.dataset=='Cora':
    colors = 'r', 'g', 'b', 'c', 'm', 'y', 'k'
elif args.dataset=='CiteSeer':
    colors = 'r', 'g', 'b', 'c', 'm', 'y'
plt.axis('off')

for i, c in zip(target_ids, colors):
    plt.scatter(X_2d[node_label == i, 0], X_2d[node_label == i, 1], c=c)

plt.savefig(args.dataset+"_dgi_emb512_epoch100_tsne.jpg")