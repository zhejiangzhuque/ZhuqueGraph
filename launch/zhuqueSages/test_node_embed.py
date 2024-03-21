# import os.path as osp
import argparse
import pickle

from launch.zhuqueSages.planetoid import Planetoid


from zhuque_graph.utils.sages_classifier import Classifier



parser = argparse.ArgumentParser()

parser.add_argument('--dataset', type=str, default='PubMed')

args = parser.parse_args()


root_file = '/INPUT/datasets/' + args.dataset
dataset = Planetoid(root=root_file, name=args.dataset)
data = dataset[0]

idx_train, idx_val, idx_test = [], [], []
for i in range(data.x.size(0)):
    if data.train_mask[i]:
        idx_train.append(i)
    elif data.val_mask[i]:
        idx_val.append(i)
    elif data.test_mask[i]:
        idx_test.append(i)

args.dataset = args.dataset.lower()

with open('../data/node_embed/'+args.dataset+'_embed.pkl','rb') as file:
    save_embed = pickle.load(file)

node_embed = save_embed[args.dataset]

classifier = Classifier(vectors=node_embed)
acc, mic_f1s, macro_f1s = classifier(idx_train, idx_test, idx_val, data.y.numpy(), seed=0)
# db, cluster_acc, cluster_nmi, cluster_adj = clustering(Cluster, node_embed, data.y.numpy())
# return acc, mic_f1s, macro_f1s, db, cluster_acc, cluster_nmi, cluster_adj, node_embed

print('acc: {:.4f}, mic_f1: {:.4f}'.format(acc, mic_f1s))

