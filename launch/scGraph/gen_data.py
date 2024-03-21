import numpy as np
import pandas as pd
import argparse
import time
from sklearn import preprocessing
import torch
from torch_geometric.utils import  add_remaining_self_loops


def get_path(data):
    expr_dir = f'./{data}/data_df.csv'
    net_dir = f'./{data}/graph_df.csv'
    label_dir = f'./{data}/label_df.csv'
    out_dir = f'./{data}/dataset.npz'
    return expr_dir, net_dir, label_dir, out_dir


if __name__ == '__main__':
    start =time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', default='Demo', type=str, help='Baron_Human, Demo')
    parser.add_argument('--quantile', default='0.99', type=float)
    args = parser.parse_args()

    expr_dir, net_dir, label_dir, out_dir = get_path(args.data)

    # load data：expr_mat， inter_network， cell_label
    # 1. expr_mat
    data_df = pd.read_csv(expr_dir, index_col=0) #(1000, 23459) cell x gene
    # normalize by dividing it by its total expression & multiply by a scale factor 10e6
    data_df = data_df.apply(lambda x: x/x.sum() * 1e6, axis=1)
    # add a pseudo count & apply log2 transformation
    data_df = np.log2(data_df+1)

    gene = data_df.columns.values
    gene = [int(i) for i in gene]
    barcode  = data_df.index.values #cell (1000,)


    # 2. inter_network
    graph_df = pd.read_csv(net_dir, index_col=0)  # (10831388, 3)
    # filter out edges have low score
    thr = graph_df.score.quantile(args.quantile) #957
    graph_df = graph_df[graph_df.score >= thr] #(109915, 3)
    edge_index = graph_df[['node1', 'node2']].values.T #(2, 109916)
    # add remaining self-interaction edge (convert gene to gene idx)
    edge_index = torch.tensor(edge_index)
    edge_index = add_remaining_self_loops(edge_index)[0].numpy() #(2, 133373)

    # remove isolate gene from data_df
    isolate_gene = list(set(gene) - set(edge_index[0]) - set(edge_index[1]))
    isolate_gene = [str(i) for i in isolate_gene]
    data_df = data_df.drop(isolate_gene, axis=1)


    # 3. label
    if args.data == 'Demo':
        label_df = pd.read_csv(label_dir, header=0, index_col=0)  # cell_label (1000, 1)
    else:
        label_df = pd.read_csv(label_dir, header=0, index_col=0, names=['label'])  # cell_label (8569, 1)
    label_cat = list(np.unique(label_df.values))  #cell type 5
    # convert label to number
    le = preprocessing.LabelEncoder()
    le.fit(label_df.label)
    label_df.label = le.transform(label_df.label)
    label = label_df.label.tolist()


    # save pre-processed data
    data_dict = {}
    data_dict['logExpr'] = data_df.values  #(1000, 23459)
    data_dict['gene'] = gene  #gene_idx (23459,)
    data_dict['barcode'] = barcode  #(1000,)
    data_dict['edge_index'] = edge_index  #(2, 133373)
    data_dict['label_cat'] = label_cat  # 5
    data_dict['label'] = label  # 1000

    np.savez(out_dir, **data_dict)
    print('Finished.')
    end = time.time()
    print("Total time = ", end - start) # 31s -> 6s





