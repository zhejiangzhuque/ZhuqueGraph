import numpy as np
import pandas as pd
import networkx as nx
from texttable import Texttable
from scipy.sparse import coo_matrix
import sys
import pickle as pkl
import scipy.sparse as sp
from networkx.readwrite import json_graph
import json
import os

def tab_printer(args):
    """
    Function to print the logs in a nice tabular format.
    :param args: Parameters used for the model.
    """
    args = vars(args)
    keys = sorted(args.keys())
    t = Texttable() 
    t.add_rows([["Parameter", "Value"]] +  [[k.replace("_"," ").capitalize(),args[k]] for k in keys])
    print(t.draw())

def construct_graph(edge_index):
    #graph = nx.Graph()
    #for i in range(edge_index.shape[-1]):
    #    graph.add_edge(edge_index[0][i], edge_index[1][i])
    edge_index=edge_index.swapaxes(0,1)
    graph = nx.from_edgelist(edge_index.tolist())
    return graph

def parse_index_file(filename):
    """Parse index file."""
    index = []
    for line in open(filename):
        index.append(int(line.strip()))
    return index


def graph_reader1(dataset_str):

    names = ['graph']
    objects = []
    for i in range(len(names)):
        with open("/INPUT/datasets/{0}/ind.{2}.{1}".format(dataset_str, names[i],dataset_str.replace("_ind","")), 'rb') as f:
            if sys.version_info > (3, 0):
                objects.append(pkl.load(f, encoding='latin1'))
            else:
                objects.append(pkl.load(f))

    #graph = tuple(objects)
    adj = nx.from_dict_of_lists(objects[0])
    return adj


def feature_reader1(dataset_str, compression=0):

    names = ['x',  'tx',  'allx']
    objects = []
    for i in range(len(names)):
        with open("/INPUT/datasets/{0}/ind.{2}.{1}".format(dataset_str, names[i],dataset_str.replace("_ind","")), 'rb') as f:
            if sys.version_info > (3, 0):
                objects.append(pkl.load(f, encoding='latin1'))
            else:
                objects.append(pkl.load(f))

    x, tx, allx  = tuple(objects)
    test_idx_reorder = parse_index_file("/INPUT/datasets/{0}/ind.{1}.test.index".format(dataset_str,dataset_str.replace("_ind","")))
    test_idx_range = np.sort(test_idx_reorder)

    if dataset_str == 'citeseer':
        # Fix citeseer dataset (there are some isolated nodes in the graph)
        # Find isolated nodes, add them as zero-vecs into the right position
        test_idx_range_full = range(min(test_idx_reorder), max(test_idx_reorder)+1)
        tx_extended = sp.lil_matrix((len(test_idx_range_full), x.shape[1]))
        tx_extended[test_idx_range-min(test_idx_range), :] = tx
        tx = tx_extended

    features = sp.vstack((allx, tx)).tolil()
    features[test_idx_reorder, :] = features[test_idx_range, :]
    preprocess_features(features)
    features = features.tocoo()
    features = features.toarray()

    feature_list = []
    feature_list.append(features)



    #features = sp.vstack((allx, tx)).tocoo()
    """
    features = sp.vstack((allx, tx)).tocoo()
    preprocess_features(features)
    features = features.toarray()
    """
    #features[test_idx_reorder, :] = features[test_idx_range, :]


    #features = coo_matrix((feature_values, (node_index, feature_index)), shape=(node_count, feature_count)).toarray()

    return features

def preprocess_features(features):
    """Row-normalize feature matrix and convert to tuple representation"""
    rowsum = np.array(features.sum(1))
    r_inv = np.power(rowsum, -1).flatten()
    r_inv[np.isinf(r_inv)] = 0.
    r_mat_inv = sp.diags(r_inv)
    features = r_mat_inv.dot(features)
    return features

def target_reader1(dataset_str):
    """
    Loads input data from gcn/data directory

    ind.dataset_str.x => the feature vectors of the training instances as scipy.sparse.csr.csr_matrix object;
    ind.dataset_str.tx => the feature vectors of the test instances as scipy.sparse.csr.csr_matrix object;
    ind.dataset_str.allx => the feature vectors of both labeled and unlabeled training instances
        (a superset of ind.dataset_str.x) as scipy.sparse.csr.csr_matrix object;
    ind.dataset_str.y => the one-hot labels of the labeled training instances as numpy.ndarray object;
    ind.dataset_str.ty => the one-hot labels of the test instances as numpy.ndarray object;
    ind.dataset_str.ally => the labels for instances in ind.dataset_str.allx as numpy.ndarray object;
    ind.dataset_str.graph => a dict in the format {index: [index_of_neighbor_nodes]} as collections.defaultdict
        object;
    ind.dataset_str.test.index => the indices of test instances in graph, for the inductive setting as list object.

    All objects above must be saved using python pickle module.

    :param dataset_str: Dataset name
    :return: All data input files loaded (as well the training/test data).
    """

    names = [ 'y',  'ty',  'ally']
    objects = []
    for i in range(len(names)):
        with open("/INPUT/datasets/{0}/ind.{2}.{1}".format(dataset_str, names[i],dataset_str.replace("_ind","")), 'rb') as f:
            if sys.version_info > (3, 0):
                objects.append(pkl.load(f, encoding='latin1'))
            else:
                objects.append(pkl.load(f))

    y, ty,  ally = tuple(objects)
    test_idx_reorder = parse_index_file("/INPUT/datasets/{0}/ind.{1}.test.index".format(dataset_str,dataset_str.replace("_ind","")))
    test_idx_range = np.sort(test_idx_reorder)

    #target = np.array(pd.read_csv(path)["target"]).reshape(-1, 1)

    if dataset_str == 'citeseer':
        # Fix citeseer dataset (there are some isolated nodes in the graph)
        # Find isolated nodes, add them as zero-vecs into the right position
        test_idx_range_full = range(min(test_idx_reorder), max(test_idx_reorder)+1)

        ty_extended = np.zeros((len(test_idx_range_full), y.shape[1]))
        ty_extended[test_idx_range-min(test_idx_range), :] = ty
        ty = ty_extended

    #labels = np.vstack((ally, ty)).reshape(-1,1)

    #labels = np.vstack((ally, ty))
    #labels = labels.reshape(-1, 1)
    #labels = np.append(ally, ty)

    ally = np.argmax(ally, axis=1)
    ty = np.argmax(ty, axis=1)
    labels = np.concatenate((ally, ty))
    labels[test_idx_reorder] = labels[test_idx_range]
    labels = labels.reshape(-1,1)


    """
    ally = np.argmax(ally, axis=1)
    ty = np.argmax(ty, axis=1)
    labels = np.concatenate((ally, ty), axis=0).reshape(-1,1)
    """
    return labels

