# coding=utf-8
# Copyright 2019 The Google Research Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Collections of preprocessing functions for different graph formats."""

import json
import time

from networkx.readwrite import json_graph
import numpy as np
import metis
import scipy.sparse as sp
import sklearn.metrics
import sklearn.preprocessing
import tensorflow.compat.v1 as tf
from tensorflow.compat.v1 import gfile
import networkx as nx
import sys
import pickle as pkl
import pandas as pd


def partition_graph(adj, idx_nodes, num_clusters):
    """partition a graph by METIS."""

    start_time = time.time()
    num_nodes = len(idx_nodes)  # number of visiable nodes
    num_all_nodes = adj.shape[0]

    neighbor_intervals = []  # indptr for csr-matrix
    neighbors = []  # record neighbors for each node in idx , 1-dimension
    edge_cnt = 0        # total number of neighbors of all nodes in idx
    neighbor_intervals.append(0)
    # lookup sub-adj-matrix for idx, then convert to lil
    train_adj_lil = adj[idx_nodes, :][:, idx_nodes].tolil()
    train_ord_map = dict()  # dict map original-id to sub-graph-id
    # record neighbors for each node in idx , 2-dimension ; indice for
    # csr-matrix
    train_adj_lists = [[] for _ in range(num_nodes)]
    for i in range(num_nodes):
        rows = train_adj_lil[i].rows[0]
        # self-edge needs to be removed for valid format of METIS
        if i in rows:
            rows.remove(i)
        train_adj_lists[i] = rows
        neighbors += rows
        edge_cnt += len(rows)
        neighbor_intervals.append(edge_cnt)
        train_ord_map[idx_nodes[i]] = i

    if num_clusters > 1:
        # group contains group-id for each node
        _, groups = metis.part_graph(train_adj_lists, num_clusters, seed=1)
    else:
        groups = [0] * num_nodes  # all nodes in the same group:"0"

    part_row = []
    part_col = []
    part_data = []
    # a list for a cluster, contains the original-ids in this cluster
    parts = [[] for _ in range(num_clusters)]
    for nd_idx in range(num_nodes):
        gp_idx = groups[nd_idx]
        nd_orig_idx = idx_nodes[nd_idx]
        parts[gp_idx].append(nd_orig_idx)
        for nb_orig_idx in adj[nd_orig_idx].indices:
            nb_idx = train_ord_map[nb_orig_idx]
            # only if the neighbor is in the same group with current node that
            # this edge will be added into part_adj
            if groups[nb_idx] == gp_idx:
                part_data.append(1)
                part_row.append(nd_orig_idx)
                part_col.append(nb_orig_idx)
    part_data.append(0)
    part_row.append(num_all_nodes - 1)
    part_col.append(num_all_nodes - 1)
    part_adj = sp.coo_matrix((part_data, (part_row, part_col))).tocsr()
    # part_adj =sp.csr_matrix((part_data,(part_row,part_col)),shape=(num_all_nodes,num_all_nodes))

    tf.logging.info('Partitioning done. %f seconds.', time.time() - start_time)
    return part_adj, parts


def partition_graph_weight(adj, idx_nodes, num_clusters,weight_reg = 5,sym_prop=True,sample_rate = -1):
    """partition a graph by METIS."""
    # weight_regularization: convert edge weight from (0,1] to 0~ weight_reg

    # regular the full adj
    # 12.26Ver
    if sym_prop:
        adj_prop = adj+adj.transpose()
    else:
        adj_prop =adj
    adj_coo =(adj+adj.transpose()).tocoo()


    # if sym_prop:
    #     adj_coo = (adj+adj.transpose()).tocoo()
    # else:
    #     adj_coo = adj.tocoo()

    adj_data_raw = np.array(adj_coo.data,dtype=np.float32)
    if np.max(adj_data_raw)>np.min(adj_data_raw):
        adj_data_raw = (adj_data_raw-np.min(adj_data_raw)+0.0)/(np.max(adj_data_raw)-np.min(adj_data_raw)+0.0)
        adj_data = np.array(adj_data_raw*weight_reg+1,dtype = np.int32)
    else:
        adj_data=np.ones_like(adj_data_raw,dtype=np.int)
    adj_row = np.array(adj_coo.row)
    adj_col = np.array(adj_coo.col)

    if sample_rate>0:
        # Naive sample without considering the symmetric matrix
        # sample_pos = np.random.choice(np.size(adj_data),int(sample_rate*np.size(adj_data)),replace=False)
        # adj_data = adj_data[sample_pos]
        # adj_row = adj_row[sample_pos]
        # adj_col = adj_col[sample_pos]

        adj_all = np.vstack([adj_row,adj_col,adj_data]) # 3* 2M
        adj_all = adj_all.transpose()[adj_all[0]>adj_all[1]] # M*3
        sample_pos = np.random.choice(adj_all.shape[0],int(sample_rate*adj_all.shape[0]),replace=False)
        adj_all = adj_all[sample_pos].transpose() #  3 * 0.9M
        adj_rev = adj_all.copy()
        adj_rev[[0,1],:] = adj_rev[[1,0],:]
        adj_all = np.hstack([adj_all,adj_rev])
        adj_row = adj_all[0].flatten()
        adj_col = adj_all[1].flatten()
        adj_data = adj_all[2].flatten()

    adj_reg = sp.csr_matrix((adj_data,(adj_row,adj_col)),shape=adj_coo.shape)

    start_time = time.time()
    num_nodes = len(idx_nodes)  # number of visiable nodes
    num_all_nodes = adj.shape[0]

    neighbor_intervals = []  # indptr for csr-matrix
    neighbors = []  # record neighbors for each node in idx , 1-dimension
    edge_cnt = 0        # total number of neighbors of all nodes in idx
    neighbor_intervals.append(0)
    # lookup sub-adj-matrix for idx,  then convert to lil
    train_adj_lil = adj_reg[idx_nodes, :][:, idx_nodes].tolil()
    # dict map original-id to sub-graph-id
    train_ord_map = dict()
    # record neighbors for each node in idx , 2-dimension ; indice for
    # csr-matrix
    train_adj_lists_weight = [[] for _ in range(num_nodes)]
    for i in range(num_nodes):
        rows = train_adj_lil[i].rows[0]
        weights = train_adj_lil[i].data[0]
        # self-edge needs to be removed for valid format of METIS

        if i in rows:
            weights.pop(rows.index(i))
            rows.remove(i)
        row_w_zip = zip(rows,weights)
        train_adj_lists_weight[i] = list(row_w_zip)
        neighbors += rows
        edge_cnt += len(rows)
        neighbor_intervals.append(edge_cnt)
        train_ord_map[idx_nodes[i]] = i

    if num_clusters > 1:
        # group contains group-id for each node
        _, groups = metis.part_graph(train_adj_lists_weight, num_clusters, seed=1)
    else:
        groups = [0] * num_nodes  # all nodes in the same group:"0"

    part_row = []
    part_col = []
    part_data = []
    # a list for a cluster, contains the original-ids in this cluster
    parts = [[] for _ in range(num_clusters)]

    #adj_prop=adj_reg
    # if sym_prop:
    #     adj_prop = adj_reg
    # else:
    #     adj_prop = adj
    # for each visiable node
    for nd_idx in range(num_nodes):
        gp_idx = groups[nd_idx]
        nd_orig_idx = idx_nodes[nd_idx]
        parts[gp_idx].append(nd_orig_idx)

        indices_now = adj_prop[nd_orig_idx].indices
        data_now = adj_prop[nd_orig_idx].data
        for ii in range(len(indices_now)):
            nb_orig_idx = indices_now[ii]
            nb_idx = train_ord_map[nb_orig_idx]
            if groups[nb_idx] == gp_idx:
                part_data.append(data_now[ii])
                part_row.append(nd_orig_idx)
                part_col.append(nb_orig_idx)

    # part_data.append(0)
    # part_row.append(num_all_nodes - 1)
    # part_col.append(num_all_nodes - 1)
    # part_adj = sp.coo_matrix((part_data, (part_row, part_col))).tocsr()
    part_adj =sp.csr_matrix((part_data,(part_row,part_col)),shape=(num_all_nodes,num_all_nodes))
    # part_adj is also symmetric

    tf.logging.info('Partitioning done. %f seconds.', time.time() - start_time)
    return part_adj, parts


def parse_index_file(filename):
    """Parse index file."""
    index = []
    for line in gfile.Open(filename):
        index.append(int(line.strip()))
    return index


def sample_mask(idx, l):
    """Create mask."""
    mask = np.zeros(l)
    mask[idx] = 1
    return np.array(mask, dtype=np.bool)


def sym_normalize_adj(adj):
    """Normalization by D^{-1/2} (A+I) D^{-1/2}."""
    adj = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj.sum(1)) + 1e-20
    d_inv_sqrt = np.power(rowsum, -0.5).flatten()
    d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.
    d_mat_inv_sqrt = sp.diags(d_inv_sqrt, 0)
    adj = adj.dot(d_mat_inv_sqrt).transpose().dot(d_mat_inv_sqrt)
    return adj


def normalize_adj(adj):
    rowsum = np.array(adj.sum(1)).flatten()
    d_inv = 1.0 / (np.maximum(1.0, rowsum))
    d_mat_inv = sp.diags(d_inv, 0)
    adj = d_mat_inv.dot(adj)
    return adj


def normalize_adj_diag_enhance(adj, diag_lambda):
    """Normalization by  A'=(D+I)^{-1}(A+I), A'=A'+lambda*diag(A')."""
    adj = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj.sum(1)).flatten()
    d_inv = 1.0 / (rowsum + 1e-20)
    d_mat_inv = sp.diags(d_inv, 0)
    adj = d_mat_inv.dot(adj)
    adj = adj + diag_lambda * sp.diags(adj.diagonal(), 0)
    return adj


def sparse_to_tuple(sparse_mx):
    """Convert sparse matrix to tuple representation."""

    def to_tuple(mx):
        if not sp.isspmatrix_coo(mx):
            mx = mx.tocoo()
        coords = np.vstack((mx.row, mx.col)).transpose()
        values = mx.data
        shape = mx.shape
        return coords, values, shape

    if isinstance(sparse_mx, list):
        for i in range(len(sparse_mx)):
            sparse_mx[i] = to_tuple(sparse_mx[i])
    else:
        sparse_mx = to_tuple(sparse_mx)

    return sparse_mx


def tuple_to_sparse(a):
    if isinstance(a, sp.coo_matrix):
        return a
    elif isinstance(a, tuple):
        coords, values, shape = a
        row = coords[:, 0]
        col = coords[:, 1]
        return sp.coo_matrix((values, (row, col)), shape=shape)
    else:
        assert 0, "Invalid Input"


def calc_f1(y_pred, y_true, multilabel):
    if multilabel:
        y_pred[y_pred > 0] = 1
        y_pred[y_pred <= 0] = 0
    else:
        y_true = np.argmax(y_true, axis=1)
        y_pred = np.argmax(y_pred, axis=1)
    return sklearn.metrics.f1_score(
        y_true, y_pred, average='micro'), sklearn.metrics.f1_score(
            y_true, y_pred, average='macro')

def calc_acc(y_pred, y_true, multilabel):
    if multilabel:
        y_pred[y_pred > 0] = 1
        y_pred[y_pred <= 0] = 0
    else:
        y_true = np.argmax(y_true, axis=1)
        y_pred = np.argmax(y_pred, axis=1)
    return sklearn.metrics.accuracy_score(y_true,y_pred)

def construct_feed_dict(features, support, labels, labels_mask, nodes_weight,placeholders):
    """Construct feed dictionary."""
    feed_dict = dict()
    feed_dict.update({placeholders['labels']: labels})
    feed_dict.update({placeholders['labels_mask']: labels_mask})
    feed_dict.update({placeholders['features']: features})
    feed_dict.update({placeholders['support']: support})
    feed_dict.update({placeholders['num_features_nonzero']: features[1].shape})
    feed_dict.update({placeholders['nodes_weight']:nodes_weight})
    return feed_dict


def construct_feed_dict_original(features, support, labels, labels_mask, placeholders):
    """Construct feed dictionary."""
    feed_dict = dict()
    feed_dict.update({placeholders['labels']: labels})
    feed_dict.update({placeholders['labels_mask']: labels_mask})
    feed_dict.update({placeholders['features']: features})
    feed_dict.update({placeholders['support']: support})
    feed_dict.update({placeholders['num_features_nonzero']: features[1].shape})
    return feed_dict


def construct_feed_dict_hc(
        features,
        support,
        labels,
        labels_mask,
        placeholders,
        hier_layer,
        features_list,
        support_list,
        lookup_list,
        placeholders_hc):
    """Construct feed dictionary."""
    feed_dict = dict()
    feed_dict.update({placeholders['labels']: labels})
    feed_dict.update({placeholders['labels_mask']: labels_mask})
    feed_dict.update({placeholders['features']: features})
    feed_dict.update({placeholders['support']: support})
    feed_dict.update({placeholders['num_features_nonzero']: features[1].shape})
    for ii in range(hier_layer):
        feed_dict.update({placeholders_hc[ii]['features']: features_list[ii]})
        feed_dict.update({placeholders_hc[ii]['support']: support_list[ii]})
        feed_dict.update({placeholders_hc[ii]['lookup']: lookup_list[ii]})
    return feed_dict


def preprocess_multicluster(adj,
                            parts,
                            features,
                            y_train,
                            train_mask,
                            num_clusters,
                            block_size,
                            diag_lambda=-1):
    """Generate the batch for multiple clusters."""

    features_batches = []
    support_batches = []
    y_train_batches = []
    train_mask_batches = []
    total_nnz = 0
    np.random.shuffle(parts)
    for _, st in enumerate(range(0, num_clusters, block_size)):
        pt = parts[st]
        for pt_idx in range(st + 1, min(st + block_size, num_clusters)):
            pt = np.concatenate((pt, parts[pt_idx]), axis=0)
        features_batches.append(features[pt, :])
        y_train_batches.append(y_train[pt, :])
        support_now = adj[pt, :][:, pt]
        if diag_lambda == -1:
            support_batches.append(sparse_to_tuple(normalize_adj(support_now)))
        else:
            support_batches.append(
                sparse_to_tuple(
                    normalize_adj_diag_enhance(
                        support_now,
                        diag_lambda)))
        total_nnz += support_now.count_nonzero()

        train_pt = []
        for newidx, idx in enumerate(pt):
            if train_mask[idx]:
                train_pt.append(newidx)
        train_mask_batches.append(sample_mask(train_pt, len(pt)))
    return (features_batches, support_batches, y_train_batches,
            train_mask_batches)

def preprocess_multicluster_w(adj,
                            parts,
                            features,
                            y_train,
                            train_mask,
                            num_clusters,
                            block_size,
                              nodes_weight,
                            diag_lambda=-1,
                              sym_prop=True):
    """Generate the batch for multiple clusters."""

    features_batches = []
    support_batches = []
    y_train_batches = []
    train_mask_batches = []
    nodes_weight_batches=[]
    real_parts =[]
    total_nnz = 0
    np.random.shuffle(parts)
    if sym_prop:
        adj_prop = adj+adj.transpose()
    else:
        adj_prop = adj
    for _, st in enumerate(range(0, num_clusters, block_size)):
        pt = parts[st]
        for pt_idx in range(st + 1, min(st + block_size, num_clusters)):
            pt = np.concatenate((pt, parts[pt_idx]), axis=0)
        real_parts.append(pt)
        features_batches.append(features[pt, :])
        y_train_batches.append(y_train[pt, :])
        nodes_weight_batches.append(nodes_weight[pt])
        support_now = adj_prop[pt, :][:, pt]
        if diag_lambda == -1:
            support_batches.append(sparse_to_tuple(normalize_adj(support_now)))
        else:
            support_batches.append(
                sparse_to_tuple(
                    normalize_adj_diag_enhance(
                        support_now,
                        diag_lambda)))
        total_nnz += support_now.count_nonzero()

        train_pt = []
        for newidx, idx in enumerate(pt):
            if train_mask[idx]:
                train_pt.append(newidx)
        train_mask_batches.append(sample_mask(train_pt, len(pt)))
    return (features_batches, support_batches, y_train_batches,
            train_mask_batches,nodes_weight_batches,real_parts)

def preprocess_multicluster_important(adj,
                            parts,
                            features,
                            y_train,
                            train_mask,
                            num_clusters,
                            block_size,
                              nodes_weight,
                            diag_lambda=-1,
                              sym_prop=True,
                            important_rate = 0.01):
    """Generate the batch for multiple clusters."""

    features_batches = []
    support_batches = []
    y_train_batches = []
    train_mask_batches = []
    nodes_weight_batches=[]
    real_parts =[]
    total_nnz = 0
    np.random.shuffle(parts)
    if sym_prop:
        adj_prop = adj+adj.transpose()
    else:
        adj_prop = adj

    if important_rate>0:
        important_nb = int(important_rate * len(nodes_weight))
        important_idx = np.argpartition(nodes_weight, -important_nb)[-important_nb:]

    for _, st in enumerate(range(0, num_clusters, block_size)):
        pt = parts[st]
        for pt_idx in range(st + 1, min(st + block_size, num_clusters)):
            pt = np.concatenate((pt, parts[pt_idx]), axis=0)
        if important_rate>0:
            tmpt = np.array(list(set(pt).difference(important_idx)))
            #print("Tmpt  ",tmpt.shape," important  ",important_idx.shape)
            pt = np.concatenate((tmpt,important_idx))
        real_parts.append(pt)
        features_batches.append(features[pt, :])
        y_train_batches.append(y_train[pt, :])
        nodes_weight_batches.append(nodes_weight[pt])
        support_now = adj_prop[pt, :][:, pt]
        if diag_lambda == -1:
            support_batches.append(sparse_to_tuple(normalize_adj(support_now)))
        else:
            support_batches.append(
                sparse_to_tuple(
                    normalize_adj_diag_enhance(
                        support_now,
                        diag_lambda)))
        total_nnz += support_now.count_nonzero()

        train_pt = []
        for newidx, idx in enumerate(pt):
            if train_mask[idx]:
                train_pt.append(newidx)
        train_mask_batches.append(sample_mask(train_pt, len(pt)))
    return (features_batches, support_batches, y_train_batches,
            train_mask_batches,nodes_weight_batches,real_parts)


def preprocess(adj,
               features,
               y_train,
               train_mask,
               visible_data,
               num_clusters,
               diag_lambda=-1):
    """Do graph partitioning and preprocessing for SGD training."""

    # Do graph partitioning
    part_adj, parts = partition_graph(adj, visible_data,
                                                      num_clusters)
    if diag_lambda == -1:
        part_adj = normalize_adj(part_adj)
    else:
        part_adj = normalize_adj_diag_enhance(part_adj, diag_lambda)
    parts = [np.array(pt) for pt in parts]

    features_batches = []
    support_batches = []
    y_train_batches = []
    train_mask_batches = []
    total_nnz = 0
    for pt in parts:
        features_batches.append(features[pt, :])
        now_part = part_adj[pt, :][:, pt]
        total_nnz += now_part.count_nonzero()
        support_batches.append(sparse_to_tuple(now_part))
        y_train_batches.append(y_train[pt, :])

        train_pt = []
        for newidx, idx in enumerate(pt):
            if train_mask[idx]:
                train_pt.append(newidx)
        train_mask_batches.append(sample_mask(train_pt, len(pt)))
    return (parts, features_batches, support_batches, y_train_batches,
            train_mask_batches)

def preprocess_w(adj,
               features,
               y_train,
               train_mask,
               visible_data,
               num_clusters,
                 nodes_weight,
                 weight_reg=5,
               diag_lambda=-1,
                 sym_prop = True,
                 sample_rate = -1):
    """Do graph partitioning and preprocessing for SGD training."""

    # Do graph partitioning
    part_adj, parts = partition_graph_weight(adj, visible_data,
                                                      num_clusters,weight_reg,sym_prop=sym_prop,sample_rate=sample_rate)
    if diag_lambda == -1:
        part_adj = normalize_adj(part_adj)
    else:
        part_adj = normalize_adj_diag_enhance(part_adj, diag_lambda)
    parts = [np.array(pt) for pt in parts]

    features_batches = []
    support_batches = []
    y_train_batches = []
    train_mask_batches = []
    nodes_weight_batches = []
    total_nnz = 0
    for pt in parts:
        features_batches.append(features[pt, :])
        now_part = part_adj[pt, :][:, pt]
        total_nnz += now_part.count_nonzero()
        support_batches.append(sparse_to_tuple(now_part))
        y_train_batches.append(y_train[pt, :])
        nodes_weight_batches.append(nodes_weight[pt])

        train_pt = []
        for newidx, idx in enumerate(pt):
            if train_mask[idx]:
                train_pt.append(newidx)
        train_mask_batches.append(sample_mask(train_pt, len(pt)))
    return (parts, features_batches, support_batches, y_train_batches,
            train_mask_batches,nodes_weight_batches)

def preprocess_w_important(adj,
               features,
               y_train,
               train_mask,
               visible_data,
               num_clusters,
                 nodes_weight,
                 weight_reg=5,
               diag_lambda=-1,
                 sym_prop = True,
                 sample_rate = -1,
                important_rate=0.01):
    """Do graph partitioning and preprocessing for SGD training."""

    # Do graph partitioning
    part_adj, parts = partition_graph_weight(adj, visible_data,
                                                      num_clusters,weight_reg,sym_prop=sym_prop,sample_rate=sample_rate)
    if diag_lambda == -1:
        part_adj = normalize_adj(part_adj)
    else:
        part_adj = normalize_adj_diag_enhance(part_adj, diag_lambda)

    if important_rate>0:
        important_nb = int(important_rate * len(nodes_weight))
        important_idx = np.argpartition(nodes_weight, -important_nb)[-important_nb:]

        parts = [np.concatenate((np.array(list(set(pt).difference(important_idx))),important_idx)) for pt in parts]
    else:
        parts = [np.array(pt) for pt in parts]

    features_batches = []
    support_batches = []
    y_train_batches = []
    train_mask_batches = []
    nodes_weight_batches = []
    total_nnz = 0
    for pt in parts:
        features_batches.append(features[pt, :])
        now_part = part_adj[pt, :][:, pt]
        total_nnz += now_part.count_nonzero()
        support_batches.append(sparse_to_tuple(now_part))
        y_train_batches.append(y_train[pt, :])
        nodes_weight_batches.append(nodes_weight[pt])

        train_pt = []
        for newidx, idx in enumerate(pt):
            if train_mask[idx]:
                train_pt.append(newidx)
        train_mask_batches.append(sample_mask(train_pt, len(pt)))
    return (parts, features_batches, support_batches, y_train_batches,
            train_mask_batches,nodes_weight_batches)

def load_graphsage_data(dataset_path, dataset_str, normalize=True):
    """Load GraphSAGE data."""
    start_time = time.time()

    graph_json = json.load(
        gfile.Open('{}/{}/{}-G.json'.format(dataset_path, dataset_str,
                                            dataset_str)))
    graph_nx = json_graph.node_link_graph(graph_json)

    id_map = json.load(
        gfile.Open('{}/{}/{}-id_map.json'.format(dataset_path, dataset_str,
                                                 dataset_str)))
    is_digit = list(id_map.keys())[0].isdigit()
    id_map = {(int(k) if is_digit else k): int(v) for k, v in id_map.items()}
    class_map = json.load(
        gfile.Open('{}/{}/{}-class_map.json'.format(dataset_path, dataset_str,
                                                    dataset_str)))

    is_instance = isinstance(list(class_map.values())[0], list)
    class_map = {(int(k) if is_digit else k): (v if is_instance else int(v))
                 for k, v in class_map.items()}

    broken_count = 0
    to_remove = []
    for node in graph_nx.nodes():
        if node not in id_map:
            to_remove.append(node)
            broken_count += 1
    for node in to_remove:
        graph_nx.remove_node(node)
    tf.logging.info(
        'Removed %d nodes that lacked proper annotations due to networkx versioning issues',
        broken_count)

    feats = np.load(
        gfile.Open(
            '{}/{}/{}-feats.npy'.format(dataset_path, dataset_str, dataset_str),
            'rb')).astype(np.float32)

    tf.logging.info('Loaded data (%f seconds).. now preprocessing..',
                    time.time() - start_time)
    start_time = time.time()

    edges = []
    for edge in graph_nx.edges():
        if edge[0] in id_map and edge[1] in id_map:
            edges.append((id_map[edge[0]], id_map[edge[1]]))
    num_data = len(id_map)

    val_data = np.array(
        [id_map[n] for n in graph_nx.nodes() if graph_nx.nodes[n]['val']],
        dtype=np.int32)
    test_data = np.array(
        [id_map[n] for n in graph_nx.nodes() if graph_nx.nodes[n]['test']],
        dtype=np.int32)
    is_train = np.ones((num_data), dtype=np.bool)
    is_train[val_data] = False
    is_train[test_data] = False
    train_data = np.array([n for n in range(num_data) if is_train[n]],
                          dtype=np.int32)

    train_edges = [
        (e[0], e[1]) for e in edges if is_train[e[0]] and is_train[e[1]]
    ]
    edges = np.array(edges, dtype=np.int32)
    train_edges = np.array(train_edges, dtype=np.int32)

    # Process labels
    if isinstance(list(class_map.values())[0], list):
        num_classes = len(list(class_map.values())[0])
        labels = np.zeros((num_data, num_classes), dtype=np.float32)
        for k in class_map.keys():
            labels[id_map[k], :] = np.array(class_map[k])
    else:
        num_classes = len(set(class_map.values()))
        labels = np.zeros((num_data, num_classes), dtype=np.float32)
        for k in class_map.keys():
            labels[id_map[k], class_map[k]] = 1

    if normalize:
        train_ids = np.array([
            id_map[n]
            for n in graph_nx.nodes()
            if not graph_nx.nodes[n]['val'] and not graph_nx.nodes[n]['test']
        ])
        train_feats = feats[train_ids]
        scaler = sklearn.preprocessing.StandardScaler()
        scaler.fit(train_feats)
        feats = scaler.transform(feats)

    def _construct_adj(edges):
        adj = sp.csr_matrix((np.ones(
            (edges.shape[0]), dtype=np.float32), (edges[:, 0], edges[:, 1])),
            shape=(num_data, num_data))
        adj += adj.transpose()
        return adj

    train_adj = _construct_adj(train_edges)
    full_adj = _construct_adj(edges)

    train_feats = feats[train_data]
    test_feats = feats

    tf.logging.info('Data loaded, %f seconds.', time.time() - start_time)
    return num_data, train_adj, full_adj, feats, train_feats, test_feats, labels, train_data, val_data, test_data


def process_part_to_hier(features_batch_list,
                         support_batch_list,
                         lookup_list,
                         num_hier_layer,
                         hyper_node_size=3):
    if num_hier_layer == 0:
        del features_batch_list[0]
        del support_batch_list[0]
        features_batch_list.reverse()
        support_batch_list.reverse()
        lookup_list.reverse()
        return features_batch_list, support_batch_list, lookup_list
    #print("In process hier layer: ",num_hier_layer)
    num_nodes = features_batch_list[-1].shape[0]
    part_coo = tuple_to_sparse(support_batch_list[-1])
    part_adj = part_coo.tolil()
    #part_adj = sp.coo_matrix(support_batch_list[-1]).tolil()
    edge_list = [part_adj[i].rows[0] for i in range(num_nodes)]
    num_nodes_sub = int(num_nodes / hyper_node_size)
    assert num_nodes / \
        num_nodes_sub > 2.0, "size of hyper_node is too small to do graph partition"
    _, parts = metis.part_graph(edge_list, num_nodes_sub, seed=1)
    parts = np.array(parts)
    gps, parts = np.unique(parts, return_inverse=True)
    num_nodes_sub_real = len(gps)
    assert num_nodes_sub_real > 1, "Only one group"

    part_row = []
    part_col = []
    part_value = []
    groups = [[] for _ in range(num_nodes_sub_real)]
    groups_feat = []
    for i in range(num_nodes):
        rows = part_adj[i].rows[0]
        values = part_adj[i].data[0]
        for j in range(len(rows)):
            part_row.append(parts[i])
            part_col.append(parts[rows[j]])
            part_value.append(values[j])
        groups[parts[i]].append(i)
    for i in range(num_nodes_sub_real):
        feats = features_batch_list[-1][groups[i]]
        feat_avg = np.average(feats, axis=0)
        feat_max = np.max(feats, axis=0)
        feat_sum = np.sum(feats, axis=0)
        groups_feat.append(
            np.concatenate(
                (feat_avg, feat_max, feat_sum), axis=1))
        # groups_feat.append(np.concatenate((feat_avg,feat_max)))
    groups_feat = np.array(groups_feat).reshape((num_nodes_sub_real, -1))

    part_adj_n = sp.coo_matrix(
        (part_value, (part_row, part_col)), shape=(
            num_nodes_sub_real, num_nodes_sub_real))

    features_batch_list.append(groups_feat)
    support_batch_list.append(sparse_to_tuple(sym_normalize_adj(part_adj_n)))
    lookup_list.append(parts)

    return process_part_to_hier(features_batch_list,
                                support_batch_list,
                                lookup_list,
                                num_hier_layer=num_hier_layer - 1,
                                hyper_node_size=hyper_node_size)

def update_adj_edge_weight(adj,nodes_weight,feature,outputs,num_clusters,threshold=-1):

    #feature = feature / np.sqrt(np.sum(np.square(feature),axis=-1))[:,np.newaxis]
    #outputs = outputs / np.sqrt(np.sum(np.square(outputs),axis=-1))[:,np.newaxis]
    adj_coo = adj.tocoo()

    ##
    # visible_dict = dict(zip(visible_data,list(range(len(visible_data)))))
    # row =  list(map(visible_dict.__getitem__,adj_coo.row))
    # col =  list(map(visible_dict.__getitem__,adj_coo.col))
    row = adj_coo.row
    col = adj_coo.col

    nodes_weight = np.array(nodes_weight).flatten()
    part_length = len(row)//num_clusters+1
    p_i2j=[]
    for kk in range(num_clusters):
        row_k = row[kk*part_length:(kk+1)*part_length]
        col_k = col[kk*part_length:(kk+1)*part_length]
        w_i = nodes_weight[row_k]
        w_j = nodes_weight[col_k]
        x_i = feature[row_k]
        x_j = feature[col_k]
        h_i = outputs[row_k]
        h_j = outputs[col_k]
        simijx = np.sum(x_i*x_j,axis=-1)
        simijh = np.sum(h_i*h_j,axis=-1)

        # Paper version
        p_i2j_k = 1 - (1 - w_i) * w_j * (1 - simijx) * simijh - w_i * (1 - w_j) * simijx * (1 - simijh)

        # Stable to 0.9947
        # p_i2j_k = 1 - (1-simijh*(1-w_j))*(1-w_j*simijx*(1-simijh))

        #p_i2j_k = (1-w_j)*(1-w_j*simijx*(1-simijh)) + w_i*(1-w_j)*(1-simijx*(1-simijh))
        p_i2j.append(p_i2j_k)
    p_i2j = np.concatenate(p_i2j)
    p_i2j = (p_i2j-np.min(p_i2j))/(np.max(p_i2j)-np.min(p_i2j))

    ## (col,row) or (row,col) ???
    return sp.csr_matrix((p_i2j,(col,row)),shape=adj_coo.shape)
    #return sp.csr_matrix((p_i2j,(row,col)),shape=adj_coo.shape)

def load_saint_data(prefix, normalize=True):
    adj_full = sp.load_npz('./{}/adj_full.npz'.format(prefix)).astype(np.bool)
    adj_train = sp.load_npz('./{}/adj_train.npz'.format(prefix)).astype(np.bool)
    role = json.load(open('./{}/role.json'.format(prefix)))
    feats = np.load('./{}/feats.npy'.format(prefix))
    class_map = json.load(open('./{}/class_map.json'.format(prefix)))
    class_map = {int(k):v for k,v in class_map.items()}
    assert len(class_map) == feats.shape[0]
    # ---- normalize feats ----
    train_nodes = np.array(list(set(adj_train.nonzero()[0])))
    train_feats = feats[train_nodes]
    scaler = sklearn.preprocessing.StandardScaler()
    scaler.fit(train_feats)
    feats = scaler.transform(feats)
    # -------------------------
    num_vertices = adj_full.shape[0]
    if isinstance(list(class_map.values())[0], list):
        num_classes = len(list(class_map.values())[0])
        class_arr = np.zeros((num_vertices, num_classes))
        for k, v in class_map.items():
            class_arr[k] = v
    else:
        num_classes = max(class_map.values()) - min(class_map.values()) + 1
        class_arr = np.zeros((num_vertices, num_classes))
        offset = min(class_map.values())
        for k, v in class_map.items():
            class_arr[k][v - offset] = 1


    num_data = train_nodes.shape[0]
    return adj_full, adj_train, feats, class_arr, role
    #return num_data, train_adj, full_adj, feats, train_feats, test_feats, labels, train_data, val_data, test_data

def load_gcn_data(dataset_path, dataset_str):
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
    names = ['x', 'y', 'tx', 'ty', 'allx', 'ally', 'graph']
    objects = []
    for i in range(len(names)):
        with open("{}/{}/ind.{}.{}".format(dataset_path,dataset_str,dataset_str, names[i]), 'rb') as f:
            if sys.version_info > (3, 0):
                objects.append(pkl.load(f, encoding='latin1'))
            else:
                objects.append(pkl.load(f))

    x, y, tx, ty, allx, ally, graph = tuple(objects)
    test_idx_reorder = parse_index_file("{}/{}/ind.{}.test.index".format(dataset_path,dataset_str,dataset_str))
    test_idx_range = np.sort(test_idx_reorder)

    if dataset_str == 'citeseer':
        # Fix citeseer dataset (there are some isolated nodes in the graph)
        # Find isolated nodes, add them as zero-vecs into the right position
        test_idx_range_full = range(min(test_idx_reorder), max(test_idx_reorder)+1)
        tx_extended = sp.lil_matrix((len(test_idx_range_full), x.shape[1]))
        tx_extended[test_idx_range-min(test_idx_range), :] = tx
        tx = tx_extended
        ty_extended = np.zeros((len(test_idx_range_full), y.shape[1]))
        ty_extended[test_idx_range-min(test_idx_range), :] = ty
        ty = ty_extended

    features = sp.vstack((allx, tx)).tolil()
    features[test_idx_reorder, :] = features[test_idx_range, :]
    adj = nx.adjacency_matrix(nx.from_dict_of_lists(graph))
    #laplacian = nx.laplacian_matrix(nx.from_dict_of_lists(graph))

    labels = np.vstack((ally, ty))
    labels[test_idx_reorder, :] = labels[test_idx_range, :]

    # GCN split
    # idx_test = test_idx_range.tolist()
    # idx_train = range(len(y))
    # idx_val = range(len(y), len(y)+500)

    # FastGCN Split
    idx_test = test_idx_range.tolist()
    idx_train = range(len(ally) - 500)
    idx_val = range(len(ally) - 500, len(ally))

    # train_mask = sample_mask(idx_train, labels.shape[0])
    # val_mask = sample_mask(idx_val, labels.shape[0])
    # test_mask = sample_mask(idx_test, labels.shape[0])

    # y_train = np.zeros(labels.shape)
    # y_val = np.zeros(labels.shape)
    # y_test = np.zeros(labels.shape)
    # y_train[train_mask, :] = labels[train_mask, :]
    # y_val[val_mask, :] = labels[val_mask, :]
    # y_test[test_mask, :] = labels[test_mask, :]
    #return adj, features, y_train, y_val, y_test, train_mask, val_mask, test_mask, laplacian


    num_data = labels.shape[0]
    full_adj = adj
    train_adj = full_adj[idx_train,:][:,idx_train].tocoo()
    train_adj = sp.csr_matrix((train_adj.data,(train_adj.row,train_adj.col)),shape=(num_data,num_data))

    feats = test_feats= features.todense()
    train_feats = feats[idx_train]
    train_data = idx_train
    val_data = idx_val
    test_data = idx_test
    return num_data, train_adj, full_adj, feats, train_feats, test_feats, labels, train_data, val_data, test_data

def weight_pagerank(original_weight,adj,round,alpha):
    weight0 = np.reshape(original_weight,[-1,1])
    new_weight = weight0.copy()
    for _ in range(round):
        new_weight = (1-alpha)*adj*new_weight + alpha*weight0
    return np.array(new_weight).flatten() / np.average(new_weight)

def print_configuration_op(FLAGS):
    tf.logging.info('My Configurations:')
    #pdb.set_trace()
    for name, value in FLAGS.__flags.items():
        value=value.value
        if type(value) == float:
            tf.logging.info(' %s:\t %f'%(name, value))
        elif type(value) == int:
            tf.logging.info(' %s:\t %d'%(name, value))
        elif type(value) == str:
            tf.logging.info(' %s:\t %s'%(name, value))
        elif type(value) == bool:
            tf.logging.info(' %s:\t %s'%(name, value))
        else:
            tf.logging.info('%s:\t %s' % (name, value))
    #for k, v in sorted(FLAGS.__dict__.items()):
        #print(f'{k}={v}\n')
    tf.logging.info('End of configuration')

def reindex_parts_important(parts, re_list, total_len, important_rate):
    if(important_rate<=0):
        p = np.concatenate(parts)
        ll_con = np.concatenate(re_list, axis=0)
        ll_new = np.ones(shape=(total_len, ll_con.shape[1]))
        for ii in range(len(p)):
            ll_new[p[ii], :] = ll_con[ii, :]
        new_list = np.array(ll_new)
        return new_list
    else:
        important_nb = int(important_rate * total_len)

        tf.logging.info("Reindexing Concat")
        p = np.concatenate([x[:-important_nb] for x in parts])
        p = np.concatenate((p, parts[0][-important_nb:]))

        ll_con = np.concatenate([x[:-important_nb] for x in re_list], axis=0)  # corresponding content
        ll_imp = np.average([x[-important_nb:] for x in re_list],axis=0)
        ll_con = np.concatenate((ll_con,ll_imp),axis=0)

        tf.logging.info("Constructing list")
        ll_new= np.ones(shape=(total_len,ll_con.shape[1]))
        for ii in range(len(p)):
            ll_new[p[ii],:] = ll_con[ii,:]
        new_list = np.array(ll_new)
        return new_list

        # p = np.concatenate(parts)  # indexes,not unique
        # ll_con = np.concatenate(re_list, axis=0)  # corresponding content
        #
        #
        # tf.logging.info("Reindexing Append")
        # ll_list = [[]] * total_len
        # for ii in range(len(p)):
        #     ll_list[p[ii]].append(ll_con[ii, :])
        #
        # tf.logging.info("Reindexing Average")
        # ll_new= np.ones(shape=(total_len,ll_con.shape[1]))
        # for ii in range(total_len):
        #     if(ii%100 == 0):
        #         tf.logging.info("averaging %d"%(ii))
        #     if(len(ll_list[ii])>1):
        #         ll_new[ii,:] = np.average(ll_list[ii], axis=0)
        #     else:
        #         ll_new[ii,:] = ll_list[ii][0]
        #
        # return ll_new

def load_adclick_data_pai(FLAGS):
    tables = (FLAGS.tables).split(",")
    reader = tf.python_io.TableReader(tables[0],selected_cols = "")
    total_num = reader.get_row_count()
    arr_all = np.array(reader.read(total_num))
    reader.close()

def load_adclick_data_server(dataset_path,dataset_str):
    pass




