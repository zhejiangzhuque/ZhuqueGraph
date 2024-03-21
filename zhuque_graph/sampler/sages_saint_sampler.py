import copy
import os.path as osp

import torch
import scipy.sparse as sp
from tqdm import tqdm
from torch_sparse import SparseTensor
from graph_sample.fastgae_sampling import (get_distribution, node_sampling)

class GraphSAINTSampler(object):
    r"""
    Args:
        data (torch_geometric.data.Data): The graph data object.
        batch_size (int): The approximate number of samples per batch to load.
        num_steps (int, optional): The number of iterations.
            (default: :obj:`1`)
        sample_coverage (int): How many samples per node should be used to
            compute normalization statistics. (default: :obj:`50`)
        save_dir (string, optional): If set, will save normalization
            statistics to the :obj:`save_dir` directory for faster re-use.
            (default: :obj:`None`)

        log (bool, optional): If set to :obj:`False`, will not log any
            progress. (default: :obj:`True`)
    """
    def __init__(self, data, batch_size, num_steps=1, sample_coverage=50,
                 save_dir=None, log=True, reload=False, cal_norm=True):
        """
        :param data:            pyG格式的data
        :param batch_size:      每个子图中有，多少条边，或者节点.
        :param num_steps:       重复几次等等
        :param sample_coverage: 计算采样概率， 每个节点重复几次
        :param save_dir:
        :param log:
        """
        # 保证必须有边，且未算过 node_norm 和 edge_norm 这两个参数
        assert data.edge_index is not None
        assert 'node_norm' not in data
        assert 'edge_norm' not in data
        # N是节点个数，E是边个数
        self.N = N = data.num_nodes
        self.E = data.num_edges
        # 构建稀疏矩阵，用的是SparseTensor，矩阵大小为（N, N）
        self.adj = SparseTensor(row=data.edge_index[0], col=data.edge_index[1],
                                value=data.edge_attr, sparse_sizes=(N, N))
        # 浅层拷贝data，同时把结构信息，边，边值清空
        self.data = copy.copy(data)
        self.data.edge_index = None
        self.data.edge_attr = None

        self.batch_size = batch_size
        self.num_steps = num_steps
        self.sample_coverage = sample_coverage
        self.log = log
        self.__count__ = 0
        self.cal_norm = cal_norm

        if cal_norm:
            path = osp.join(save_dir or '', self.__filename__)
            if save_dir is not None and osp.exists(path) and reload==False:  # pragma: no cover
                self.node_norm, self.edge_norm = torch.load(path)
            else:
                self.node_norm, self.edge_norm, self.edge_sample_norm = self.__compute_norm__()
                if save_dir is not None:  # pragma: no cover
                    torch.save((self.node_norm, self.edge_norm), path)
        else:
            self.node_norm, self.edge_norm = None, None

    @property
    def __filename__(self):
        return f'{self.__class__.__name__.lower()}_{self.sample_coverage}.pt'

    # 给定参数num_examples，得到多个子图的节点，存为一个list，每个元素是子图节点列表
    def __sample_nodes__(self, num_examples):
        raise NotImplementedError

    # 解答; num_examples是子图个数，正常flickr中，其为1
    def __sample__(self, num_examples):
        node_samples = self.__sample_nodes__(num_examples)
        #     samples存储每个子图返回的结果，为节点列表，边列表，adj
        # 问题： 节点编号是否从0开始（对子图来讲），edge边中两个节点编号是否从0开始？ adj格式
        # 解答： 节点编号是全局id，edge编号也是全局id，可见compute中计数部分
        samples = []
        for node_idx in node_samples:
            # 这里去重
            node_idx = node_idx.unique()
            adj, edge_idx = self.adj.saint_subgraph(node_idx)
            samples.append((node_idx, edge_idx, adj))
        return samples

    def __compute_norm__(self):
        node_count = torch.zeros(self.N, dtype=torch.float)
        edge_count = torch.zeros(self.E, dtype=torch.float)

        if self.log:  # pragma: no cover
            pbar = tqdm(total=self.N * self.sample_coverage)
            pbar.set_description('Compute GraphSAINT normalization')
        # num_samples含义， total_sampled_nodes含义（所有被sample的节点的个数和）
        num_samples = total_sampled_nodes = 0
        # sample_converage可能表示，预期每个node被sample多少次，如50，需要sample的所有图的node个数为 50 * N
        while total_sampled_nodes < self.N * self.sample_coverage:
            num_sampled_nodes = 0
            samples = self.__sample__(200)
            for node_idx, edge_idx, _ in samples:
                node_count[node_idx] += 1
                edge_count[edge_idx] += 1
                num_sampled_nodes += node_idx.size(0)
            total_sampled_nodes += num_sampled_nodes
            num_samples += 200

            if self.log:  # pragma: no cover
                pbar.update(num_sampled_nodes)

        if self.log:  # pragma: no cover
            pbar.close()

        row, col, _ = self.adj.coo()
        # 采样的所有子图中，v节点出现的次数/所有边的次数
        # 问题，弄清楚node_count是对源节点，还是目标节点，根据论文是源节点
        edge_norm = (node_count[col] / edge_count).clamp_(0, 1e4)
        edge_norm[torch.isnan(edge_norm)] = 0.1
        # 问题：查看node_norm和edge_norm的大小，大于1还是小于1，弄明白含义
        node_count[node_count == 0] = 0.1
        node_norm = num_samples / node_count / self.N
        # 新加入的edge norm
        edge_count[edge_count == 0] = 0.1
        edge_sample_norm = num_samples / edge_count / self.E
        return node_norm, edge_norm, edge_sample_norm

    # 这里可以看出node_idx和edge_idx都是全局的id
    def __get_data_from_sample__(self, sample):
        node_idx, edge_idx, adj = sample

        data = self.data.__class__()
        data.num_nodes = node_idx.size(0)
        row, col, value = adj.coo()
        data.edge_index = torch.stack([row, col], dim=0)
        data.edge_attr = value

        for key, item in self.data:
            if item.size(0) == self.N:
                data[key] = item[node_idx]
            elif item.size(0) == self.E:
                data[key] = item[edge_idx]
            else:
                data[key] = item

        if self.cal_norm:
            data.node_norm = self.node_norm[node_idx]
            data.edge_norm = self.edge_norm[edge_idx]

        return data

    def __next__(self):
        if self.__count__ < len(self):
            self.__count__ += 1
            sample = self.__sample__(1)[0]
            data = self.__get_data_from_sample__(sample)
            return data
        else:
            raise StopIteration

    def __len__(self):
        return self.num_steps

    def __iter__(self):
        self.__count__ = 0
        return self


class GraphSAINTNodeSampler(GraphSAINTSampler):
    r"""The GraphSAINT node sampler class (see
    :class:`torch_geometric.data.GraphSAINTSampler`).
    Args:
        batch_size (int): The number of nodes to sample per batch.
    """
    def __sample_nodes__(self, num_examples):
        # 该sampler下，num_examples表示有多少张子图，batch_size表示每张子图下有多少条边
        edge_sample = torch.randint(0, self.E, (num_examples, self.batch_size),
                                    dtype=torch.long)
        node_sample = self.adj.storage.row()[edge_sample]
        # unbind是tensor到函数，转为tuple，原来node_sample是[1, 18000],现在变为tuple了
        return node_sample.unbind(dim=0)


class GraphSAINTEdgeSampler(GraphSAINTSampler):
    r"""The GraphSAINT edge sampler class (see
    :class:`torch_geometric.data.GraphSAINTSampler`).
    Args:
        batch_size (int): The number of edges to sample per batch.
    """
    def __sample_nodes__(self, num_examples):
        # This function corresponds to the `Edge2` sampler in the official
        # code repository that weights all edges as equally important.
        # This is the default configuration in the GraphSAINT implementation.
        edge_sample = torch.randint(0, self.E, (num_examples, self.batch_size),
                                    dtype=torch.long)

        source_node_sample = self.adj.storage.row()[edge_sample]
        target_node_sample = self.adj.storage.col()[edge_sample]

        node_sample = torch.cat([source_node_sample, target_node_sample], -1)
        # node_sample.unbind是去掉某一维度的操作
        return node_sample.unbind(dim=0)


class GraphSAINTRandomWalkSampler(GraphSAINTSampler):
    r"""The GraphSAINT random walk sampler class (see
    :class:`torch_geometric.data.GraphSAINTSampler`).
    Args:
        batch_size (int): The number of walks to sample per batch.
        walk_length (int): The length of each random walk.
    """
    def __init__(self, data, batch_size, walk_length, num_steps=1,
                 sample_coverage=50, save_dir=None, log=True, reload=False, cal_norm=True):
        self.walk_length = walk_length
        super(GraphSAINTRandomWalkSampler,
              self).__init__(data, batch_size, num_steps, sample_coverage,
                             save_dir, log, reload, cal_norm=cal_norm)

    @property
    def __filename__(self):
        return (f'{self.__class__.__name__.lower()}_{self.walk_length}_'
                f'{self.sample_coverage}.pt')

    def __sample_nodes__(self, num_examples):
        start = torch.randint(0, self.N, (num_examples, self.batch_size),
                              dtype=torch.long)
        node_sample = self.adj.random_walk(start.flatten(), self.walk_length)
        node_sample = node_sample.view(
            num_examples, self.batch_size * (self.walk_length + 1))
        return node_sample.unbind(dim=0)

class FastGAENodeSampler(GraphSAINTSampler):
    def __init__(self, data, measure='degree', alpha=2, batch_size=5000,
                 num_steps=3, sample_coverage=20, save_dir=None, log=True, reload=False, cal_norm=True):
        """
        :param data:        data of pyG
        :param alpha:       alpha scalar hyperparameter for degree and core sampling
        :param batch_size:  node numbers of each subgraph(batch)
        :param num_steps:
        :param sample_coverage:
        :param save_dir:
        :param log:
        """
        self.measure = measure
        self.alpha = alpha
        assert data.edge_index is not None
        # N是节点个数，E是边个数
        N = data.num_nodes
        # SparseTensor(row=data.edge_index[0], col=data.edge_index[1],
        #              value=data.edge_attr, sparse_sizes=(N, N))
        self.sp_adj = sp.coo_matrix((data.edge_attr, (data.edge_index[0], data.edge_index[1])), shape=(N, N))
        self.node_distribution = get_distribution(self.measure, self.alpha, self.sp_adj)
        assert len(self.node_distribution)==N, 'node dirstirbution is illegal'
        super(FastGAENodeSampler, self).__init__(data, batch_size, num_steps, sample_coverage, save_dir, log, reload, cal_norm=cal_norm)

    @property
    def __filename__(self):
        return (f'{self.__class__.__name__.lower()}_{self.measure}_{self.alpha}'
                f'{self.sample_coverage}.pt')

    def __sample_nodes__(self, num_examples):
        node_sample = []
        for _ in range(num_examples):
            nodes = node_sampling(adj=self.sp_adj, distribution=self.node_distribution,
                                  nb_node_samples=self.batch_size, replace=False)
            node_sample.append(torch.from_numpy(nodes))
        return node_sample

