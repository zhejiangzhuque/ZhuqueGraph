import numpy as np
from collections import Counter
import torch
import torch.utils
from torch_geometric.data import Data


class ExprDataset(torch.utils.data.Dataset):
    def __init__(self, Expr, edge, y):
        super(ExprDataset, self).__init__()
        self.Expr = Expr  # (1000, 23459)
        self.edge = edge
        self.y = y  # (1000,)
        self.num_samples = len(self.y)  # 1000
        self.sample_idx = np.arange(self.num_samples)  # (1000,)
        # the dimension of features
        if len(self.Expr.shape) == 2: self.num_expr_feature = 1
        else: self.num_expr_feature = self.Expr.shape[2]

    def __getitem__(self, idx):
        if not isinstance(idx, int):
            raise IndexError(
                'Only integers are valid '
                'indices (got {}).'.format(type(idx).__name__))
        idx = self.sample_idx[idx]
        data = Data()
        data['x'] = torch.tensor(self.Expr[idx, :].reshape([-1, self.num_expr_feature])).float()
        data['y'] = torch.tensor(self.y[idx].reshape([1, 1])).long()
        data['edge_index'] = torch.tensor(self.edge)
        return data

    def __len__(self):
        return len(self.sample_idx)

    def split(self, idx):
        return ExprDataset(self.Expr[idx, :], self.edge, self.y[idx])

    def data_aug_minor_types(self, thr, random_seed):
        np.random.seed(random_seed)
        ''' data augmentation is performed on small classes '''
        cat_count = Counter(self.y)
        max_cat_count = max(cat_count.values())
        for label in cat_count.keys():
            # if minor cell types
            # if max_cat_count / np.sum(self.y == label) > dup_odds:
            #     add_size = int(max_cat_count / dup_odds) - np.sum(self.y == label)
            #     print('For minor cell type {}, add {} more samples'.format(label, add_size))
            #     add_idx = np.random.choice(np.where(self.y == label)[0], size=add_size, replace=True)
            #     # update
            #     self.sample_idx = np.concatenate((self.sample_idx, add_idx))

            if cat_count[label]/max_cat_count < thr:
                add_size = int(max_cat_count * thr) - np.sum(self.y == label)
                print('For minor cell type {}, add {} more samples'.format(label, add_size))
                add_idx = np.random.choice(np.where(self.y == label)[0], size=add_size, replace=True)
                # update
                self.sample_idx = np.concatenate((self.sample_idx, add_idx))
        np.random.shuffle(self.sample_idx)
        print('org/updated #cells:', self.num_samples, len(self.sample_idx))
        print('org/updated # of each cell types', Counter(self.y[self.sample_idx]))

    def cal_type_weight(self):
        alpha = list(Counter(self.y).values())
        alpha /= np.max(alpha)
        alpha = np.clip(alpha, 0.02, 1)
        alpha /= np.sum(alpha)
        return torch.tensor(alpha).float()


class DataLoader(torch.utils.data.DataLoader):
    def __init__(self, dataset, batch_size=64, shuffle=False, **kwargs):
        if 'collate_fn' not in kwargs.keys():
            raise
        super(DataLoader, self).__init__(dataset, batch_size, shuffle, **kwargs)


def collate_func(batch):
    data0 = batch[0]
    if not isinstance(data0, Data):
        raise
    tmp_data = Data()
    tmp_data['x'] = torch.stack([i['x'] for i in batch], dim=1) #torch.Size([23459, 64, 1])
    tmp_data['y'] = torch.cat([i['y'] for i in batch]) #torch.Size([64, 1])
    tmp_data['edge_index'] = data0.edge_index #torch.Size([2, 109914])
    tmp_data['batch'] = torch.zeros_like(tmp_data['y']) #torch.Size([64, 1])
    return tmp_data
