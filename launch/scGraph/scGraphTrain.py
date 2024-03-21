import os
import time
import argparse
import numpy as np
import torch
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.metrics import precision_score,f1_score
from zhuque_graph.utils.scGraph_scheduler  import (CosineAnnealingWarmRestarts_Decay)
from dataset import (collate_func, DataLoader, ExprDataset)
from zhuque_graph.nn.pytorch.model.scGraphmodel import (scGraph)
from torch_geometric.utils import remove_self_loops


def train(model, optimizer, loader, epoch, scheduler_type=None, verbose=False):
    model.train()
    iters = len(loader)
    loss_all = 0
    # create empty tensor
    embedding = torch.tensor([]).to(device)
    true_label = torch.tensor([], dtype=torch.long).to(device)

    for idx,data in enumerate(loader):
        data = data.to(device)
        optimizer.zero_grad() #Sets the gradients of all optimized torch.Tensors to zero
        output, tmp_embedding = model(data)
        loss = loss_fn(output, data.y.reshape(-1))
        # L1 penalization of edge importance
        if model.edge_weight is not None:
            l2_loss = 0.1 * torch.mean(model.edge_weight ** 2)
            loss += l2_loss
        loss.backward()
        loss_all += loss.item()
        optimizer.step() #parameter update

        if scheduler_type == 'CosineAnnealingWarmRestarts_Decay':
            scheduler.step((epoch -1) + idx/iters) # let "epoch" begin from 0
            lr = optimizer.param_groups[0]['lr']
            step = (epoch-1)*iters+1+idx
            # viz.line(X=[step], Y=[lr], win=viz_lr, update='append')

        # combine embedding and true_label
        embedding = torch.cat((embedding, tmp_embedding), 0) #0 indicate concat vertical
        tmp_true_label = data.y
        true_label = torch.cat((true_label, tmp_true_label), 0)
    return loss_all/iters, embedding, true_label


def test(model, loader, predicts=False):
    model.eval()
    # create empty tensor
    pred_label = torch.tensor([], dtype=torch.long).to(device)
    true_label = torch.tensor([], dtype=torch.long).to(device)

    for data in loader:
        data = data.to(device)
        output, tmp_embedding = model(data)
        tmp_pred_label = output.max(dim=1)[1]   # .cpu().data.numpy()
        tmp_true_label = data.y
        # combine embedding and true_label
        pred_label = torch.cat((pred_label, tmp_pred_label), 0)
        true_label = torch.cat((true_label, tmp_true_label), 0)

    pred_label = pred_label.cpu().detach().numpy()
    true_label = true_label.cpu().detach().numpy()
    acc = precision_score(true_label, pred_label, average='macro')
    f1 = f1_score(true_label, pred_label, average='macro')
    if predicts:
        return acc, f1, pred_label, true_label
    else:
        return acc, f1



if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=str, default='./results/models')
    parser.add_argument('--init_lr', type=float, default= 0.01)
    parser.add_argument('--min_lr', type=float, default= 0.00001)
    parser.add_argument('--max_epoch', type=int, default= 16)
    parser.add_argument('--extend_epoch', type=int, default= 50)
    parser.add_argument('--batch_size', type=int, default= 64)
    parser.add_argument('--weight_decay', type=float, default= 1e-4)
    parser.add_argument('--dropout_ratio', type=float, default= 0.1)
    parser.add_argument('--num_workers', type=int, default= 0) # if debug too slow, change num_workers to 0
    args = parser.parse_args()
    device = torch.device('cuda:1') if torch.cuda.is_available() else torch.device('cpu')
    os.makedirs(args.outdir, exist_ok=True)

    # load data
    data = np.load(f'/INPUT/datasets/scgraph/dataset.npz', allow_pickle=True)
    logExpr = data['logExpr']  #cell x gene (1000, 23459)
    edge = data['edge_index']  #(2, 109916)
    # remove self loop
    edge = remove_self_loops(torch.tensor(edge))[0].numpy() #(2, 109914)
    label_cat = data['label_cat'] #5
    label = data['label'] #1000

    # statistic
    num_class = len(label_cat)
    num_gene = logExpr.shape[1]
    num_edge = edge.shape[1]

    # StratifiedKFold: preserve relative class frequencies in each train and validation fold
    kf = StratifiedKFold(n_splits=5, shuffle=True)
    for tr, ts in kf.split(X=label, y=label):
        train_idx = tr
        test_idx = ts

    # dataloader
    dataset = ExprDataset(Expr=logExpr,edge=edge,y=label)
    train_dataset = dataset.split(torch.tensor(train_idx))
    test_dataset = dataset.split(torch.tensor(test_idx))

    ''' Unbalance data strategy 1: data augmentation '''
    # for minor cell types (<threshold)
    train_dataset.data_aug_minor_types(thr=0.02, random_seed=2240)

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=True, collate_fn=collate_func, drop_last=True)
    test_loader = DataLoader(test_dataset, batch_size=1, num_workers=args.num_workers, collate_fn=collate_func)
    iters = len(train_loader)

    ''' Unbalance data strategy 2: weighted cross-entropy '''
    # assign minor class larger weights
    alpha = train_dataset.cal_type_weight()
    loss_fn = torch.nn.CrossEntropyLoss(weight = alpha)

    model = scGraph(in_channel= dataset.num_expr_feature,  #1
                    mid_channel=8,
                    out_channel=num_class,  # 5
                    num_nodes=num_gene,  #gene as nodes 23459
                    num_edge=num_edge,  #109914
                    dropout_ratio=args.dropout_ratio, #0.1
                    global_conv1_dim=12,
                    global_conv2_dim=4,
                    FC1_dim=256,
                    FC2_dim=64,)

    model.to(device)
    loss_fn = loss_fn.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.init_lr, weight_decay=args.weight_decay)
    scheduler =CosineAnnealingWarmRestarts_Decay(optimizer, T_0=2, T_mult=2, eta_min=args.min_lr, lr_max_decay=0.9)

    # visualize live data
    # viz = visdom.Visdom(env='scGraph', port=7788)
    # x, y = 0, 0
    # viz_loss = viz.line(X=[0], Y=[0], opts=dict(title='Train Loss', xlabel='Epoch'))
    # viz_lr = viz.line(X=[0], Y=[0], opts=dict(title='Learning Rate', xlabel='Step'))
    # viz_f1 = viz.line(X=np.array([[0, np.nan]]), Y=np.array([[0, np.nan]]), opts=dict(title='F1', xlabel='Epoch', legend=['Train', 'Valid']))
    # viz_acc = viz.line(X=np.array([[0, np.nan]]), Y=np.array([[0, np.nan]]), opts=dict(title='Accuracy', xlabel='Epoch', legend=['Train', 'Valid']))
    # nan_embed = np.zeros([len(label_cat), 2])
    # nan_embed[:] = np.nan
    # viz_embed = viz.scatter(X=nan_embed,
    #                         Y=np.arange(1, len(label_cat) + 1),  # Must include all label index
    #                         opts=dict(title='Train Embedding', markersize=3, legend=list(label_cat)))

    for epoch in range(1, args.max_epoch):
        train_loss, train_embedding, train_true_label = train(model, optimizer, train_loader, epoch,
                                                              scheduler_type = 'CosineAnnealingWarmRestarts_Decay')
        train_acc, train_f1 = test(model, train_loader)

        # visualization
        # viz.line(X=[epoch], Y=[train_loss], win=viz_loss, update='append')
        # viz.line(X=np.array([[epoch, np.nan]]), Y=np.array([[train_f1, np.nan]]), win=viz_f1, update='append')
        # viz.line(X=np.array([[epoch, np.nan]]), Y=np.array([[train_acc, np.nan]]), win=viz_acc, update='append')
        # 2D embedding (tsne)
        # train_embedding = train_embedding.cpu().detach().numpy()
        # ts = TSNE(n_components=2)
        # train_embedding = ts.fit_transform(train_embedding)
        # viz.scatter(X=train_embedding, Y=train_true_label+1, win=viz_embed, update='replace')

        lr = optimizer.param_groups[0]['lr']
        print('[Epoch {}] lr:{}, loss:{}, T-acc:{}, T-f1:{}'.format(
            epoch, round(lr,6), round(train_loss,4), round(train_acc,4), round(train_f1,4)))


    # stage two ————————————————————————————————————————————————————————————————————————————————————————————————————————
    print('\n stage 2 training... \n')
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=0)
    for tr, ts in sss.split(dataset.y[train_idx], dataset.y[train_idx]):
        train_idx2 = train_idx[tr]
        valid_idx2 = train_idx[ts]
    train_dataset = dataset.split(torch.tensor(train_idx2).long())
    valid_dataset = dataset.split(torch.tensor(valid_idx2).long())
    # add more samples for small cell types
    train_dataset.data_aug_minor_types(thr=0.02, random_seed=2240)

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=True, collate_fn = collate_func, drop_last=True)
    valid_loader = DataLoader(valid_dataset, batch_size=args.batch_size, num_workers=args.num_workers, shuffle=True, collate_fn = collate_func)

    lr = optimizer.param_groups[0]['lr']
    old_lr = lr
    print('\n stage2 initilize lr: \n',lr)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=args.weight_decay,)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'max', factor=0.1, patience=2, min_lr=0.00001, verbose=True)

    max_metric = float(0)
    max_metric_count = 0
    for epoch_idx, epoch in enumerate(range(args.max_epoch,(args.max_epoch + args.extend_epoch))):
        step = (epoch-1) * iters + epoch_idx * len(train_loader)
        # viz.line(X=[step], Y=[lr], win=viz_lr, update='append')

        if old_lr != lr: # if metric not improve anymore
            max_metric_count = 0
            print('reset max_metric_count to 0 due to updating lr from %.19f to %.19f \n'%(old_lr,lr))
            old_lr = lr

        train_loss, train_embedding, train_true_label = train(model, optimizer, train_loader, epoch)
        train_acc, train_f1 = test(model, train_loader)
        valid_acc, valid_f1 = test(model, valid_loader)

        # visualization
        # viz.line(X=[epoch], Y=[train_loss], win=viz_loss, update='append')
        # viz.line(X=np.array([[epoch, epoch]]), Y=np.array([[train_f1, valid_f1]]), win=viz_f1, update='append')
        # viz.line(X=np.array([[epoch, epoch]]), Y=np.array([[train_acc, valid_acc]]), win=viz_acc, update='append')
        #
        # # Train 2D embedding (tsne)
        # train_embedding = train_embedding.cpu().detach().numpy()
        # ts = TSNE(n_components=2)
        # train_embedding = ts.fit_transform(train_embedding)
        # viz.scatter(X=train_embedding, Y=train_true_label+1, win=viz_embed, update='replace')

        lr = optimizer.param_groups[0]['lr']
        print('[Epoch {}] lr:{}, loss:{}, T_acc:{}, T_f1:{}, V_acc:{}, V_f1:{}'.format(
            epoch, round(lr,6), round(train_loss,4), round(train_acc,4), round(train_f1,4), round(valid_acc,4), round(valid_f1,4)))

        scheduler.step(valid_f1)
        lr = optimizer.param_groups[0]['lr']


        if valid_f1 > max_metric:
            max_metric=valid_f1
            # save best model
            torch.save(model, os.path.join(args.outdir, 'model.pth'))
            max_metric_count = 0
            max_metric = valid_f1
        else:
            if epoch_idx >=2: #ignore first two epochs
                max_metric_count += 1
            if max_metric_count >3:
                print('if max_metric_count >3, break at epoch',epoch)
                break

        if lr <= 0.00001:
            print('if lr <= 0.00001, break at epoch', epoch)
            break

    # Test
    test_acc,test_f1, test_pred_label, test_true_label = test(model, test_loader, predicts=True)
    print('Test F1 = {}, Test Acc = {}'.format(round(test_acc,4), round(test_f1,4)))

    # save final model
    torch.save(model, os.path.join(args.outdir, 'final_model.pth'))

    end = time.time()
    print("Total time = ", end - start) #similary time-consuming
