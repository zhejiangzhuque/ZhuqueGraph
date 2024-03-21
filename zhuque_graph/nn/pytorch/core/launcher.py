import datetime
import os
import time
import deepspeed
import torch
from torch.utils.data import DataLoader
from torch.utils.data.distributed import DistributedSampler
from wg_torch import comm as comm
from wg_torch import embedding_ops as embedding_ops
from wg_torch import graph_ops as graph_ops
from wg_torch.wm_tensor import *
from wholegraph.torch import wholegraph_pytorch as wg
from zhuque_graph.nn.pytorch.model.LinkPredictionWgGNNModel import LinkPredictionWgGNNModel



class WgBaseLauncher(object):
    def __init__(self, args):
        super().__init__()
        self.args = args
    
    def getModel(self, dist_homo_graph):
        pass

    def train_nodeproppred(self, train_data, model):
        train_dataset = graph_ops.NodeClassificationDataset(
            train_data, comm.get_rank(), comm.get_world_size()
        )
        train_sampler = DistributedSampler(
            train_dataset,
            num_replicas=comm.get_world_size(),
            rank=comm.get_rank(),
            shuffle=True,
            drop_last=True,
        )
        train_dataloader = DataLoader(
            train_dataset,
            batch_size=self.args.batchsize,
            num_workers=self.args.dataloaderworkers,
            pin_memory=True,
            sampler=train_sampler,
        )
        if self.args.loss == "ce":
            loss_fcn = torch.nn.CrossEntropyLoss()
        else:
            loss_fcn = torch.nn.BCEWighLogitsLoss()
        epoch = 0
        step = 0
        if comm.is_main_process():
            print("Start training...")
        train_start_time = time.time()
        while epoch < self.args.epochs:
            for i, (idx, label) in enumerate(train_dataloader):
                label = torch.reshape(label, (-1,)).cuda()
                logits = model(idx)
                loss = loss_fcn(logits, label)
                model.backward(loss)
                model.step()
                if comm.is_main_process() and step % 100 == 0:
                    print(
                        "[%s] [TRAIN] Epoch=%d, Step=%d, Loss=%f"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            epoch,
                            step,
                            loss.cpu().item(),
                        )
                    )
                step = step + 1
            epoch = epoch + 1
        comm.synchronize()
        train_end_time = time.time()
        train_time = train_end_time - train_start_time
        epoch_time = train_time / self.args.epochs
        if comm.is_main_process():
            print(
                "[%s] [TRAIN_TIME] Total_train_time=%.2fs, Epoch_time=%.2fs"
                % (
                    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    train_time,
                    epoch_time
                )
            )
    
    @torch.no_grad()
    def valid_test_nodeproppred(self, dataloader, model, name):
        total_correct = 0
        total_valid_sample = 0
        if comm.is_main_process():
            print("%s..." % (name,))
        model.eval()
        from ogb.nodeproppred import Evaluator
        if self.args.source == "ogb":
            evaluator = Evaluator(name=self.args.graph_name)
        else:
            if self.args.metric == "acc":
                evaluator = Evaluator(name="ogbn-products")
            else:
                evaluator = Evaluator(name="ogbn-proteins")
        labels = []
        preds = []
        for i, (idx, label) in enumerate(dataloader):
            labels.append(label)
            logits = model(idx)
            preds.append(logits)
        labels = torch.cat(labels).cpu()
        preds = torch.cat(preds).cpu()
        if self.args.metric == "acc":
            preds = torch.argmax(preds, 1).reshape(-1, 1)
            acc = evaluator.eval(
                {
                    "y_true": labels,
                    "y_pred": preds,
                }
            )["acc"]
            if comm.is_main_process():
                print(
                    "[%s] [%s] Accuracy=%.2f%%"
                    % (
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        name,
                        acc,
                    )
                )
        else:
            rocauc = evaluator.eval(
                {
                    "y_true": labels,
                    "y_pred": preds,
                }
            )["rocauc"]
            if comm.is_main_process():
                print(
                    "[%s] [%s] ROCAUC=%.2f%%"
                    % (
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        name,
                        rocauc,
                    )
                )
    
    def valid_nodeproppred(self, valid_data, model):
        valid_dataset = graph_ops.NodeClassificationDataset(
            valid_data, 0, 1
        )
        valid_dataloader = DataLoader(
            valid_dataset,
            batch_size=self.args.batchsize,
            shuffle=False,
            pin_memory=True
        )
        self.valid_test_nodeproppred(valid_dataloader, model, "VALID")
    
    def test_nodeproppred(self, test_data, model):
        test_dataset = graph_ops.NodeClassificationDataset(
            test_data, 0, 1
        )
        test_dataloader = DataLoader(
            test_dataset,
            batch_size=self.args.batchsize,
            shuffle=False,
            pin_memory=True
        )
        self.valid_test_nodeproppred(test_dataloader, model, "TEST")
    
    def train_linkproppred(self, edge_split, model:LinkPredictionWgGNNModel, dist_homo_graph):
        if self.args.loss == "ce":
            loss_fcn = torch.nn.CrossEntropyLoss()
        else:
            loss_fcn = torch.nn.BCEWithLogitsLoss()
        epoch = 0
        step = 0
        if comm.is_main_process():
            print("Start training...")
        train_start_time = time.time()
        iters_in_epoch = dist_homo_graph.start_iter(self.args.batchsize)
        while epoch < self.args.epochs:
            for i in range(iters_in_epoch):
                src_nid, pos_dst_nid = dist_homo_graph.get_train_edge_batch(i)
                neg_dst_nid = dist_homo_graph.per_source_negative_sample(src_nid)
                pos_score, neg_score = model(src_nid, pos_dst_nid, neg_dst_nid)
                pos_label = torch.ones_like(pos_score)
                neg_label = torch.zeros_like(neg_score)
                score = torch.cat([pos_score, neg_score])
                label = torch.cat([pos_label, neg_label])
                loss = loss_fcn(score, label)
                model.backward(loss)
                model.step()
                if comm.is_main_process() and step % 100 == 0:
                    print(
                        "[%s] [TRAIN] Epoch=%d, Step=%d, Loss=%f"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            epoch,
                            step,
                            loss.cpu().item(),
                        )
                    )
                step = step + 1
            epoch = epoch + 1
        comm.synchronize()
        train_end_time = time.time()
        train_time = train_end_time - train_start_time
        epoch_time = train_time / self.args.epochs
        if comm.is_main_process():
            print(
                "[%s] [TRAIN_TIME] Total_train_time=%.2fs, Epoch_time=%.2fs"
                % (
                    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    train_time,
                    epoch_time
                )
            )
    
    @torch.no_grad()
    def valid_test_linkproppred(self, edge_split, model:LinkPredictionWgGNNModel, dist_homo_graph, name):
        if comm.is_main_process():
            print("%s..." % (name,))
        model.eval()
        embedding = dist_homo_graph.node_feat
        node_feats = [None, None]
        use_chunked = True
        use_host_memory = False
        if self.args.use_host_mem:
            use_chunked = False
            use_host_memory = True
        wm_tensor_type = get_intra_node_wm_tensor_type(
            use_chunked,
            use_host_memory
        )
        node_feats[0] = create_wm_tensor(
            dist_homo_graph.wm_comm,
            [embedding.shape[0], self.args.hiddensize],
            [],
            embedding.dtype,
            wm_tensor_type,
        )
        if self.args.layernum > 1:
            node_feats[1] = create_wm_tensor(
                dist_homo_graph.wm_comm,
                [embedding.shape[0], self.args.hiddensize],
                [],
                embedding.dtype,
                wm_tensor_type,
            )
        output_feat = node_feats[0]
        input_feat = embedding
        del embedding
        for i in range(self.args.layernum):
            model.fullbatch_single_layer_forward(
                dist_homo_graph, i, input_feat, output_feat, self.args.batchsize
            )
            wg.barrier(dist_homo_graph.wm_comm)
            input_feat = output_feat
            output_feat = node_feats[(i + 1) % 2]
        del output_feat
        del node_feats[1]
        del node_feats[0]
        del node_feats
        embedding_lookup_fn = embedding_ops.EmbeddingLookupFn.apply
        from ogb.linkproppred import Evaluator
        if self.args.source == "ogb":
            evaluator = Evaluator(name=self.args.graph_name)
        else:
            if self.args.metric == "hits100":
                evaluator = Evaluator(name="ogbl-ppa")
            elif self.args.metric == "hits50":
                evaluator = Evaluator(name="ogbl-collab")
            elif self.args.metric == "hits20":
                evaluator = Evaluator(name="ogbl-ddi")
            elif self.args.metric == "mrr":
                evaluator = Evaluator(name="ogbl-citation2")
            else:
                evaluator = Evaluator(name="ogbl-vessel")
        if self.args.graph_name == "ogbl-citation2":
            src = torch.from_numpy(edge_split["source_node"]).cuda()
            dst = torch.from_numpy(edge_split["target_node"]).cuda()
            neg_dst = torch.from_numpy(edge_split["target_node_neg"]).cuda()
            preds = []
            for start in range(0, src.shape[0], self.args.batchsize):
                end = min(start + self.args.batchsize, src.shape[0])
                all_dst = torch.cat([dst[start:end, None], neg_dst[start:end]], 1)
                h_src = embedding_lookup_fn(src[start:end], input_feat)[:, None, :]
                h_dst = embedding_lookup_fn(all_dst.view(-1), input_feat).view(*all_dst.shape, -1)
                pred = model.predict(h_src, h_dst).squeeze(-1)
                preds += [pred]
            all_preds = torch.cat(preds)
            pos_preds = all_preds[:, :1].squeeze(1)
            neg_preds = all_preds[:, 1:]
            rr = evaluator.eval(
                    {
                        "y_pred_pos": pos_preds,
                        "y_pred_neg": neg_preds,
                    }
                )["mrr_list"]
            mrr = rr.mean().item()
            if comm.is_main_process():
                print(
                    "[%s] [%s] MRR=%.2f%%"
                    % (
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        name,
                        mrr
                    )
                )
        else:
            if self.args.graph_name == "ogbl-vessel":
                pos_edge = edge_split["edge"].cuda()
                neg_edge = edge_split["edge_neg"].cuda()
            else:
                pos_edge = torch.from_numpy(edge_split["edge"]).cuda()
                neg_edge = torch.from_numpy(edge_split["edge_neg"]).cuda()
            pos_preds = []
            neg_preds = []
            for start in range(0, pos_edge.shape[0], self.args.batchsize):
                end = min(start + self.args.batchsize, pos_edge.shape[0])
                h_pos_src = embedding_lookup_fn(pos_edge[start:end, 0], input_feat)
                h_pos_dst = embedding_lookup_fn(pos_edge[start:end, 1], input_feat)
                h_neg_src = embedding_lookup_fn(neg_edge[start:end, 0], input_feat)
                h_neg_dst = embedding_lookup_fn(neg_edge[start:end, 1], input_feat)
                pos_pred = model.predict(h_pos_src, h_pos_dst).squeeze(-1)
                neg_pred = model.predict(h_neg_src, h_neg_dst).squeeze(-1)
                pos_preds += [pos_pred]
                neg_preds += [neg_pred]
            all_pos_preds = torch.cat(pos_preds)
            all_neg_preds = torch.cat(neg_preds)
            if self.args.metric == "hits100":
                hits = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["hits@100"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] Hits@100=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            hits
                        )
                    )
            elif self.args.metric == "hits50":
                hits = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["hits@50"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] Hits@50=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            hits
                        )
                    )
            elif self.args.metric == "hits20":
                hits = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["hits@20"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] Hits@20=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            hits
                        )
                    )
            else:
                rocauc = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["rocauc"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] ROCAUC=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            rocauc
                        )
                    )
    
    def valid_linkproppred(self, valid_edge, model, dist_homo_graph):
        self.valid_test_linkproppred(valid_edge, model, dist_homo_graph, "VALID")
    
    def test_linkproppred(self, test_edge, model, dist_homo_graph):
        self.valid_test_linkproppred(test_edge, model, dist_homo_graph, "TEST")

    def train_and_test(self):
        print(self.args)
        print(self.args.output_model_name)
        wg.init_lib()
        torch.set_num_threads(1)
        
        if "MASTER_ADDR" not in os.environ:
            os.environ["MASTER_ADDR"] = "localhost"
        if "MASTER_PORT" not in os.environ:
            os.environ["MASTER_PORT"] = "12345"
        master_addr = os.environ["MASTER_ADDR"]
        master_port = os.environ["MASTER_PORT"]
        world_rank = int(os.environ["RANK"])
        world_size = int(os.environ['WORLD_SIZE'])
        local_rank = int(os.environ["LOCAL_RANK"])
        local_size = int(os.environ["LOCAL_SIZE"])
        torch.cuda.set_device(local_rank)
        deepspeed.init_distributed(dist_backend="nccl", init_method=f'tcp://{master_addr}:{master_port}', world_size=world_size, rank=world_rank)

        wm_comm = create_intra_node_communicator(world_rank, world_size, local_size)
        wm_embedding_comm = None
        if self.args.use_nccl:
            if comm.is_main_process():
                print("Using nccl embeddings.")
            wm_embedding_comm = create_global_communicator(world_rank, world_size)

        train_data = None
        valid_data = None
        test_data = None
        edge_split = None
        if self.args.task == "node":
            train_data, valid_data, test_data = graph_ops.load_pickle_data(self.args.root_dir, self.args.graph_name, True)
            if comm.is_main_process():
                print("Task=node-classification")
                print("Framework=%s, Model=%s" % (self.args.framework, self.args.model))
        elif self.args.task == "link":
            edge_split = graph_ops.load_pickle_link_pred_data(self.args.root_dir, self.args.graph_name, True)
            train_edge = edge_split["train"]
            valid_edge = edge_split["valid"]
            test_edge = edge_split["test"]
            if comm.is_main_process():
                print("Task=link-prediction")
                print("Framework=%s, Model=%s" % (self.args.framework, self.args.model))
        
        dist_homo_graph = graph_ops.HomoGraph()
        use_chunked = True
        use_host_memory = False
        if self.args.use_host_mem:
            use_chunked = False
            use_host_memory = True
            if world_rank == 0:
                print("Using host memory.")
        dist_homo_graph.load(
            self.args.root_dir,
            self.args.graph_name,
            wm_comm,
            use_chunked,
            use_host_memory,
            wm_embedding_comm,
            feat_dtype=None,
            id_dtype=None,
            ignore_embeddings=None,
        )
        print("Rank=%d, Graph loaded." % (world_rank,))
        raw_model = self.getModel(dist_homo_graph)
        print("Rank=%d, model created." % (world_rank,))
        raw_model.cuda()
        print("Rank=%d, model moved to cuda." % (world_rank,))
        parameters = filter(lambda p: p.requires_grad, raw_model.parameters())
        model_engine, optimizer, _, _ = deepspeed.initialize(
                args=self.args, model=raw_model, model_parameters=parameters)
        model = model_engine
        print("Rank=%d, ddp model created." % (world_rank,))

        if self.args.task == "node":
            self.train_nodeproppred(train_data, model)
            self.valid_nodeproppred(valid_data, model)
            self.test_nodeproppred(test_data, model)
        elif self.args.task == "link":
            self.train_linkproppred(train_edge, model, dist_homo_graph)
            self.valid_linkproppred(valid_edge, model, dist_homo_graph)
            self.test_linkproppred(test_edge, model, dist_homo_graph)
        print(self.args)
        torch.save(model.state_dict(), self.args.output_model_name)
        wg.finalize_lib()
        print("Rank=%d, wholegraph shutdown." % (world_rank,))


class NodeClassificationWgLauncher(object):
    def __init__(self, args):
        super().__init__()
        self.args = args
    
    def getModel(self, dist_homo_graph):
        pass
        
    def train(self, train_data, model):
        train_dataset = graph_ops.NodeClassificationDataset(
            train_data, comm.get_rank(), comm.get_world_size()
        )
        train_sampler = DistributedSampler(
            train_dataset,
            num_replicas=comm.get_world_size(),
            rank=comm.get_rank(),
            shuffle=True,
            drop_last=True,
        )
        train_dataloader = DataLoader(
            train_dataset,
            batch_size=self.args.batchsize,
            num_workers=self.args.dataloaderworkers,
            pin_memory=True,
            sampler=train_sampler,
        )
        if self.args.loss == "ce":
            loss_fcn = torch.nn.CrossEntropyLoss()
        elif self.args.loss == "bcelogit":
            loss_fcn = torch.nn.BCEWighLogitsLoss()
        epoch = 0
        step = 0
        if comm.is_main_process():
            print("Start training...")
        train_start_time = time.time()
        while epoch < self.args.epochs:
            for i, (idx, label) in enumerate(train_dataloader):
                label = torch.reshape(label, (-1,)).cuda()
                logits = model(idx)
                loss = loss_fcn(logits, label)
                model.backward(loss)
                model.step()
                if comm.is_main_process() and step % 100 == 0:
                    print(
                        "[%s] [TRAIN] Epoch=%d, Step=%d, Loss=%f"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            epoch,
                            step,
                            loss.cpu().item(),
                        )
                    )
                step = step + 1
            epoch = epoch + 1
        comm.synchronize()
        train_end_time = time.time()
        train_time = train_end_time - train_start_time
        epoch_time = train_time / self.args.epochs
        if comm.is_main_process():
            print(
                "[%s] [TRAIN_TIME] Total_train_time=%.2fs, Epoch_time=%.2fs"
                % (
                    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    train_time,
                    epoch_time
                )
            )
    
    @torch.no_grad()
    def valid_test(self, dataloader, model, name):
        total_correct = 0
        total_valid_sample = 0
        if comm.is_main_process():
            print("%s..." % (name,))
        model.eval()
        from ogb.nodeproppred import Evaluator
        if self.args.source == "ogb":
            evaluator = Evaluator(name=self.args.graph_name)
        else:
            if self.args.metric == "acc":
                evaluator = Evaluator(name="ogbn-products")
            elif self.args.metric == "rocauc":
                evaluator = Evaluator(name="ogbn-proteins")
        labels = []
        preds = []
        for i, (idx, label) in enumerate(dataloader):
            labels.append(label)
            logits = model(idx)
            preds.append(logits)
        labels = torch.cat(labels).cpu()
        preds = torch.cat(preds).cpu()
        if self.args.metric == "acc":
            preds = torch.argmax(preds, 1).reshape(-1, 1)
            acc = evaluator.eval(
                {
                    "y_true": labels,
                    "y_pred": preds,
                }
            )["acc"]
            if comm.is_main_process():
                print(
                    "[%s] [%s] Accuracy=%.2f%%"
                    % (
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        name,
                        acc,
                    )
                )
        elif self.args.metric == "rocauc":
            rocauc = evaluator.eval(
                {
                    "y_true": labels,
                    "y_pred": preds,
                }
            )["rocauc"]
            if comm.is_main_process():
                print(
                    "[%s] [%s] ROCAUC=%.2f%%"
                    % (
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        name,
                        rocauc,
                    )
                )
    
    def valid(self, valid_data, model):
        valid_dataset = graph_ops.NodeClassificationDataset(
            valid_data, 0, 1
        )
        valid_dataloader = DataLoader(
            valid_dataset,
            batch_size=self.args.batchsize,
            shuffle=False,
            pin_memory=True
        )
        self.valid_test(valid_dataloader, model, "VALID")
    
    def test(self, test_data, model):
        test_dataset = graph_ops.NodeClassificationDataset(
            test_data, 0, 1
        )
        test_dataloader = DataLoader(
            test_dataset,
            batch_size=self.args.batchsize,
            shuffle=False,
            pin_memory=True
        )
        self.valid_test(test_dataloader, model, "TEST")
    
    def train_and_test(self):
        wg.init_lib()
        torch.set_num_threads(1)
        
        if "MASTER_ADDR" not in os.environ:
            os.environ["MASTER_ADDR"] = "localhost"
        if "MASTER_PORT" not in os.environ:
            os.environ["MASTER_PORT"] = "12345"
        master_addr = os.environ["MASTER_ADDR"]
        master_port = os.environ["MASTER_PORT"]
        world_rank = int(os.environ["RANK"])
        world_size = int(os.environ['WORLD_SIZE'])
        local_rank = int(os.environ["LOCAL_RANK"])
        local_size = int(os.environ["LOCAL_SIZE"])
        torch.cuda.set_device(local_rank)
        deepspeed.init_distributed(dist_backend="nccl", init_method=f'tcp://{master_addr}:{master_port}', world_size=world_size, rank=world_rank)

        wm_comm = create_intra_node_communicator(world_rank, world_size, local_size)
        wm_embedding_comm = None
        if self.args.use_nccl:
            if comm.is_main_process():
                print("Using nccl embeddings.")
            wm_embedding_comm = create_global_communicator(world_rank, world_size)

        train_data, valid_data, test_data = graph_ops.load_pickle_data(self.args.root_dir, self.args.graph_name, True)
        if comm.is_main_process():
            print("Task=node-classification")
            print("Framework=%s, Model=%s" % (self.args.framework, self.args.model))
        
        dist_homo_graph = graph_ops.HomoGraph()
        use_chunked = True
        use_host_memory = False
        if self.args.use_host_mem:
            use_chunked = False
            use_host_memory = True
            if world_rank == 0:
                print("Using host memory.")
        dist_homo_graph.load(
            self.args.root_dir,
            self.args.graph_name,
            wm_comm,
            use_chunked,
            use_host_memory,
            wm_embedding_comm,
            feat_dtype=None,
            id_dtype=None,
            ignore_embeddings=None,
        )
        print("Rank=%d, Graph loaded." % (world_rank,))
        model = self.getModel(dist_homo_graph)
        print("Rank=%d, model created." % (world_rank,))
        model.cuda()
        print("Rank=%d, model moved to cuda." % (world_rank,))
        model_parameters = filter(lambda p: p.requires_grad, model.parameters())
        model_engine, optimizer, _, _ = deepspeed.initialize(
                args=self.args, model=model, model_parameters=model_parameters)
        print("Rank=%d, ddp model created." % (world_rank,))

        self.train(train_data, model_engine)
        self.valid(valid_data, model_engine)
        self.test(test_data, model_engine)

        #torch.save(model.state_dict(), 'gnn_model.pth')
        wg.finalize_lib()
        print("Rank=%d, wholegraph shutdown." % (world_rank,))


class LinkPredictionWgLauncher(object):
    def __init__(self, args):
        super().__init__()
        self.args = args
    
    def getModel(self, dist_homo_graph):
        pass
    
    def train(self, edge_split, model: LinkPredictionWgGNNModel, dist_homo_graph):
        if self.args.loss == "ce":
            loss_fcn = torch.nn.CrossEntropyLoss()
        else:
            loss_fcn = torch.nn.BCEWithLogitsLoss()
        epoch = 0
        step = 0
        if comm.is_main_process():
            print("Start training...")
        train_start_time = time.time()
        iters_in_epoch = dist_homo_graph.start_iter(self.args.batchsize)
        while epoch < self.args.epochs:
            for i in range(iters_in_epoch):
                src_nid, pos_dst_nid = dist_homo_graph.get_train_edge_batch(i)
                neg_dst_nid = dist_homo_graph.per_source_negative_sample(src_nid)
                pos_score, neg_score = model(src_nid, pos_dst_nid, neg_dst_nid)
                pos_label = torch.ones_like(pos_score)
                neg_label = torch.zeros_like(neg_score)
                score = torch.cat([pos_score, neg_score])
                label = torch.cat([pos_label, neg_label])
                loss = loss_fcn(score, label)
                model.backward(loss)
                model.step()
                if comm.is_main_process() and step % 100 == 0:
                    print(
                        "[%s] [TRAIN] Epoch=%d, Step=%d, Loss=%f"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            epoch,
                            step,
                            loss.cpu().item(),
                        )
                    )
                step = step + 1
            epoch = epoch + 1
        comm.synchronize()
        train_end_time = time.time()
        train_time = train_end_time - train_start_time
        epoch_time = train_time / self.args.epochs
        if comm.is_main_process():
            print(
                "[%s] [TRAIN_TIME] Total_train_time=%.2fs, Epoch_time=%.2fs"
                % (
                    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    train_time,
                    epoch_time
                )
            )
    
    @torch.no_grad()
    def valid_test(self, edge_split, model:LinkPredictionWgGNNModel, dist_homo_graph, name):
        if comm.is_main_process():
            print("%s..." % (name,))
        model.eval()
        embedding = dist_homo_graph.node_feat
        node_feats = [None, None]
        use_chunked = True
        use_host_memory = False
        if self.args.use_host_mem:
            use_chunked = False
            use_host_memory = True
        wm_tensor_type = get_intra_node_wm_tensor_type(
            use_chunked,
            use_host_memory
        )
        node_feats[0] = create_wm_tensor(
            dist_homo_graph.wm_comm,
            [embedding.shape[0], self.args.hiddensize],
            [],
            embedding.dtype,
            wm_tensor_type,
        )
        if self.args.layernum > 1:
            node_feats[1] = create_wm_tensor(
                dist_homo_graph.wm_comm,
                [embedding.shape[0], self.args.hiddensize],
                [],
                embedding.dtype,
                wm_tensor_type,
            )
        output_feat = node_feats[0]
        input_feat = embedding
        del embedding
        for i in range(self.args.layernum):
            model.fullbatch_single_layer_forward(
                dist_homo_graph, i, input_feat, output_feat, self.args.batchsize
            )
            wg.barrier(dist_homo_graph.wm_comm)
            input_feat = output_feat
            output_feat = node_feats[(i + 1) % 2]
        del output_feat
        del node_feats[1]
        del node_feats[0]
        del node_feats
        embedding_lookup_fn = embedding_ops.EmbeddingLookupFn.apply
        from ogb.linkproppred import Evaluator
        if self.args.source == "ogb":
            evaluator = Evaluator(name=self.args.graph_name)
        else:
            if self.args.metric == "hits100":
                evaluator = Evaluator(name="ogbl-ppa")
            elif self.args.metric == "hits50":
                evaluator = Evaluator(name="ogbl-collab")
            elif self.args.metric == "hits20":
                evaluator = Evaluator(name="ogbl-ddi")
            elif self.args.metric == "mrr":
                evaluator = Evaluator(name="ogbl-citation2")
            elif self.args.metric == "rocauc":
                evaluator = Evaluator(name="ogbl-vessel")
        if self.args.graph_name == "ogbl-citation2":
            src = torch.from_numpy(edge_split["source_node"]).cuda()
            dst = torch.from_numpy(edge_split["target_node"]).cuda()
            neg_dst = torch.from_numpy(edge_split["target_node_neg"]).cuda()
            preds = []
            for start in range(0, src.shape[0], self.args.batchsize):
                end = min(start + self.args.batchsize, src.shape[0])
                all_dst = torch.cat([dst[start:end, None], neg_dst[start:end]], 1)
                h_src = embedding_lookup_fn(src[start:end], input_feat)[:, None, :]
                h_dst = embedding_lookup_fn(all_dst.view(-1), input_feat).view(*all_dst.shape, -1)
                pred = model.predict(h_src, h_dst).squeeze(-1)
                preds += [pred]
            all_preds = torch.cat(preds)
            pos_preds = all_preds[:, :1].squeeze(1)
            neg_preds = all_preds[:, 1:]
            rr = evaluator.eval(
                    {
                        "y_pred_pos": pos_preds,
                        "y_pred_neg": neg_preds,
                    }
                )["mrr_list"]
            mrr = rr.mean().item()
            if comm.is_main_process():
                print(
                    "[%s] [%s] MRR=%.2f%%"
                    % (
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        name,
                        mrr
                    )
                )
        else:
            if self.args.graph_name == "ogbl-vessel":
                pos_edge = edge_split["edge"].cuda()
                neg_edge = edge_split["edge_neg"].cuda()
            else:
                pos_edge = torch.from_numpy(edge_split["edge"]).cuda()
                neg_edge = torch.from_numpy(edge_split["edge_neg"]).cuda()
            pos_preds = []
            neg_preds = []
            for start in range(0, pos_edge.shape[0], self.args.batchsize):
                end = min(start + self.args.batchsize, pos_edge.shape[0])
                h_pos_src = embedding_lookup_fn(pos_edge[start:end, 0], input_feat)
                h_pos_dst = embedding_lookup_fn(pos_edge[start:end, 1], input_feat)
                h_neg_src = embedding_lookup_fn(neg_edge[start:end, 0], input_feat)
                h_neg_dst = embedding_lookup_fn(neg_edge[start:end, 1], input_feat)
                pos_pred = model.predict(h_pos_src, h_pos_dst).squeeze(-1)
                neg_pred = model.predict(h_neg_src, h_neg_dst).squeeze(-1)
                pos_preds += [pos_pred]
                neg_preds += [neg_pred]
            all_pos_preds = torch.cat(pos_preds)
            all_neg_preds = torch.cat(neg_preds)
            if self.args.metric == "hits100":
                hits = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["hits@100"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] Hits@100=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            hits
                        )
                    )
            elif self.args.metric == "hits50":
                hits = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["hits@50"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] Hits@50=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            hits
                        )
                    )
            elif self.args.metric == "hits20":
                hits = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["hits@20"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] Hits@20=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            hits
                        )
                    )
            elif self.args.metric == "rocauc":
                rocauc = evaluator.eval(
                        {
                            "y_pred_pos": all_pos_preds,
                            "y_pred_neg": all_neg_preds,
                        }
                    )["rocauc"]
                if comm.is_main_process():
                    print(
                        "[%s] [%s] ROCAUC=%.2f%%"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            name,
                            rocauc
                        )
                    )
    
    def valid(self, valid_edge, model, dist_homo_graph):
        self.valid_test(valid_edge, model, dist_homo_graph, "VALID")
    
    def test(self, test_edge, model, dist_homo_graph):
        self.valid_test(test_edge, model, dist_homo_graph, "TEST")
    
    def train_and_test(self):
        wg.init_lib()
        torch.set_num_threads(1)
        
        if "MASTER_ADDR" not in os.environ:
            os.environ["MASTER_ADDR"] = "localhost"
        if "MASTER_PORT" not in os.environ:
            os.environ["MASTER_PORT"] = "12345"
        master_addr = os.environ["MASTER_ADDR"]
        master_port = os.environ["MASTER_PORT"]
        world_rank = int(os.environ["RANK"])
        world_size = int(os.environ['WORLD_SIZE'])
        local_rank = int(os.environ["LOCAL_RANK"])
        local_size = int(os.environ["LOCAL_SIZE"])
        torch.cuda.set_device(local_rank)
        deepspeed.init_distributed(dist_backend="nccl", init_method=f'tcp://{master_addr}:{master_port}', world_size=world_size, rank=world_rank)

        wm_comm = create_intra_node_communicator(world_rank, world_size, local_size)
        wm_embedding_comm = None
        if self.args.use_nccl:
            if comm.is_main_process():
                print("Using nccl embeddings.")
            wm_embedding_comm = create_global_communicator(world_rank, world_size)

        edge_split = graph_ops.load_pickle_link_pred_data(self.args.root_dir, self.args.graph_name, True)
        train_edge = edge_split["train"]
        valid_edge = edge_split["valid"]
        test_edge = edge_split["test"]
        if comm.is_main_process():
            print("Task=link-prediction")
            print("Framework=%s, Model=%s" % (self.args.framework, self.args.model))
        
        dist_homo_graph = graph_ops.HomoGraph()
        use_chunked = True
        use_host_memory = False
        if self.args.use_host_mem:
            use_chunked = False
            use_host_memory = True
            if world_rank == 0:
                print("Using host memory.")
        dist_homo_graph.load(
            self.args.root_dir,
            self.args.graph_name,
            wm_comm,
            use_chunked,
            use_host_memory,
            wm_embedding_comm,
            feat_dtype=None,
            id_dtype=None,
            ignore_embeddings=None,
        )
        print("Rank=%d, Graph loaded." % (world_rank,))
        model = self.getModel(dist_homo_graph)
        print("Rank=%d, model created." % (world_rank,))
        model.cuda()
        print("Rank=%d, model moved to cuda." % (world_rank,))
        model_parameters = filter(lambda p: p.requires_grad, model.parameters())
        model_engine, optimizer, _, _ = deepspeed.initialize(
                args=self.args, model=model, model_parameters=model_parameters)
        print("Rank=%d, ddp model created." % (world_rank,))

        self.train(train_edge, model_engine, dist_homo_graph)
        self.valid(valid_edge, model_engine, dist_homo_graph)
        self.test(test_edge, model_engine, dist_homo_graph)

        wg.finalize_lib()
        print("Rank=%d, wholegraph shutdown." % (world_rank,))