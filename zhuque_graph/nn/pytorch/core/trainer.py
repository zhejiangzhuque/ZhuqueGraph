import datetime
import os
import time

import apex
import torch
from apex.parallel import DistributedDataParallel as DDP
from mpi4py import MPI
from torch.utils.data import DataLoader
from wg_torch import comm as comm
from wg_torch import graph_ops as graph_ops
from wg_torch.wm_tensor import *

from wholegraph.torch import wholegraph_pytorch as wg

class WgBaseTrainer(object):
    def __init__(self, options):
        super().__init__()
        self.options = options

    def get_train_step(self, sample_count, epochs, batch_size, global_size):
        return sample_count * epochs // (batch_size * global_size)

    def create_test_dataset(self, data_tensor_dict):
        return DataLoader(
            dataset=graph_ops.NodeClassificationDataset(data_tensor_dict, 0, 1),
            batch_size=(self.options.batchsize + 3) // 4,
            shuffle=False,
            pin_memory=True,
        )

    def valid_test(self, dataloader, model, name):
        total_correct = 0
        total_valid_sample = 0
        if comm.get_rank() == 0:
            print("%s..." % (name,))
        for i, (idx, label) in enumerate(dataloader):
            label = torch.reshape(label, (-1,)).cuda()
            model.eval()
            logits = model(idx)
            pred = torch.argmax(logits, 1)
            correct = (pred == label).sum()
            total_correct += correct.cpu()
            total_valid_sample += label.shape[0]
        if comm.get_rank() == 0:
            print(
                "[%s] [%s] accuracy=%5.2f%%"
                % (
                    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    name,
                    100.0 * total_correct / total_valid_sample,
                )
            )

    def valid(self, valid_dataloader, model):
        self.valid_test(valid_dataloader, model, "VALID")

    def test(self, test_data, model):
        test_dataloader = self.create_test_dataset(data_tensor_dict=test_data)
        self.valid_test(test_dataloader, model, "TEST")

    def train_torch_sampler(self, train_data, valid_data, model, optimizer):
        print("start training...")
        train_dataset = graph_ops.NodeClassificationDataset(
            train_data, comm.get_rank(), comm.get_world_size()
        )
        valid_dataset = graph_ops.NodeClassificationDataset(
            valid_data, comm.get_rank(), comm.get_world_size()
        )
        train_sampler = torch.utils.data.distributed.DistributedSampler(
            train_dataset,
            num_replicas=comm.get_world_size(),
            rank=comm.get_rank(),
            shuffle=True,
            drop_last=True,
        )
        valid_sampler = torch.utils.data.distributed.DistributedSampler(
            valid_dataset, num_replicas=1, rank=0, shuffle=False, drop_last=False
        )

        train_dataloader = torch.utils.data.DataLoader(
            train_dataset,
            batch_size=self.options.batchsize,
            num_workers=self.options.dataloaderworkers,
            pin_memory=True,
            sampler=train_sampler,
        )
        valid_dataloader = torch.utils.data.DataLoader(
            valid_dataset,
            batch_size=self.options.batchsize,
            num_workers=self.options.dataloaderworkers,
            pin_memory=True,
            sampler=valid_sampler,
        )
        est_epoch_steps = self.get_train_step(
            len(train_data["idx"]), 1, self.options.batchsize, comm.get_world_size()
        )
        if comm.get_rank() == 0:
            print(
                "Estimated epoch steps=%d, total steps=%d"
                % (est_epoch_steps, est_epoch_steps * self.options.epochs)
            )
        train_step = 0
        epoch = 0
        loss_fcn = torch.nn.CrossEntropyLoss()
        skip_count = 8
        skip_epoch_time = time.time()
        train_start_time = time.time()
        while epoch < self.options.epochs:
            for i, (idx, label) in enumerate(train_dataloader):
                label = torch.reshape(label, (-1,)).cuda()
                optimizer.zero_grad()
                model.train()
                logits = model(idx)
                loss = loss_fcn(logits, label)
                loss.backward()
                optimizer.step()
                if comm.get_rank() == 0 and train_step % 100 == 0:
                    print(
                        "[%s] [LOSS] step=%d, loss=%f"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            train_step,
                            loss.cpu().item(),
                        )
                    )
                train_step = train_step + 1
            epoch = epoch + 1
            if epoch == skip_count:
                skip_epoch_time = time.time()
        comm.synchronize()
        train_end_time = time.time()
        train_time = train_end_time - train_start_time
        if comm.get_rank() == 0:
            print(
                "[%s] [TRAIN_TIME] core time is %.2f seconds"
                % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), train_time)
            )
            print(
                "[EPOCH_TIME] %.2f seconds"
                % ((train_end_time - skip_epoch_time) / (self.options.epochs - skip_count),)
            )
        self.valid(valid_dataloader, model)

    def create_train_dataset(self,data_tensor_dict, rank, size):
        return DataLoader(
            dataset=graph_ops.NodeClassificationDataset(data_tensor_dict, rank, size),
            batch_size=self.options.batchsize,
            shuffle=True,
            num_workers=self.options.dataloaderworkers,
            pin_memory=True,
        )

    def create_valid_dataset(self,data_tensor_dict):
        return DataLoader(
            dataset=graph_ops.NodeClassificationDataset(data_tensor_dict, 0, 1),
            batch_size=(self.options.batchsize + 3) // 4,
            shuffle=False,
            pin_memory=True,
        )

    def getModel(self, dist_homo_graph):
        pass

    def train(self, train_data, valid_data, model, optimizer):
        if comm.get_rank() == 0:
            print("start training...")
        train_dataloader = self.create_train_dataset(
            data_tensor_dict=train_data, rank=comm.get_rank(), size=comm.get_world_size()
        )
        valid_dataloader = self.create_valid_dataset(data_tensor_dict=valid_data)
        total_steps = self.get_train_step(
            len(train_data["idx"]), self.options.epochs, self.options.batchsize, comm.get_world_size()
        )
        if comm.get_rank() == 0:
            print(
                "epoch=%d total_steps=%d"
                % (
                    self.options.epochs,
                    total_steps,
                )
            )
        train_step = 0
        epoch = 0
        loss_fcn = torch.nn.CrossEntropyLoss()
        skip_world_size_epoch_time = 0
        train_start_time = time.time()
        while train_step < total_steps:
            if epoch == 1:
                skip_world_size_epoch_time = time.time()
            for i, (idx, label) in enumerate(train_dataloader):
                if train_step >= total_steps:
                    break
                label = torch.reshape(label, (-1,)).cuda()
                optimizer.zero_grad()
                model.train()
                logits = model(idx)
                loss = loss_fcn(logits, label)
                loss.backward()
                optimizer.step()
                if comm.get_rank() == 0 and train_step % 100 == 0:
                    print(
                        "[%s] [LOSS] step=%d, loss=%f"
                        % (
                            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            train_step,
                            loss.cpu().item(),
                        )
                    )
                train_step = train_step + 1
            epoch = epoch + 1
        comm.synchronize()
        train_end_time = time.time()
        train_time = train_end_time - train_start_time
        if comm.get_rank() == 0:
            print(
                "[%s] [TRAIN_TIME] core time is %.2f seconds"
                % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), train_time)
            )
            if self.options.epochs <= comm.get_world_size():
                print(
                    "[EPOCH_TIME] %.2f seconds, maybe large due to not enough epoch skipped."
                    % ((train_end_time - train_start_time) / self.options.epochs,)
                )
            else:
                print(
                    "[EPOCH_TIME] %.2f seconds"
                    % (
                        (train_end_time - skip_world_size_epoch_time)
                        / (self.options.epochs - comm.get_world_size()),
                    )
                )
        self.valid(valid_dataloader, model)


    def train_and_test(self):
        wg.init_lib()
        torch.set_num_threads(1)
        comma = MPI.COMM_WORLD
        shared_comma = comma.Split_type(MPI.COMM_TYPE_SHARED)
        os.environ["RANK"] = str(comma.Get_rank())
        os.environ["WORLD_SIZE"] = str(comma.Get_size())
        # slurm in Selene has MASTER_ADDR env
        if "MASTER_ADDR" not in os.environ:
            os.environ["MASTER_ADDR"] = "localhost"
        if "MASTER_PORT" not in os.environ:
            os.environ["MASTER_PORT"] = "12335"
        local_rank = shared_comma.Get_rank()
        local_size = shared_comma.Get_size()
        dev_count = torch.cuda.device_count()
        print("Rank=%d,  local_rank=%d, local_size=%d, dev_count=%d" % (comma.Get_rank(), local_rank, local_size, dev_count))
        assert dev_count > 0
        assert local_size <= dev_count
        torch.cuda.set_device(local_rank)
        torch.distributed.init_process_group(backend="nccl", init_method="env://")
        wm_comm = create_intra_node_communicator(
            comma.Get_rank(), comma.Get_size(), local_size
        )
        wm_embedding_comm = None
        if self.options.use_nccl:
            if comma.Get_rank() == 0:
                print("Using nccl embeddings.")
            wm_embedding_comm = create_global_communicator(
                comma.Get_rank(), comma.Get_size()
            )
        if comma.Get_rank() == 0:
            print("Framework=%s, Model=%s" % (self.options.framework, self.options.model))

        train_data, valid_data, test_data = graph_ops.load_pickle_data(
            self.options.root_dir, self.options.graph_name, True
        )

        dist_homo_graph = graph_ops.HomoGraph()
        use_chunked = True
        use_host_memory = False
        if self.options.use_host_mem:
            use_chunked = False
            use_host_memory = True
            print("use_host_mem="+self.options.use_host_mem)
        dist_homo_graph.load(
            self.options.root_dir,
            self.options.graph_name,
            wm_comm,
            use_chunked,
            use_host_memory,
            wm_embedding_comm,
        )
        print("Rank=%d, Graph loaded." % (comma.Get_rank(),))
        model = self.getModel(dist_homo_graph)
        print("Rank=%d, model created." % (comma.Get_rank(),))
        model.cuda()
        print("Rank=%d, model movded to cuda." % (comma.Get_rank(),))
        model = DDP(model, delay_allreduce=True)
        optimizer = apex.optimizers.FusedAdam(model.parameters(), lr=self.options.lr)
        print("Rank=%d, optimizer created." % (comma.Get_rank(),))

        self.train(train_data, valid_data, model, optimizer)
        self.test(test_data, model)

        wg.finalize_lib()
        print("Rank=%d, wholegraph shutdown." % (comma.Get_rank(),))