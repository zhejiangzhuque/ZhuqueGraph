import torch
import random
import numpy as np
import networkx as nx
from zhuque_graph.nn.pytorch.model.ClusterGcnModel import ClusterGCNTrainer
from zhuque_graph.utils.icsgnn_community import LocalCommunity
import time
import datetime


class SubGraph(object):
    def __init__(self, args, graph, features, target):
        self.args = args
        self.graph = graph
        self.features = features
        self.target = target
        self._set_sizes()
        self.methods = {}
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.time_map = {}
        self.rankloss = 0
        self.posforrank = []
        self.negforrank = []
        self.cntmap = {}

    def _set_sizes(self):
        self.feature_count = self.features.shape[1]
        self.class_count = 2
        self.clusters = [0]
        self.cluster_membership = {node: 0 for node in self.graph.nodes()}

    def build_local_candidate_iteration(self, seed):
        '''
        Build subgraphs with iteration
        '''
        posNodes = []
        negNodes = []
        allNodes = []
        length = self.args.subgraph_size
        if (len(self.allnode) == 1):
            allNodes.append(seed)
            numLabel = int(length * self.args.train_ratio / 2)
            pos = 0
            while pos < len(allNodes) and pos < length and len(allNodes) < length:
                cnode = allNodes[pos]
                for nb in self.graph.neighbors(cnode):
                    if nb not in allNodes and len(allNodes) < length:
                        allNodes.append(nb)
                pos = pos + 1
        else:
            numLabel = self.args.possize
            allNodes = self.allnode[:]
        print("The length of list is %d" % len(allNodes))
        print("The degree of seed is %d" % self.graph.degree(seed))
        for i in self.oldpos:
            if i not in allNodes:
                allNodes.append(i)
                print("pos not in subgraph")
            cnt = self.args.upsize
            for nb in self.graph.neighbors(i):
                if (cnt == 0): break
                if (nb not in allNodes):
                    allNodes.append(nb)
                    cnt -= 1
        self.allnode = allNodes
        seedLabel = self.target[seed]
        for node in allNodes:
            if (node == seed):
                continue
            if self.target[node] == seedLabel and node not in self.oldpos:
                posNodes.append(node)
            elif self.target[node] != seedLabel and node not in self.oldneg:
                negNodes.append(node)

        random.shuffle(posNodes)
        random.shuffle(negNodes)
        print('extern pos size', len(posNodes))
        print('extern neg size', len(negNodes))
        posNodes = posNodes[:numLabel]
        negNodes = negNodes[:numLabel]
        if (len(posNodes) < numLabel):
            print('e1')
            return 0
        if (seed not in self.oldpos):
            posNodes[0] = seed
        if (len(negNodes) < numLabel):
            print('e2')
            return 0

        posNodes = self.oldpos + posNodes
        negNodes = self.oldneg + negNodes
        print("Positive Nodes are ")
        print(posNodes)
        print("Negative Nodes are ")
        print(negNodes)
        self.oldpos = posNodes[:]
        self.oldneg = negNodes[:]
        negNodes = []
        for i in self.oldneg:
            if (i in allNodes):
                negNodes.append(i)
            else:
                print('neg不在里面')


        self.sg_nodes = {}
        self.sg_edges = {}
        self.sg_train_nodes = {}
        self.sg_test_nodes = {}
        self.sg_features = {}
        self.sg_targets = {}

        self.subgraph = nx.Graph()
        for i in range(len(allNodes)):
            for j in range(i):
                if ((allNodes[i], allNodes[j]) in self.graph.edges) or ((allNodes[j], allNodes[i]) in self.graph.edges):
                    self.subgraph.add_edge(allNodes[i], allNodes[j])

        print("size of nodes %d size of edges %d" % (len(self.subgraph.nodes), len(self.subgraph.edges)))
        self.sg_nodes[0] = [node for node in sorted(self.subgraph.nodes())]
        self.sg_predProbs = [0.0] * len(self.sg_nodes[0])
        self.sg_predLabels = [0] * len(self.sg_nodes[0])
        self.mapper = {node: i for i, node in enumerate(sorted(self.sg_nodes[0]))}
        self.rmapper = {i: node for i, node in enumerate(sorted(self.sg_nodes[0]))}
        self.sg_edges[0] = [[self.mapper[edge[0]], self.mapper[edge[1]]] for edge in self.subgraph.edges()] + [
            [self.mapper[edge[1]], self.mapper[edge[0]]] for edge in self.subgraph.edges()]
        self.sg_posNodes = [self.mapper[node] for node in posNodes]
        self.sg_negNodes = [self.mapper[node] for node in negNodes]
        allNodes1 = [self.mapper[node] for node in allNodes]
        self.sg_train_nodes[0] = self.sg_posNodes + self.sg_negNodes
        self.sg_test_nodes[0] = list(set(allNodes1).difference(set(self.sg_train_nodes[0])))
        self.sg_test_nodes[0] = sorted(self.sg_test_nodes[0])
        self.sg_train_nodes[0] = sorted(self.sg_train_nodes[0])
        self.sg_features[0] = self.features[self.sg_nodes[0], :]
        self.sg_targets[0] = self.target[self.sg_nodes[0], :]
        self.sg_targets[0] = self.sg_targets[0] == seedLabel
        self.sg_targets[0] = self.sg_targets[0].astype(int)


        print("Value 0 %d, Value 1 %d" % (sum(self.sg_targets[0] == 0), sum(self.sg_targets[0] == 1)))
        for x in self.sg_posNodes:
            self.sg_predProbs[x] = 1.0
            self.sg_predLabels[x] = 1
            if self.sg_targets[0][x] != 1.0:
                print("wrong1")
        for x in self.sg_negNodes:
            self.sg_predProbs[x] = 0.0
            self.sg_predLabels[x] = 0
            if self.sg_targets[0][x] != 0:
                print("wrong0")

        self.transfer_edges_and_nodes()
        self.TOPK_SIZE = int(self.args.community_size)
    def build_local_candidate(self, seed, trian_node, label):
        '''
        Build subgraphs
        '''
        allNodes = []
        allNodes.append(seed)
        posNodes = set()
        negNodes = set()

        length = self.args.subgraph_size
        numLabel = int(length * self.args.train_ratio / 2)
        pos = 0
        while pos < len(allNodes) and pos < length and len(allNodes) < length:
            cnode = allNodes[pos]
            for nb in self.graph.neighbors(cnode):
                if nb not in allNodes and len(allNodes) < length:
                    allNodes.append(nb)
                    if(nb!=seed and self.target is not None):
                        if(self.target[nb] == self.target[seed]):
                            posNodes.add(nb)
                        else:
                            negNodes.add(nb)
            pos = pos + 1
        posNodes=list(posNodes)
        negNodes=list(negNodes)
        print("The length of list is %d" % len(allNodes))
        print("The degree of seed is %d" % self.graph.degree(seed))
        if (trian_node is not None):
            posNodes = trian_node[:numLabel]
            negNodes = trian_node[numLabel:]
        else:
            seedLabel = self.target[seed]
            posNodes.append(seed)
            if(len(posNodes+[seed])<numLabel or len(negNodes)<numLabel):
                return 0
            random.shuffle(posNodes)
            random.shuffle(negNodes)
            posNodes=[seed]+posNodes[:numLabel-1]
            negNodes=negNodes[:numLabel]
        print("Positive Nodes are ")
        print(posNodes)
        print("Negative Nodes are ")
        print(negNodes)


        self.sg_nodes = {}
        self.sg_edges = {}
        self.sg_train_nodes = {}
        self.sg_test_nodes = {}
        self.sg_features = {}
        self.sg_targets = {}

        self.subgraph = nx.Graph()
        for i in range(len(allNodes)):
            for j in range(i):
                if ((allNodes[i], allNodes[j]) in self.graph.edges) or ((allNodes[j], allNodes[i]) in self.graph.edges):
                    self.subgraph.add_edge(allNodes[i], allNodes[j])

        print("size of nodes %d size of edges %d" % (len(self.subgraph.nodes), len(self.subgraph.edges)))
        self.sg_nodes[0] = [node for node in sorted(self.subgraph.nodes())]
        self.sg_predProbs = [0.0] * len(self.sg_nodes[0])
        self.sg_predLabels = [0] * len(self.sg_nodes[0])
        self.mapper = {node: i for i, node in enumerate(sorted(self.sg_nodes[0]))}
        self.rmapper = {i: node for i, node in enumerate(sorted(self.sg_nodes[0]))}
        self.sg_edges[0] = [[self.mapper[edge[0]], self.mapper[edge[1]]] for edge in self.subgraph.edges()] + [
            [self.mapper[edge[1]], self.mapper[edge[0]]] for edge in self.subgraph.edges()]
        self.sg_posNodes = [self.mapper[node] for node in posNodes]
        self.sg_negNodes = [self.mapper[node] for node in negNodes]

        allNodes1 = [self.mapper[node] for node in allNodes]
        self.sg_train_nodes[0] = self.sg_posNodes + self.sg_negNodes
        self.sg_test_nodes[0] = list(set(allNodes1).difference(set(self.sg_train_nodes[0])))

        self.sg_test_nodes[0] = sorted(self.sg_test_nodes[0])
        self.sg_train_nodes[0] = sorted(self.sg_train_nodes[0])

        self.sg_features[0] = self.features[self.sg_nodes[0], :]
        if (label is None):
            self.sg_targets[0] = self.target[self.sg_nodes[0], :]
            self.sg_targets[0] = self.sg_targets[0] == seedLabel
            self.sg_targets[0] = self.sg_targets[0].astype(int)
        else:
            self.sg_targets[0] = [0] * len(allNodes1)
            for i in label:
                self.sg_targets[0][self.mapper[i]] = 1
            self.sg_targets[0] = np.array(self.sg_targets[0])
            self.sg_targets[0] = self.sg_targets[0][:, np.newaxis]
            self.sg_targets[0] = self.sg_targets[0].astype(int)

        print("Value 0 %d, Value 1 %d" % (sum(self.sg_targets[0] == 0), sum(self.sg_targets[0] == 1)))
        for x in self.sg_posNodes:
            self.sg_predProbs[x] = 1.0
            self.sg_predLabels[x] = 1
            if self.sg_targets[0][x] != 1.0:
                print("wrong")
        for x in self.sg_negNodes:
            self.sg_predProbs[x] = 0.0
            self.sg_predLabels[x] = 0
            if self.sg_targets[0][x] != 0:
                print("wrong")

        self.transfer_edges_and_nodes()
        self.TOPK_SIZE = self.args.community_size


    def transfer_edges_and_nodes(self):
        '''
        Transfering the data to PyTorch format.
        '''
        for cluster in self.clusters:
            self.sg_nodes[cluster] = torch.LongTensor(self.sg_nodes[cluster]).to(self.device)
            self.sg_edges[cluster] = torch.LongTensor(self.sg_edges[cluster]).t().to(self.device)
            self.sg_train_nodes[cluster] = torch.LongTensor(self.sg_train_nodes[cluster]).to(self.device)
            self.sg_test_nodes[cluster] = torch.LongTensor(self.sg_test_nodes[cluster]).to(self.device)
            self.sg_features[cluster] = torch.FloatTensor(self.sg_features[cluster]).to(self.device)
            self.sg_targets[cluster] = torch.LongTensor(self.sg_targets[cluster]).to(self.device)


    def community_search(self, seed, trian_node, label):
        '''
        GNN training subgraph, heuristic search community without/with rking loss
        '''
        self.rankloss = 0
        isOK = self.build_local_candidate(seed, trian_node, label)
        if isOK == 0:
            print("cannot build a local subgraph")
            return 0
        for round in range(2):
            keepLayers = self.args.layers.copy()
            gcn_trainer = ClusterGCNTrainer(self.args, self)
            begin_time = time.time()
            nodeweight, predlabels, f1score = gcn_trainer.train_test_community()
            if 'gcn' not in self.time_map:
                self.time_map['gcn'] = time.time() - begin_time
                self.cntmap['gcn']=1
            else:
                self.time_map['gcn'] = time.time() - begin_time + self.time_map['gcn']
                self.cntmap['gcn'] +=1
            self.args.layers = keepLayers
            lc = LocalCommunity(self.args, self)

            for i in range(len(self.sg_test_nodes[0])):
                self.sg_predProbs[self.sg_test_nodes[0][i]] = nodeweight[i].item()
                self.sg_predLabels[self.sg_test_nodes[0][i]] = predlabels[i].item()

            if(self.rankloss == 1):
                prefix='With rking loss'
            else:
                prefix="Without rking loss"


            begin_time = time.time()
            topk = lc.locate_community_BFS_only(seed)
            lc.evaluate_community(topk, prefix + " BSF Only",time.time() - begin_time)


            begin_time = time.time()
            topk = lc.locate_community_BFS(seed)
            lc.evaluate_community(topk, prefix + " BSF Swap", time.time() - begin_time)




            begin_time = time.time()
            topk = lc.locate_community_greedy(seed)
            lc.evaluate_community(topk, prefix + " Greedy-T", time.time() - begin_time)



            begin_time = time.time()
            topk = lc.locate_community_greedy_graph_prepath(seed)
            lc.evaluate_community(topk, prefix + " Greedy-G", time.time() - begin_time)

            self.posforrank, self.negforrank = self.getPNpairs()
            self.rankloss = 1


        return 1


    def getPNpairs(self):
        '''
        Get rking loss pair
        '''
        probs = self.sg_predProbs.copy()
        for x in self.sg_train_nodes[0]:
            probs[x] = 2
        for i in range(len(self.sg_targets[0])):
            if self.sg_targets[0][i] == 0:
                probs[i] = 2
        posIdx = np.argsort(np.array(probs))[0:int(self.args.train_ratio * self.args.subgraph_size / 2)]
        probs = self.sg_predProbs.copy()
        for x in self.sg_train_nodes[0]:
            probs[x] = -2
        for i in range(len(self.sg_targets[0])):
            if self.sg_targets[0][i] == 1:
                probs[i] = -2
        negIdx = np.argsort(-np.array(probs))[0:int(self.args.train_ratio * self.args.subgraph_size / 2)]
        return posIdx, negIdx

    def community_search_iteration(self,seed):
        '''
        GNN training subgraph, heuristic search community  with iteration without rking loss
        '''
        self.seed = seed
        self.oldpos = []
        self.oldneg = []
        self.allnode = [seed]
        for round in range(self.args.round):
            seed = self.seed
            isOK = self.build_local_candidate_iteration(seed)
            if isOK == 0:
                print("cannot build a local subgraph")
                return 0
            keepLayers = self.args.layers.copy()
            gcn_trainer = ClusterGCNTrainer(self.args, self)
            begin_time = time.time()
            nodeweight, predlabels, f1score = gcn_trainer.train_test_community()
            if 'gcn' not in self.time_map:
                self.time_map['gcn'] = time.time() - begin_time
                self.cntmap['gcn']=1
            else:
                self.time_map['gcn'] = time.time() - begin_time + self.time_map['gcn']
                self.cntmap['gcn'] +=1
            self.args.layers = keepLayers
            lc = LocalCommunity(self.args, self)

            for i in range(len(self.sg_test_nodes[0])):
                self.sg_predProbs[self.sg_test_nodes[0][i]] = nodeweight[i].item()
                self.sg_predLabels[self.sg_test_nodes[0][i]] = predlabels[i].item()

            prefix=str(round)+' Round'
            begin_time = time.time()
            topk = lc.locate_community_BFS_only(seed)
            lc.evaluate_community(topk, prefix + " BSF Only",time.time() - begin_time)


            begin_time = time.time()
            topk = lc.locate_community_BFS(seed)
            lc.evaluate_community(topk, prefix + " BSF Swap", time.time() - begin_time)




            begin_time = time.time()
            topk = lc.locate_community_greedy(seed)
            lc.evaluate_community(topk, prefix + " Greedy-T", time.time() - begin_time)



            begin_time = time.time()
            topk = lc.locate_community_greedy_graph_prepath(seed)
            lc.evaluate_community(topk, prefix + " Greedy-G", time.time() - begin_time)

        return 1
    def methods_result(self):
        '''
        save result
        '''
        file_handle = open('results.txt', mode='a+')
        now = datetime.datetime.now()
        sTime = now.strftime("%Y-%m-%d %H:%M:%S")
        file_handle.write(sTime + "\n")
        args = vars(self.args)
        keys = sorted(args.keys())
        keycontent = [[k.replace("_", " ").capitalize(), args[k]] for k in keys]
        for x in keycontent:
            file_handle.writelines(str(x) + "\n")
        file_handle.write("gcn time %f \n" % (self.time_map['gcn'] / self.cntmap['gcn']))
        print("gcn time %f " % (self.time_map['gcn'] / self.cntmap['gcn']))
        for method in self.methods:
            if isinstance(self.methods[method], list):
                pre = self.methods[method][0] / self.cntmap[method]
                rpre= self.methods[method][1] / self.cntmap[method]
                times = self.time_map[method] / self.cntmap[method]
                print(
                    "%s Method achieve the average precision= %f precision without posnode = %f using %d seeds with avgtime=%f s " % (
                        method, pre, rpre,self.cntmap[method], times))
                file_handle.writelines(
                    "%s Method achieve the average precision= %f precision without posnode = %f using %d seeds with time=%f s\n" % (
                        method, pre, rpre,self.cntmap[method], times))

        file_handle.writelines("\n")
        file_handle.close()
