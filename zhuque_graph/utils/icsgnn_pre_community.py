import networkx as nx
import numpy as np
import random
import os, sys
import gzip
import pathlib
import  tarfile
p = os.path.dirname(os.path.dirname((os.path.abspath('__file__'))))
if p not in sys.path:
    sys.path.append(p)
import os.path as osp

def pre_com(data_set='com_dblp',subgraph_list=[400], train_ratio=0.02,seed_cnt=20,cmty_size=30):
    path = osp.join(osp.dirname(osp.realpath(__file__)), '..','data',data_set)
    print(f"Load {data_set} edges")
    if(os.path.exists(path + '//edge.npy') == False):
        untar_snap_data(data_set[4:])
    new_edge=np.load(path+'//edges.npy').tolist()
    graph = nx.from_edgelist(new_edge)
    print(f"Load {data_set} cmty")
    com_list=np.load(path+'//comms.npy',allow_pickle=True).tolist()

    com_len=[(i,len(line)) for i,line in enumerate(com_list)]
    com_len.sort(key=lambda x:x[1],reverse=True) 
    for subgraph_size in subgraph_list:
        numlabel = int(subgraph_size * train_ratio / 2) 
        ok_com_len=[(i,lens) for i,lens in com_len if lens>=(numlabel+cmty_size) ]
        seed_list=[]
        train_node=[]
        labels=[]
        error_seed=[]
        time=0
        while len(seed_list)<seed_cnt:
            time+=1
            seed_com_index=random.randint(0,len(ok_com_len)-1)
            seed_com=com_list[ok_com_len[seed_com_index][0]]
            seed_com_suff=seed_com[:]
            random.shuffle(seed_com_suff)
            seed_index=0
            seed=seed_com_suff[seed_index]
            while (seed in seed_list or seed in error_seed ) and (seed_index+1)<len(seed_com_suff):
                seed_index+=1
                seed=seed_com_suff[seed_index]
            if(seed in seed_list or seed in error_seed ): 
                continue
            allNodes=[] 
            allNodes.append(seed)
            pos = 0
            while pos < len(allNodes) and pos < subgraph_size and len(allNodes) < subgraph_size:
                cnode = allNodes[pos]
                for nb in graph.neighbors(cnode):
                    if nb not in allNodes and len(allNodes) < subgraph_size:
                        allNodes.append(nb)
                pos += 1
            posNodes = []
            posNodes.append(seed)
            seed_com_intersection=list(set(seed_com).intersection(set(allNodes)))
            if(len(seed_com_intersection)< numlabel+cmty_size):
                error_seed.append(seed)
                continue
            seed_com_intersection_noseed=seed_com_intersection[:]
            seed_com_intersection_noseed.remove(seed)
            random.shuffle(seed_com_intersection_noseed)
            posNodes.extend(seed_com_intersection_noseed[:numlabel-1])
            negNodes=list(set(allNodes).difference(set(seed_com)))
            if(len(negNodes)< numlabel):
                error_seed.append(seed)
                continue
            random.shuffle(negNodes)
            negNodes=negNodes[:numlabel]
            seed_list.append(seed)
            train_node.append(posNodes+negNodes)
            labels.append(seed_com_intersection)
        print('error:',len(error_seed),"seed_list:",seed_list)
    return new_edge,seed_list,train_node,labels

def load_facebook(seed):
    path = osp.join(osp.dirname(osp.realpath(__file__)), '..', 'data','facebook')
    print('Load facebook data')
    if(os.path.exists(path + f'//{str(seed)}.circles') == False):
        untart_facebook()
    file_circle= path + f'//{str(seed)}.circles'
    file_edges=path + f'//{str(seed)}.edges'
    file_egofeat=path + f'//{str(seed)}.egofeat'
    file_feat=path + f'//{str(seed)}.feat'
    edges=[]
    node=[]
    feature = {}
    with open(file_egofeat) as f:
        feature[seed] = [int(i) for i in f.readline().split()]
    with open(file_feat) as f:
        for line in f:
            line = [int(i) for i in line.split()]
            feature[int(line[0])] = line[1:]
            node.append(int(line[0]))
    with open(file_edges,'r') as f:
        for line in f:
            u,v=line.split()
            u=int(u)
            v=int(v)
            if(u in feature.keys() and v in feature.keys()):
                edges.append((u,v))

    for i in node:
        edges.append((seed, i))
    node=sorted(node+[seed])
    mapper = {n: i for i, n in enumerate(node)}
    edges=[(mapper[u],mapper[v]) for u,v in edges]
    node=[mapper[u] for u in node]

    features=[0]*len(node)
    for i in list(feature.keys()):
        features[mapper[i]]=feature[i]
    circle=[]
    with open(file_circle) as f:
        for line in f:
            line=line.split()
            line=[ mapper[int(i)] for i  in line[1:]]
            if(len(line)<8):continue
            circle.append(line)

    seed=mapper[seed]

    return edges,features,circle,seed

def load_snap(data_set,com_size):
    print(f'Load {data_set} edge')
    path = osp.join(osp.dirname(osp.realpath(__file__)), '..', 'data', data_set)
    if(os.path.exists(path + '//edge.npy') == False):
        untar_snap_data(data_set[4:])
    edges=np.load(path + '//edge.npy').tolist()
    print(f'Load {data_set} cmty')
    com_list = np.load(path + '//comms.npy', allow_pickle=True).tolist()
    com_list=[i for i in com_list if len(i)>=com_size]
    return edges,com_list

def untar_snap_data(name):
    """Load the snap comm datasets."""
    print(f'Untar {name} edge')
    root = pathlib.Path('raw')
    #print(root)
    with gzip.open(root / f'com-{name}.ungraph.txt.gz', 'rt') as fh:
        edges = fh.read().strip().split('\n')[4:]
    edges = [[int(i) for i in e.split()] for e in edges]
    nodes = {i for x in edges for i in x}
    mapping = {u: i for i, u in enumerate(sorted(nodes))}
    edges = [[mapping[u], mapping[v]] for u, v in edges]
    print(f'Untar {name} cmty')
    with gzip.open(root / f'com-{name}.top5000.cmty.txt.gz', 'rt') as fh:
        comms = fh.readlines()
    comms = [[mapping[int(i)] for i in x.split()] for x in comms]
    root = pathlib.Path()/'data'/f'com_{name}'
    root.mkdir(exist_ok=True, parents=True)
    np.save(root/'edges',edges)
    np.save(root/'comms',comms,allow_pickle=True)
    np.save(root/'map',mapping,allow_pickle=True)

def untart_facebook():
    print(f'Untar  facebook')
    tar = tarfile.open(osp.join(osp.dirname(osp.realpath(__file__)), '..', 'raw','facebook.tar.gz'))
    names = tar.getnames()
    path =osp.join(osp.dirname(osp.realpath(__file__)), '..', 'data')
    for name in names:
        tar.extract(name,path)
    tar.close()