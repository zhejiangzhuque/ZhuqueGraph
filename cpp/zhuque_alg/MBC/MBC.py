# -*- coding: utf-8 -*-
# import matplotlib.pyplot as plt
import datetime
import readfile
import gc


class Graph:
    def __init__(self, path):
        path, tmin, tmax = path[0], path[1], path[2]
        self.DSS, self.tmin, self.tmax, self.adj = readfile.readGraph(path, tmin, tmax)

        print("n:" + str(len(self.adj)))
        msum = 0
        m2sum = 0
        for i in self.adj:
            msum += len(self.adj[i])
            for j in self.adj[i]:
                m2sum += len(self.adj[i][j])
        print("m:" + str(msum / 2))
        print("temporal edges:" + str(m2sum / 2))
        print("T:" + str(self.tmax))
        print("Graph loaded!")

    def node_by_degree(self):
        # ranktemp = {}
        # for node_temp in self.adj:
        #     ranktemp[node_temp] = len(self.adj[node_temp])
        # rank = sorted(ranktemp.items(), key=lambda item: item[1], reverse=False)
        #
        # nodes_by_degree = []
        # for item in rank:
        #     if item[1] > 0:
        #         nodes_by_degree.append(item[0])
        # return nodes_by_degree
        return set(self.adj.keys())

    def adj_by_degree(self, nodes_by_degree):
        adj_by_degree = {}
        for node_id in nodes_by_degree:
            # adjtemp = {}
            # for item in self.adj[node_id]:
            #     adjtemp[item] = len(self.adj[item])
            # adjtemp = sorted(adjtemp.items(), key=lambda item: item[1], reverse=False)
            #
            # adj_by_degree[node_id] = []
            # for item in adjtemp:
            #     adj_by_degree[node_id].append(item[0])

            adj_by_degree[node_id] = self.adj[node_id].keys()
        return adj_by_degree

    def core_decomposition(self, nodelist, adj, degeneracy):
        q = []
        deg = {}
        D = set()
        for node_id in nodelist:
            deg[node_id] = len(adj[node_id])

        for k in range(3, len(nodelist)):
            for node_id in nodelist:
                if deg[node_id] < k:
                    q.append(node_id)

            while (len(q) > 0):
                v = q.pop()
                D.add(v)
                for w in adj[v]:
                    if deg[w] >= k:
                        deg[w] = deg[w] - 1
                        if deg[w] < k:
                            q.append(w)

            allowed_nodes = set(nodelist) - D
            if len(allowed_nodes) == 0:
                return k - 1
            for node in allowed_nodes:
                degeneracy[node] = k
        return 0

    def kcore(self, nodelist, adj, k):
        q = []
        deg = {}
        D = []

        for node_id in nodelist:
            deg[node_id] = len(adj[node_id])

        for node_id in nodelist:
            if deg[node_id] < k:
                q.append(node_id)

        while (len(q) > 0):
            v = q.pop()
            D.append(v)
            for w in adj[v]:
                if deg[w] >= k:
                    deg[w] = deg[w] - 1
                    if deg[w] < k:
                        q.append(w)

        adjnew = {}
        Vc = set(nodelist) - set(D)
        for item in Vc:
            adjnew[item] = set(adj[item]) & Vc
        return (Vc, adjnew)

    def MBC(self, Vc, adj, l, delta):
        starttime = datetime.datetime.now()
        (Vc, adj) = self.kcore(Vc, adj, delta)
        endtime = datetime.datetime.now()
        interval = (endtime - starttime)
        print("kcore time:" + str(interval))

        Q = []
        D = set()
        deg = {}
        msd = {}
        DS = {}
        self.initDSS(Vc)

        print("length(Vc):" + str(len(Vc)))

        for node_id in adj:
            deg[node_id] = len(adj[node_id])
            if deg[node_id] < delta:
                msd[node_id] = 0
                Q.append(node_id)
            else:
                msd[node_id] = max(self.computeMSD(l, node_id, 0, DS))
                if msd[node_id] < delta:
                    deg[node_id] = 0
                    Q.append(node_id)

        print("length(first Q):" + str(len(Q)))

        while (len(Q) > 0):
            v = Q.pop()
            D.add(v)
            Vc.remove(v)
            self.maintainDS(v, DS)

            for w in adj[v]:
                if deg[w] >= delta and msd[w] >= delta:
                    deg[w] -= 1
                    if deg[w] < delta:
                        Q.append(w)
                    else:
                        msd[w] = max(self.computeMSD(l, w, Vc, DS))
                        if msd[w] < delta:
                            deg[w] = 0
                            Q.append(w)
        return Vc

    def initDSS(self, Vc):
        del (self.DSS)
        self.DSS = {}
        for from_id in self.adj:
            if from_id in Vc:
                for to_id in self.adj[from_id]:
                    if to_id in Vc:
                        for time_id in self.adj[from_id][to_id]:
                            if from_id in self.DSS:
                                if time_id in self.DSS[from_id]:
                                    if (to_id not in self.DSS[from_id][time_id]):
                                        self.DSS[from_id][time_id].append(to_id)
                                else:
                                    self.DSS[from_id][time_id] = [to_id]
                            else:
                                self.DSS[from_id] = {}
                                self.DSS[from_id][time_id] = [to_id]

    def initDS(self, node_id, Vc, DS):
        if node_id in DS:
            pass
        else:
            if Vc == 0:
                DS[node_id] = [0] * self.tmax
                for t in self.DSS[node_id]:
                    DS[node_id][t] = len(self.DSS[node_id][t])
            else:
                # print("can not happen")
                DS[node_id] = [0] * self.tmax
                for t in self.DSS[node_id]:
                    vnew = Vc & set(self.DSS[node_id][t])
                    DS[node_id][t] = len(vnew)
        return DS

    def maintainDS(self, deleted_node, DS):
        for t in self.DSS[deleted_node]:
            for to_node in self.DSS[deleted_node][t]:
                if to_node in DS:
                    DS[to_node][t] -= 1
        del (DS[deleted_node])
        return DS

    def computeMSD(self, l, node_id, Vc, DS):
        DS = self.initDS(node_id, Vc, DS)
        CSC = [0] * (self.tmax + 1)
        for t in range(self.tmax):
            CSC[t + 1] = CSC[t] + DS[node_id][t]

        if l == self.tmax:
            return [float(CSC[-1]) / l]

        CH = [0] * self.tmax
        i_s = 0
        i_r = -1
        MTS = []

        for t in range(l, len(CSC)):
            while (i_s < i_r) and (self.slop(CH[i_r], t - l, CSC) <= self.slop(CH[i_r - 1], CH[i_r], CSC)):
                i_r -= 1
            i_r += 1
            CH[i_r] = t - l

            while (i_s < i_r) and (self.slop(CH[i_s], t, CSC) >= self.slop(CH[i_s], CH[i_s + 1], CSC)):
                i_s += 1

            if t - CH[i_s] >= 2*l:
                temp = 0
            else:
                temp = self.slop(CH[i_s], t, CSC)
            MTS.append(temp)
        return (MTS)

    def MBC_PLUS(self, Vc, adj, l, delta):
        starttime = datetime.datetime.now()
        (Vc, adj) = self.kcore(Vc, adj, delta)
        endtime = datetime.datetime.now()
        interval = (endtime - starttime)
        print("kCore Time:" + str(interval))
        print("length(Vc):" + str(len(Vc)))

        Q = []
        D = set()
        deg = {}
        MSD = {}
        MTS = {}
        DS = {}
        for node_id in adj:
            deg[node_id] = len(adj[node_id])

        for u in adj:
            if u in D:
                continue
            MTS[u] = self.computeMSD(l, u, Vc, DS)
            MSD[u] = max(MTS[u])
            if (MSD[u] < delta):
                deg[u] = 0
                Q.append(u)
                del (DS[u])

            while (len(Q) > 0):
                v = Q.pop()
                D.add(v)

                Vc.remove(v)
                if v in DS:
                    del (DS[v])
                    del (MTS[v])
                    del (MSD[v])

                for w in adj[v]:
                    if deg[w] >= delta:
                        deg[w] -= 1
                        if deg[w] < delta:
                            Q.append(w)
                            continue

                        if w in MSD:
                            time = []
                            for time_temp in self.adj[w][v]:
                                DS[w][time_temp] -= 1
                                time.append(time_temp)
                            for time_temp in time:
                                self.updateMSD(w, time_temp, l, MTS, DS)
                                MSD[w] = max(MTS[w])
                                if MSD[w] < delta:
                                    Q.append(w)
                                    deg[w] = 0
                                    break

        return Vc

    def updateMSD(self, w, time_temp, l, MTS, DS):

        t_s = max(0, time_temp - 2 * l)
        t_e = min(time_temp + 2 *l , self.tmax)

        CSC = [0] * (t_e - t_s + 1)
        length = t_e - t_s
        for t in range(length):
            CSC[t + 1] = CSC[t] + DS[w][t + t_s]

        CH = [0] * self.tmax
        i_s = 0
        i_r = -1

        for t in range(l, len(CSC)):
            while (i_s < i_r) and (self.slop(CH[i_r], t - l, CSC) <= self.slop(CH[i_r - 1], CH[i_r], CSC)):
                i_r -= 1
            i_r += 1
            CH[i_r] = t - l

            while (i_s < i_r) and (self.slop(CH[i_s], t, CSC) >= self.slop(CH[i_s], CH[i_s + 1], CSC)):
                i_s += 1

            if t - CH[i_s] >= 2*l:
                temp = 0
            else:
                temp = self.slop(CH[i_s], t, CSC)

            if t >= time_temp - t_s:
                MTS[w][t + t_s - l] = temp

    def POMBC(self, C, adj):
        # starttime = datetime.datetime.now()
        # (C, adj) = self.kcore(C, adj, 20)
        # endtime = datetime.datetime.now()
        # interval = (endtime - starttime)
        # print("kCore Time:" + str(interval) + "; 20-core Size:"+ str(len(C)))

        t = self.tmax
        R = []
        MTS = {}
        MSD = {}
        Gc = {}
        deg = {}
        gc.collect()

        l = 2
        for item in adj:
            Gc[item] = set(adj[item])

        while l <= t:
            DS = {}
            for u in C:
                MTS[u] = self.computeMSD(l, u, C, DS)
                MSD[u] = max(MTS[u])
                deg[u] = len(Gc[u])

            (delta, C, Gc) = self.MaxDelta(Gc, l, C, DS, MSD, MTS, deg)
            (l, C) = self.MaxL(Gc, l+1, delta, C)
            R.append([l, delta, C])
            print([l, delta, C])
            k = int(float(delta * l) / (l + 1))
            (C, Gc) = self.kcore(set(adj.keys()), adj, k)
            l += 1

    def MaxDelta(self, adj, l, Vc, DS, MSD, MTS, deg):
        while (1):
            Q = []
            D = set()
            deltalist = sorted(list(MSD.values()))
            delta_min = min(deltalist)
            delta = -1

            for i in deltalist:
                if i > delta_min:
                    delta = i
                    break

            if delta < delta_min:
                return (delta_min, Vc, adj)

            for u in Vc:
                if deg[u] < delta or MSD[u] < delta:
                    Q.append(u)
                    deg[u] = 0
                    MSD[u] = 0
            # delta = int(delta)

            while (len(Q) > 0):
                v = Q.pop()
                D.add(v)
                Vc.remove(v)
                if v in DS:
                    del (DS[v])
                    del (MTS[v])
                    # del (MSD[v])

                for w in adj[v]:
                    if deg[w] >= delta and MSD[w] >= delta:
                        deg[w] -= 1
                        if deg[w] < delta:
                            Q.append(w)
                            continue
                        time = []
                        for time_temp in self.adj[w][v]:
                            DS[w][time_temp] -= 1
                            time.append(time_temp)
                        for time_temp in time:
                            self.updateMSD(w, time_temp, l, MTS, DS)
                            MSD[w] = max(MTS[w])
                            if MSD[w] < delta:
                                Q.append(w)
                                deg[w] = 0
                                break

            if len(Vc) == 0:
                return (delta_min, D, adj)
            else:
                for v in D:
                    del(MSD[v])
                adjtemp = {}
                for node in Vc:
                    adjtemp[node] = set(adj[node]) & Vc
                del(adj)
                adj = adjtemp


    def MaxL(self, adj, l, delta, Vc):
        while (l <= self.tmax):
            Q = []
            D = set()
            MSD = {}
            MTS = {}
            DS = {}
            deg = {}
            for u in adj:
                deg[u] = len(adj[u])

            for u in adj:
                if u in D:
                    continue
                MTS[u] = self.computeMSD(l, u, Vc, DS)
                MSD[u] = max(MTS[u])
                if (MSD[u] < delta):
                    deg[u] = 0
                    Q.append(u)
                    del (DS[u])

                while (len(Q) > 0):
                    v = Q.pop()
                    D.add(v)
                    Vc.remove(v)
                    if v in DS:
                        del (DS[v])
                        del (MTS[v])
                        del (MSD[v])

                    for w in adj[v]:
                        if deg[w] >= delta:
                            deg[w] -= 1
                            if deg[w] < delta:
                                Q.append(w)
                                continue

                            if w in MSD:
                                time = []
                                for time_temp in self.adj[w][v]:
                                    DS[w][time_temp] -= 1
                                    time.append(time_temp)
                                for time_temp in time:
                                    self.updateMSD(w, time_temp, l, MTS, DS)
                                    MSD[w] = max(MTS[w])
                                    if MSD[w] < delta:
                                        Q.append(w)
                                        deg[w] = 0
                                        break

            if len(Vc) == 0:
                return (l - 1, D)
            else:
                if l == self.tmax:
                    return (l, Vc)
                adjtemp = {}
                for node in Vc:
                    adjtemp[node] = set(adj[node]) & Vc
                del(adj)
                adj = adjtemp
                l += 1
        return(l-1, Vc)

    def slop(self, i, j, list):
        result = float(list[j] - list[i]) / (j - i)
        return result

    def AS(self, Vc):
        value1 = 0
        value2 = 0
        nodeset = set()

        for node in Vc:
            for to_node in self.adj[node]:
                if to_node in Vc:
                    value1 += len(self.adj[node][to_node])
                else:
                    nodeset.add(to_node)
                    value2 += len(self.adj[node][to_node])
        return (float(value1)/len(Vc)) / (float(value2)/len(nodeset))

    def AD(self,Vc):
        value1 = 0
        value2 = 0
        for node in Vc:
            for to_node in self.adj[node]:
                if to_node in Vc:
                    value1 += len(self.adj[node][to_node])
                else:
                    value2 += len(self.adj[node][to_node])
        return float(value1)/len(Vc)




if __name__ == '__main__':
    filename = ["00chess_month",0 ,101]
    # filename = ["01lkml_month", -1 ,98]
    # filename = ["02enron_month", 121, 542]
    # filename = ["dblp_year", 0, 79]
    # filename = ["flickr_day",0 ,197]
    # filename = ["youtube_day", 0, 225]

    # filename = ["mathoverflow_day", 0, 2351]
    # filename = ["askubuntu_day", 0, 2614]
    # filename = ["wikitalk_day", 0, 2321]

    # filename = ["mathoverflow_hour", 0, 56409]
    # filename = ["askubuntu_hour", 0, 62732]
    # filename = ["wikitalk_hour", 0, 55690]

    G = Graph(filename)
    nodes = G.node_by_degree()
    adj = G.adj_by_degree(nodes)

    l = 3
    delta = 3

    starttime = datetime.datetime.now()
    Vc = G.MBC(nodes, adj, l, delta)
    endtime = datetime.datetime.now()
    interval = (endtime - starttime)
    print(interval)

    starttime = datetime.datetime.now()
    Vc = G.MBC_PLUS(nodes, adj, l, delta)
    print(len(Vc))
    print(Vc)
    endtime = datetime.datetime.now()
    interval = (endtime - starttime)
    print(interval)
    print(G.AS(Vc))
    print(G.AD(Vc))



