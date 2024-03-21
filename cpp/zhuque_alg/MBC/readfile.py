# -*- coding: utf-8 -*-
# import networkx as nx
import random
import json


def readGraph(path, tmin, tmax):
    pathNODE = path + "#" + str(tmin) + "#" + str(tmax) + "#NODE"
    pathDS = path + "#" + str(tmin) + "#" + str(tmax) + "#DS"
    pathadj = path + "#" + str(tmin) + "#" + str(tmax) + "#ADJ"

    file_obj1 = open(pathNODE)
    nodes = []
    # while 1:
        # lines = file_obj1.readlines(10000)
        # if not lines:
        #     break
        # for line in lines:
    for line in file_obj1:
        line_temp = int(line)
        nodes.append(line_temp)
    file_obj1.close()

    file_obj2 = open(pathDS)
    gDS = {}
    i = 0
    for line in file_obj2:
        line_temp = json.loads(line)
        gDS[nodes[i]] = {}
        for item in line_temp:
            item_int = int(item)
            gDS[nodes[i]][item_int] = line_temp[item]
        i += 1
    file_obj2.close()

    file_obj3 = open(pathadj)
    gadj = {}
    i = 0
    for line in file_obj3:
        line_temp = json.loads(line)
        gadj[nodes[i]] = {}
        for item in line_temp:
            item_int = int(item)
            gadj[nodes[i]][item_int] = line_temp[item]
        i += 1
    file_obj3.close()

    maxdegree = 0
    for i in gadj[nodes[-1]]:
        maxdegree += len(gadj[nodes[-1]][i])
    print("max degree:"+str(maxdegree))
    return (gDS, tmin, tmax, gadj)


def format(path):
    tmin = 1000000000000
    tmax = -1
    gDS = {}
    gadj = {}

    file = open(path)
    while 1:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            # for line in open(path):
            if line != '\n':
                line = line.strip('\n').split('\t')
                time_id = int(line[2])
                if time_id > tmax:
                    tmax = time_id
                if time_id < tmin:
                    tmin = time_id
    file.close()

    tmax = tmax - tmin + 1

    file = open(path)
    filename = path + "-" + str(tmin) + "-" + str(tmax)
    file_obj = open(filename, 'w')
    while 1:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            # for line in open(path):
            if line != '\n':
                line = line.strip('\n').split('\t')
                from_id, to_id, time_id = int(line[0]), int(line[1]), int(line[2])
                if from_id == to_id:
                    continue
                time_id -= tmin

                line_temp = [from_id, to_id, time_id]
                json_temp = json.dumps(line_temp)
                file_obj.write(json_temp)
                file_obj.write("\n")


def formatfile(path):
    tmin = 1000000000000
    tmax = -1
    gDS = {}
    gadj = {}

    file = open(path)
    while 1:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            if line != '\n':
                line = line.strip('\n').split('\t')
                time_id = int(line[2])
                if time_id > tmax:
                    tmax = time_id
                if time_id < tmin:
                    tmin = time_id
    file.close()

    tmax = tmax - tmin + 1

    file = open(path)
    while 1:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            if line != '\n':
                line = line.strip('\n').split('\t')
                from_id, to_id, time_id = int(line[0]), int(line[1]), int(line[2])
                if from_id == to_id:
                    continue
                time_id -= tmin

                for (f_id, t_id) in [(from_id, to_id), (to_id, from_id)]:
                    if f_id in gadj:
                        if t_id not in gadj[f_id]:
                            gadj[f_id][t_id] = [time_id]
                        else:
                            if time_id not in gadj[f_id][t_id]:
                                gadj[f_id][t_id].append(time_id)
                    else:
                        gadj[f_id] = {}
                        gadj[f_id][t_id] = [time_id]
    file.close()

    ranktemp = {}
    for node_temp in gadj:
        ranktemp[node_temp] = len(gadj[node_temp])
    rank = sorted(ranktemp.items(), key=lambda item: item[1], reverse=False)
    nodes_by_degree = []
    for item in rank:
        if item[1] > 0:
            nodes_by_degree.append(item[0])

    filename = path + "#" + str(tmin) + "#" + str(tmax) + "#NODE"
    file_obj = open(filename, 'w')
    for line in nodes_by_degree:
        file_obj.write((str(line)))
        file_obj.write("\n")
    file_obj.close()
    filename = path + "#" + str(tmin) + "#" + str(tmax) + "#ADJ"
    file_obj = open(filename, 'w')
    for line in nodes_by_degree:
        file_obj.write(json.dumps(gadj[line]))
        file_obj.write("\n")
    file_obj.close()

    del (gadj)

    file = open(path)
    while 1:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            # for line in open(path):
            if line != '\n':
                line = line.strip('\n').split('\t')
                from_id, to_id, time_id = int(line[0]), int(line[1]), int(line[2])
                if from_id == to_id:
                    continue
                time_id -= tmin

                for (f_id, t_id) in [(from_id, to_id), (to_id, from_id)]:
                    if f_id in gDS:
                        if time_id in gDS[f_id]:
                            if (t_id not in gDS[f_id][time_id]):
                                gDS[f_id][time_id].append(t_id)
                        else:
                            gDS[f_id][time_id] = [t_id]
                    else:
                        gDS[f_id] = {}
                        gDS[f_id][time_id] = [t_id]
    file.close()

    filename = path + "#" + str(tmin) + "#" + str(tmax) + "#DS"
    file_obj = open(filename, 'w')
    for line in nodes_by_degree:
        file_obj.write(json.dumps(gDS[line]))
        file_obj.write("\n")
    file_obj.close()

    return tmin, tmax


if __name__ == '__main__':
    filename = "askubuntu_day"
    # filename = "askubuntu_hour"
    tmin, tmax = formatfile(filename)
    # tmin, tmax = 0, 79
    # (DS, tmin, tmax, adj) = readGraph(filename, tmin, tmax)
