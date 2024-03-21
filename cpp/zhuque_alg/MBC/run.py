import click
import MBC
import datetime
import readfile

@click.command()
@click.option('--name', prompt='Type one number to chose the dataset: [1]Chess; [2]Lkml; [3]Enron. (int)', help='The name of the loaded dataset')
@click.option('--l', prompt='l(int)', help='The value of the parameter l')
@click.option('--delta', prompt='Delta(int)', help='The value of the parameter Delta')
@click.option('--method', prompt='Type one number to chose the algorithm: [1]MBC; [2]MBCPLUS; [3]POMBC. (int)', help='Three bursting community detection algorithms')
def doit(name, l, delta, method):
    if name is "1":
        filename = ["00chess_month",0 ,101]
    if name is "2":
        filename = ["01lkml_month", -1 ,98]
    if name is "3":
        filename = ["02enron_month", 0, 87]

    l = int(l)
    readfile.formatfile(filename[0])
    delta = float(delta)
    G = MBC.Graph(filename)
    nodes = G.node_by_degree()
    adj = G.adj_by_degree(nodes)

    if method is "1":
        starttime = datetime.datetime.now()
        Vc = G.MBC(nodes, adj, l, delta)
        print(Vc)
        endtime = datetime.datetime.now()
        interval = (endtime - starttime)
        print(interval)

    if method is "2":
        starttime = datetime.datetime.now()
        Vc = G.MBC_PLUS(nodes, adj, l, delta)
        print(Vc)
        endtime = datetime.datetime.now()
        interval = (endtime - starttime)
        print(interval)

    if method is "3":
        starttime = datetime.datetime.now()
        G.POMBC(nodes, adj)
        endtime = datetime.datetime.now()
        interval = (endtime - starttime)
        print(interval)

if __name__ == '__main__':
    doit()