#lib TopPPR
import sys
sys.path.append("../../cpp/zhuque_alg/pvr/")
from optparse import OptionParser
import TopPPR

parser = OptionParser()
parser.add_option(
    "-a",
    "--algorithm_name",
    dest="algorithm_name",
    default="PC",
    help="algorithm name, GEN_QUERY|GEN_QUERY_DIS|GEN_GROUND_TRUTH|GEN_GROUND_TRUTH_CENTRALITY|SS|COMPARE|PC|DIS",
)
parser.add_option(
    "--graph_name",
    dest="graph_name",
    default="youtube_u",
    help="input graph name",
)
parser.add_option(
    "-p",
    "--file_path",
    dest="file_path",
    default="../../cpp/zhuque_alg/pvr_c2py/pvr_cpp/dataset/",
    help="input data file path",
)
parser.add_option(
    "--alpha",
    dest="alpha",
    type="float",
    default=0.2,
    help="teleport probability",
)
parser.add_option(
    "--node_count",
    dest="node_count",
    type="int",
    default=20,
    help="number of nodes to run",
)
parser.add_option(
    "-m",
    "--method",
    dest="method_name",
    default="PW",
    help="method_name,PF|PW|PPF|PPW|FORA",
)
parser.add_option(
    "-g",
    "--graphd",
    dest="graph_type",
    default="directed",
    help="graph type,undirected|directed",
)
parser.add_option(
    "--power_iterations",
    dest="power_iterations",
    type="int",
    default=3,
    help="number of power iterations",
)
parser.add_option(
    "-r",
    "--rd_ratio",
    dest="rd_ratio",
    type="float",
    default=1.0,
    help="ratio of random walk to degree",
)
parser.add_option(
    "-s",
    "--samples",
    dest="samples",
    type="float",
    default=0.004,
    help="number of samples",
)
parser.add_option(
    "-e",
    "--epsilon",
    dest="epsilon",
    type="float",
    default=0.5,
    help="error tolerance",
)
parser.add_option(
    "--l1_error",
    dest="l1_error",
    type="float",
    default=1e-5,
    help="l1 error tolerance",
)
parser.add_option(
    "--num_forests",
    dest="num_forests",
    type="int",
    default=6,
    help="number of forests",
)
parser.add_option(
    "-b",
    "--batch_size",
    dest="batch_size",
    type="int",
    default=1,
    help="batch size",
)
(options, args) = parser.parse_args()
TopPPR.algorithm(algo=options.algorithm_name, filename=options.graph_name, datapath=options.file_path,
                 alpha=options.alpha, node_count=options.node_count, method=options.method_name,
                 graphd=options.graph_type,
                 epsilon=options.epsilon)

#test set
# eps_arr=[0.5,0.4,0.3,0.2,0.1]
# for eps in eps_arr:
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",method="MCW",alpha=0.2,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",method="MCF",alpha=0.2,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="FORA",alpha=0.2,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PF",alpha=0.2,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PW",alpha=0.2,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PPF",alpha=0.2,batch_size=2,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PPW",alpha=0.2,batch_size=2,epsilon=eps)
# for eps in eps_arr:
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",method="MCF",alpha=0.2,epsilon=eps,graph="undirected")
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PF",alpha=0.2,epsilon=eps,graph="undirected")
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PPF",alpha=0.2,batch_size=2,epsilon=eps,graph="undirected")
# TopPPR.algorithm(filename="youtube_u",algo="GEN_GROUND_TRUTH",datapath=options.file_path,node_count=50,alpha=0.01)
# TopPPR.algorithm(filename="youtube_u",algo="COMPARE",datapath=options.file_path,alpha=0.01,graph="undirected")
# for eps in eps_arr:
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="SS",node_count=50,method="PF",alpha=0.01,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="SS",node_count=50,method="PW",alpha=0.01,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="SS",node_count=50,method="PPF",alpha=0.01,batch_size=5,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="SS",node_count=50,method="PPW",alpha=0.01,batch_size=5,epsilon=eps)
# for eps in eps_arr:
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="SS",node_count=50,method="PF",alpha=0.01,epsilon=eps,graph="undirected")
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="SS",node_count=50,method="PPF",alpha=0.01,batch_size=5,epsilon=eps,graph="undirected")
# for eps in eps_arr:
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PF",alpha=0.01,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PW",alpha=0.01,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PPF",alpha=0.01,batch_size=5,epsilon=eps)
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PPW",alpha=0.01,batch_size=5,epsilon=eps)
# for eps in eps_arr:
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PF",alpha=0.01,epsilon=eps,graph="undirected")
#     TopPPR.algorithm(filename="youtube_u",datapath=options.file_path,algo="PC",node_count=50,method="PPF",alpha=0.01,batch_size=5,epsilon=eps,graph="undirected")


