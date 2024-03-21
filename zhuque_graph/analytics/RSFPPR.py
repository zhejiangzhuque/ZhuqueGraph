#lib rsf
import sys
sys.path.append("../../cpp/zhuque_alg/RSFPPR/")
from optparse import OptionParser
import rsf

parser = OptionParser()
parser.add_option(
    "-a",
    "--algorithm_name",
    dest="algorithm_name",
    default="GROUND_TRUTH",
    help="algorithm name, GEN_QUERY|GROUND_TRUTH|GROUND_TRUTH_BACK|ITERATIVE_METHOD|RANDOM_WALK|FORA|SPEEDPPR|BACK|RBACK|COMPARE_RESULTS_BACK|COMPARE_RESULTS|VERY_SMALL_ALPHA|PLOT_RESULTS",
)
parser.add_option(
    "--graph_name",
    dest="graph_name",
    default="dblp_weighted",
    help="input graph name",
)
parser.add_option(
    "-p",
    "--file_path",
    dest="file_path",
    default="../../cpp/zhuque_alg/RSFPPR/data/",
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
    dest="num_query_nodes",
    type="int",
    default=20,
    help="number of query nodes to run",
)
parser.add_option(
    "--l1_error",
    dest="l1_error",
    type="float",
    default=1e-12,
    help="l1 error tolerance",
)
parser.add_option(
    "--epsilon",
    dest="epsilon",
    type="float",
    default=0.1,
    help="error tolerance for RSFPPR",
)
parser.add_option(
    "--r_epsilon",
    dest="r_epsilon",
    type="float",
    default=1e-12,
    help="error tolerance for RBACK",
)
(options, args) = parser.parse_args()
rsf.algorithm(algo=options.algorithm_name, filename=options.graph_name,
              datapath=options.file_path, num_query_nodes=options.num_query_nodes, alpha=options.alpha,
              l1_error=options.l1_error, eps=options.epsilon)


