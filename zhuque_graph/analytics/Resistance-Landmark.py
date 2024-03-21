#lib landmark
import sys
sys.path.append("../../cpp/zhuque_alg/Resistance-Landmark/")
from optparse import OptionParser
import landmark

parser = OptionParser()
parser.add_option(
    "-p",
    "--file_path",
    type="str",
    dest="file_path",
    default="../../cpp/zhuque_alg/Resistance-Landmark/graphs/facebook.txt",
    help="input data file path",
)
parser.add_option(
    "--num_trials",
    dest="num_trials",
    type="int",
    default=1000,
    help="number of trials to run",
)
parser.add_option(
    "-N",
    "--sample_num",
    dest="sample_num",
    type="int",
    default=10000,
    help="number of samples to use ",
)
parser.add_option(
    "--rmax",
    dest="rmax",
    type="float",
    default=1e-4,
    help="maximum radius to use",
)
parser.add_option(
    "-m",
    "--method",
    type="str",
    dest="method",
    default="source_landmark",
    help="method to use,source_landmark|source|monte-carlo|akp|commute|abwalk|localtree|bipush|push",
)
parser.add_option(
    "--in_pairs",
    dest="in_pairs",
    type="str",
    default="node",
    help="whether to use vertex-pairs or node,vertex-pairs|edges|node",
)
parser.add_option(
    "-g",
    "--graph_name",
    type="str",
    dest="graph_name",
    default="facebook",
    help="input graph name",
)
parser.add_option(
    "-l",
    "--landmark",
    dest="landmark",
    type="str",
    default="degree",
    help="landmark to use,degree|core|pagerank|ecc|random",
)
(options, args) = parser.parse_args()
landmark.algorithm(file_name=options.file_path, num_trials=options.num_trials, N=options.sample_num,
                   rmax=options.rmax, method=options.method,
                   in_pairs=options.in_pairs, filename=options.graph_name, landmark=options.landmark)

#tests
#exp facebook
# landmark.algorithm("facebook", "./graphs/facebook.txt", "source_landmark", "node")
# landmark.algorithm("facebook", "./graphs/facebook.txt", "source", "node", 5)
# landmark.algorithm("facebook", "./graphs/facebook.txt", "akp", "vertex-pairs")
# landmark.algorithm("facebook", "./graphs/facebook.txt", "commute", "vertex-pairs")
# landmark.algorithm("facebook", "./graphs/facebook.txt", "abwalk", "vertex-pairs")
# landmark.algorithm("facebook", "./graphs/facebook.txt", "localtree", "vertex-pairs")
# landmark.algorithm("facebook", "./graphs/facebook.txt", "push", "vertex-pairs")
# landmark.algorithm("facebook", "./graphs/facebook.txt", "bipush", "vertex-pairs")
# #exp landmark
# landmark_arr=["degree", "core", "pagerank", "random"]
# for landmark in landmark_arr:
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "abwalk", "vertex-pairs", landmark)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "localtree", "vertex-pairs", landmark)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "bipush", "vertex-pairs", landmark)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "source_landmark", "node", landmark)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "source", "node", landmark)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "push", "vertex-pairs", landmark)
#
# #exp parameter
# rmax_arr=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
# N_arr=[100, 1000, 10000, 100000, 1000000]
# for N in N_arr:
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "abwalk", "vertex-pairs", "degree", N)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "localtree", "vertex-pairs", "degree", N)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "bipush", "vertex-pairs", "degree", N)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "source_landmark", "node", "degree", N)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "source", "node", "degree", N)
# for rmax in rmax_arr:
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "push", "vertex-pairs", "degree", N, rmax)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "bipush", "vertex-pairs", "degree", N, rmax)
#     landmark.algorithm("facebook", "./graphs/facebook.txt", "source-rmax", "node", "degree", N, rmax)