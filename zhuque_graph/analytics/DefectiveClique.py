#lib clique
import sys
sys.path.append("../../cpp/zhuque_alg/DefectiveClique/")
from _testcapi import LONG_MAX
from optparse import OptionParser
import clique

parser = OptionParser()
parser.add_option(
    "-w",
    "--work",
    dest="work",
    default="branch",
    help="branch or delay",
)
parser.add_option(
    "-d",
    "--data_file",
    dest="data_file",
    default="../../cpp/zhuque_alg/DefectiveClique/src/datas/ca-GrQc.txt",
    help="input data file",
)
parser.add_option(
    "-q",
    "--q",
    dest="threads_num",
    type="int",
    default=10,
    help="number of threads,[2,n]",
)
parser.add_option(
    "-k",
    "--k",
    dest="clique_size",
    type="int",
    default=1,
    help="size of the clique,[1,5]",
)
parser.add_option(
    "-a",
    "--a",
    dest="iterations_num",
    type="int",
    default=6,
    help="number of iterations,[1,7]",
)
parser.add_option(
    "-r",
    "--r",
    dest="repetitions_num",
    type="long",
    default=LONG_MAX,
    help="number of repetitions,default LONG_MAX",
)
(options, args) = parser.parse_args()
clique.algorithm(work=options.work,data=options.data_file,q=options.threads_num,k=options.clique_size,
                 a=options.iterations_num,r=options.repetitions_num)

