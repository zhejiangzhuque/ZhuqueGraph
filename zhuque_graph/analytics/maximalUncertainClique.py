#lib maximalCliques
import sys
sys.path.append("../../cpp/zhuque_alg/maximalUncertainClique/")
from optparse import OptionParser
import maximalCliques

parser = OptionParser()
parser.add_option(
    "-p",
    "--file_path",
    type="str",
    dest="file_path",
    default="../../cpp/zhuque_alg/maximalUncertainClique/example/example.txt",
    help="input data file path",
)
parser.add_option(
    "-a",
    "--alg",
    dest="alg",
    type="int",
    default=2,
    help="algorithm to use, 1 topKEtaCore|2 topKEtatriangle",
)
parser.add_option(
    "-k",
    "--min_clique_size",
    dest="min_clique_size",
    type="int",
    default=3,
    help="minsize constraint, [0,n]",
)
parser.add_option(
    "-e",
    "--eta",
    dest="eta",
    type="float",
    default=0.5,
    help="clique probability constraint, [0,1]",
)
(options, args) = parser.parse_args()
maximalCliques.algorithm(filename=options.file_path, alg=options.alg,
                         k=options.min_clique_size, eta=options.eta)
