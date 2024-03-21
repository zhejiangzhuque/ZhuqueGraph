#lib fairnessclique
import sys
sys.path.append("../../cpp/zhuque_alg/fairnessclique/")
from optparse import OptionParser
import fairnessclique
parser = OptionParser()
parser.add_option(
    "-p",
    "--file_path",
    type="str",
    dest="file_path",
    default="../../cpp/zhuque_alg/fairnessclique/example.txt",
    help="input data file path",
)
parser.add_option(
    "-a",
    "--attribute_path",
    type="str",
    dest="attribute_path",
    default="../../cpp/zhuque_alg/fairnessclique/example2.txt",
    help="input attribute file path",
)
parser.add_option(
    "--algorithm",
    type="str",
    dest="algorithm",
    default="FairClique",
    help="algorithm to use, FairClique|StrongClique|RelativeStrong|relativeRelativeWeak|Baseline",
)
parser.add_option(
    "--threshold",
    type="str",
    dest="threshold",
    default="2",
    help="threshold for fair cliques ",
)
parser.add_option(
    "--order",
    type="str",
    dest="order",
    default="sorted",
    help="sort reason , [FairClique]sorted|degree|colorful|degeneracy [StrongClique]sorted|fairness|heuristic "
         "[RelativeStrong]sorted|colorful|degeneracy|relative",
)
parser.add_option(
    "--delta",
    type="str",
    dest="delta",
    default="1",
    help="temporal subgraph in which each node has an average degree no less than delta ",
)
(options, args) = parser.parse_args()
fairnessclique.algorithm(filename=options.file_path, attribute=options.attribute_path, algorithm=options.algorithm,
                         threshold=options.threshold, order=options.order, delta=options.delta)