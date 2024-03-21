#lib mbiclique
import sys
sys.path.append("../../cpp/zhuque_alg/BiHSE/")
from optparse import OptionParser
import mbiclique

parser = OptionParser()
parser.add_option(
    "-f",
    "--file_path",
    dest="file_path",
    default="../../cpp/zhuque_alg/BiHSE/datas/fjwiki.txt",
    help="input data file path",
)
parser.add_option(
    "--order",
    type="str",
    dest="order",
    default="two",
    help="order reason , core|two",
)
parser.add_option(
    "-l",
    type="int",
    dest="ls",
    default=3,
)
parser.add_option(
    "-r",
    type="int",
    dest="rs",
    default=3,
)
parser.add_option(
    "--graphMode",
    type="int",
    dest="graphMode",
    default=0,
    help="graph mode noUVM or not, 0|1",
)
(options, args) = parser.parse_args()
mbiclique.run(fPath=options.file_path, order=options.order, ls=options.ls, rs=options.rs,
              noUVM=options.graphMode)
