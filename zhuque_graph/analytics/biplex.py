#lib multiPivot
import sys
sys.path.append("../../cpp/zhuque_alg/biplex/")
from google.protobuf.internal.wire_format import INT64_MAX
from optparse import OptionParser
import multiPivot

parser = OptionParser()
parser.add_option(
    "-f",
    "--file_path",
    dest="file_path",
    default="../../cpp/zhuque_alg/biplex/datas/fjwiki.txt",
    help="input data file path",
)
parser.add_option(
    "--order",
    type="str",
    dest="order",
    default="core",
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
    default=1,
    help="graph mode noUVM or not, 0|1",
)
parser.add_option(
    "-k",
    type="int",
    dest="k_plex",
    default=1,
) # l,r >= 2k+1
parser.add_option(
    "--outPutT",
    type="int",
    dest="outPutT",
    default=INT64_MAX,
)
(options, args) = parser.parse_args()
multiPivot.run(fPath=options.file_path, order=options.order, ls=options.ls, rs=options.rs,
               graphMode=options.graphMode, k=options.k_plex, outPutT=options.outPutT)