#lib biclique
import sys
sys.path.append("../../cpp/zhuque_alg/biclique/")
from optparse import OptionParser
import biclique

parser = OptionParser()
parser.add_option(
    "-f",
    "--file_path",
    dest="file_path",
    default="../../cpp/zhuque_alg/BiHSE/datas/fjwiki.txt",
    help="input data file path",
)
parser.add_option(
    "-o",
    "--out_file_path",
    dest="out_file_path",
    default=".",
    help="output data file path",
)
parser.add_option(
    "-p",
    type="int",
    dest="p",
    default=4,
)
parser.add_option(
    "-q",
    type="int",
    dest="q",
    default=4,
)
parser.add_option(
    "-a",
    "--algorithm_name",
    dest="algorithm_name",
    default="pm",
    help="algorithm name, m|pm|apm|fpmPQ|cp|cppq|peqq|pp|pppeqq|tubcpath|bclist++(none)",
)
parser.add_option(
    "-v",
    "--algorithm_version",
    dest="version",
    default="",
    help="algorithm version, none|v5",
)
parser.add_option(
    "-r",
    dest="r",
    type="float",
    default=0.1,
)
parser.add_option(
    "--bar",
    dest="bar",
    default="1000",
)
parser.add_option(
    "-H",
    type="int",
    dest="H",
    default=11,
)
parser.add_option(
    "-T",
    type="int",
    dest="T",
    default=100000,
)
parser.add_option(
    "--realV",
    dest="realV",
    type="float",
    default=1.0,
)
parser.add_option(
    "--initialDensity",
    dest="initialDensity",
    type="float",
    default=0.0,
)
(options, args) = parser.parse_args()
biclique.run(filePath=options.file_path, outFilePath=options.out_file_path, p=options.p, q=options.q,
             algo=options.algorithm_name, v5=options.version, r=options.r, bar=options.bar,
             H=options.H, T=options.T, realV=options.realV)
#biclique.exdensest(filePath=options.file_path, outFilePath=options.out_file_path, p=options.p, q=options.q,
#                   initialDensity=options.initialDensity)
#biclique.densest(filePath=options.file_path, outFilePath=options.out_file_path, p=options.p, q=options.q)
