#lib dpcolor
import sys
sys.path.append("../../cpp/zhuque_alg/dpcolor/")
from optparse import OptionParser
import dpcolor
parser = OptionParser()
parser.add_option(
    "--deb",
    type="int",
    dest="deb",
    default=0,
    help="method to use, 0|1|3",
)
parser.add_option(
    "-f",
    "--file_path",
    type="str",
    dest="file_path",
    default="../../cpp/zhuque_alg/dpcolor/data/dblp/",
    help="input data file path",
)
parser.add_option(
    "-N",
    "--samples_num",
    type="int",
    dest="samples_num",
    default=5000000,
    help="number of samples",
)
parser.add_option(
    "--alpha",
    dest="alpha",
    type="float",
    default=1.0,
    help="teleport probability",
)
parser.add_option(
    "-k",
    "--clique_size",
    type="int",
    dest="clique_size",
    default=10,
    help="clique size",
)
parser.add_option(
    "-p",
    "--threads",
    dest="threads",
    type="int",
    default=10,
    help="number of threads",
)
parser.add_option(
    "--algo",
    dest="algo",
    type="str",
    default="",
    help="algorithm to use,cc|ccpath|none",
)
parser.add_option(
    "--f1",
    type="str",
    dest="file_path_CSR",
    default="../../cpp/zhuque_alg/dpcolor/data/dblp/dblp.txt",
    help="input data file for makeCSR",
)
parser.add_option(
    "--pEdgePath",
    type="str",
    dest="pEdgePath",
    default="tmpedge.bin",
    help="pEdgePath for makeCSR and changeToD",
)
parser.add_option(
    "--pIdxPath",
    type="str",
    dest="pIdxPath",
    default="tmpidx.bin",
    help="pIdxPath for makeCSR and changeToD",
)
parser.add_option(
    "--vcnt",
    type="int",
    dest="v_count",
    default=425957,
    help="vertex nums",
)
(options, args) = parser.parse_args()
dpcolor.run(deb=options.deb, filePath=options.file_path, N=options.samples_num, alpha=options.alpha,
            k=options.clique_size, threads=options.threads, algo=options.algo)
dpcolor.makeCSR(f1=options.file_path_CSR, pEdgePath=options.pEdgePath, pIdxPath=options.pIdxPath)
dpcolor.changeToD(pEdgePath=options.pEdgePath, pIdxPath=options.pIdxPath, vCnt_=options.v_count)