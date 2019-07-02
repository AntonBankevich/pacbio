import itertools
import os
import shutil
import subprocess
import sys
import time


sys.path.append("py")

from common.dot_parser import DotParser
from disjointig_resolve.unique_marker import UniqueMarker
from alignment.align_tools import Aligner, DirDistributor
from alignment.polishing import Polisher
from common import basic, SeqIO, sam_parser
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.tests import Tester
from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.sequences import ContigCollection, Contig
from common.alignment_storage import ReadCollection, AlignmentPiece
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.saves_io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.cl_params import Params

def CreateLog(dir):
    old_logs_dir = os.path.join(dir, "old")
    basic.ensure_dir_existance(old_logs_dir)
    log_file = os.path.join(dir, "log.info")
    if os.path.isfile(log_file):
        num = len(os.listdir(old_logs_dir))
        shutil.copy(log_file, os.path.join(old_logs_dir, str(num) + ".log"))
    log = open(log_file, "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stdout.prefix = lambda s: time.strftime("%I:%M:%S") + "  "
    sys.stderr = sys.stdout


def main(args):
    cf = args[1]
    rf = args[2]
    dir = args[3]
    CreateLog(dir)
    basic.ensure_dir_existance(dir)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    contigs=ContigCollection().loadFromFasta(open(cf, "r"), False)
    ref = ContigCollection().loadFromFasta(open(rf, "r"), False)
    good = set()
    print "Good"
    for al in aligner.dotplotAlign(contigs, ref):
        if len(al) > 20000:
            al = al.reduce(target=al.seg_to.shrink(10000))
        print al, len(al), al.percentIdentity()
        print "\n".join(al.asMatchingStrings())
        good.add(al.seg_from.contig.id)
    contigs = contigs.filter(lambda contig: contig.id not in good)
    print "Bad"
    for al in aligner.localAlign(contigs, ref):
        print al


if __name__ == "__main__":
    main(sys.argv)



















if __name__ == "__main__":
    main(sys.argv)