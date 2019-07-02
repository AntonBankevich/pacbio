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
    rf = args[2]
    dir = args[3]
    CreateLog(dir)
    disjointigs = ContigCollection().loadFromFasta(open(args[1], "r"), num_names=False)
    basic.ensure_dir_existance(dir)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    clen = 5000000
    reads = ReadCollection().loadFromFasta(open(rf, "r"))
    tlen0 = sum(map(len, reads))
    for i in range(10):
        good_reads = set()
        for al in aligner.localAlign(reads, disjointigs):
            if not al.contradictingRTC(al.seg_to.contig.asSegment(), 500):
                good_reads.add(al.seg_from.contig.id)
        rf = os.path.join(dir, "reads" + str(i) + ".fasta")
        reads = reads.filter(lambda read: read.id not in good_reads).cleanCopy()
        tlen = sum(map(len, reads))
        reads.print_fasta(open(rf, "w"))
        l = tlen * clen / tlen0
        assembly_dir = os.path.join(dir, "assembly" + str(i))
        subprocess.check_call(["./bin/flye", "-o", assembly_dir, "-t", "8", "--pacbio-raw", rf, "--genome-size", str(l), "--no-trestle"])
        df= os.path.join(assembly_dir, "10-consensus", "consensus.fasta")
        disjointigs.addAll(ContigCollection().loadFromFasta(open(df, "r"), num_names=False))
        df = os.path.join(dir, "df" + str(i) + ".fasta")
        disjointigs.print_fasta(open(df, "w"))
        # TODO: add new disjointigs

if __name__ == "__main__":
    main(sys.argv)



















if __name__ == "__main__":
    main(sys.argv)