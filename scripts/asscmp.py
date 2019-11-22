import itertools
import os
import sys

from alignment.align_tools import Aligner, DirDistributor
from common.basic import CreateLog
from common.sequences import ContigStorage


def main(args):
    flye_dir = sys.argv[1]
    contigsf = sys.argv[2]
    reff = sys.argv[3]
    dir = sys.argv[4]
    CreateLog(dir)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    print "Reading contigs and reference"
    contigs = ContigStorage().loadFromFasta(open(contigsf, "r"), False)
    ref = ContigStorage().loadFromFasta(open(reff, "r"), False)
    print "Aligning reads"
    filtered = filter(lambda al: len(al) > 50000, aligner.referenceAlign(contigs, ref))
    print "Sorting reads along reference"
    als = sorted(filtered, key = lambda al: al.seg_to.contig.id)
    filtered = []
    for contig, it in itertools.groupby(als, lambda al: al.seg_to.contig):
        print contig
        cals = sorted(it, key = lambda al: -len(al))
        clen = 0
        tmp = []
        for al in cals:
            if clen  + 10000 > len(contig):
                return
            filtered.append(al)
            tmp.append(al)
        tmp = sorted(tmp, key = lambda al: al.seg_to.left)
        print tmp

    print "Sorting reference along reads"
    als = sorted(filtered, key = lambda al: al.seg_from.contig.id)
    for contig, it in itertools.groupby(als, lambda al: al.seg_from.contig):
        print contig
        cals = sorted(it, key = lambda al: al.seg_from.left)
        print cals
