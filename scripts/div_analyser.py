import os
import shutil
import subprocess
import sys

sys.path.append("py")

from common.seq_records import NamedSequence
from common.SimpleGraph import SimpleGraph, Edge
from common.sequences import ContigStorage, Contig
from alignment.align_tools import DirDistributor, Aligner
from common import basic, params, SeqIO


def main(contigs_file, contig_name, reads_file, dir, k):
    basic.ensure_dir_existance(dir)
    basic.CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    contigs = ContigStorage().loadFromFasta(open(contigs_file, "r"), False)
    contig = contigs[contig_name]
    contigs = ContigStorage()
    contigs.add(contig)
    reads = ContigStorage().loadFromFasta(open(reads_file, "r"), False)
    als = list(aligner.localAlign(reads, contigs))
    tmp = []
    for al in als:
        if al.seg_to.contig != contig:
            al = al.rc
        tmp.append(al)
    als = tmp
    als = sorted(als, key = lambda al: al.seg_to.left)
    for al in als:
        if len(al) < k:
            continue
        for i in range(len(contig) / 100):
            seg = contig.segment(i * 100, i * 100 + 100)
            if al.seg_to.inter(seg):
                if al.seg_to.left >= seg.left and al.seg_from.left > params.bad_end_length:
                    sys.stdout.write("B")
                elif al.seg_to.right <= seg.right and al.rc.seg_from.left > params.bad_end_length:
                    sys.stdout.write("E")
                else:
                    sys.stdout.write(al.reduce(target=seg).oneCharPI())
            else:
                sys.stdout.write("*")
        print ""


if __name__ == "__main__":
    contigs_file = sys.argv[1]
    contig_name = sys.argv[2]
    reads_file = sys.argv[3]
    dir = sys.argv[4]
    k = int(sys.argv[5])
    main(contigs_file, contig_name, reads_file, dir, k)

