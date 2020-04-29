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
    als = list(aligner.localAlign(reads.unique(), contigs))
    tmp = []
    for al in als:
        if al.seg_to.contig != contig:
            al = al.rc
        tmp.append(al)
    als = tmp
    als = sorted(als, key=lambda al: al.seg_to.left / 50 * 1000000 + al.seg_to.right - al.seg_to.left)
    counts = dict()
    for al in als:
        counts[al.seg_from.contig.id] = 0
    for al in als:
        if len(al) > k:
            counts[al.seg_from.contig.id] += 1
    w = 20
    f = open(os.path.join(dir, "reads.fasta"), "w")
    over = set()
    inter = set()
    for al in als:
        if len(al) < k:
            continue
        inter.add(basic.Normalize(al.seg_from.contig.id))
        if not al.contradictingRTC():
            over.add(basic.Normalize(al.seg_from.contig.id))
        m = al.matchingSequence(True)
        tmp = []
        for i in range(len(contig) / w + 1):
            tmp.append([])
        for a, b in m.matches:
            tmp[b / w].append((a, b))
        for i in range(len(contig) / w):
            if i + 1 < len(tmp) and len(tmp[i + 1]) > 0:
                tmp[i].append(tmp[i+1][0])
        for i in range(len(contig) / w):
            seg = contig.segment(i * w, i * w + w)
            if al.seg_to.inter(seg):
                if al.seg_to.left >= seg.left and al.seg_from.left > params.bad_end_length:
                    sys.stdout.write("B")
                elif al.seg_to.right <= seg.right and al.rc.seg_from.left > params.bad_end_length:
                    sys.stdout.write("E")
                else:
                    if len(tmp[i]) == 0:
                        sys.stdout.write("*")
                    else:
                        a = tmp[i][-1][0] - tmp[i][0][0]
                        b = tmp[i][-1][1] - tmp[i][0][1]
                        if a - b > 30:
                            sys.stdout.write("I")
                        elif a - b > 15:
                            sys.stdout.write("i")
                        elif a - b < -30:
                            sys.stdout.write("D")
                        elif a - b < -15:
                            sys.stdout.write("d")
                        else:
                            sys.stdout.write(str(min(8, max(a, b) + 1 - len(tmp[i]))))
            else:
                sys.stdout.write("*")
        print " ", al.seg_from.contig.id, counts[al.seg_from.contig.id], al.contradictingRTC()
    print inter
    for rid in inter:
        SeqIO.write(reads[rid], f, "fasta")
        print rid, reads[rid]
    f.close()
    f = open(os.path.join(dir, "reads_over.fasta"), "w")
    for rid in over:
        SeqIO.write(reads[rid], f, "fasta")
    f.close()


if __name__ == "__main__":
    contigs_file = sys.argv[1]
    contig_name = sys.argv[2]
    reads_file = sys.argv[3]
    dir = sys.argv[4]
    k = int(sys.argv[5])
    main(contigs_file, contig_name, reads_file, dir, k)

