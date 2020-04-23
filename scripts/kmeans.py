import os
import shutil
import subprocess
import sys

sys.path.append("py")
from alignment.polishing import Polisher

from common.seq_records import NamedSequence
from common.SimpleGraph import SimpleGraph, Edge
from common.sequences import ContigStorage, Contig
from alignment.align_tools import DirDistributor, Aligner
from common import basic, params, SeqIO
from common.line_align import Scorer

def fixAlDir(als, contig):
    res = []
    for al in als:
        if al.seg_to.contig.id != contig.id:
            res.append(al.rc)
        else:
            res.append(al)
    return res

def prolong(aligner, polisher, contig, reads):
    als = list(aligner.overlapAlign(reads.unique(), ContigStorage([contig])))
    contig, als = polisher.polishEnd(fixAlDir(als, contig), min_cov=5)
    contig, als = polisher.polishEnd([al.rc for al in als], min_cov=5)
    return contig


def main(contigs_file, contig_name, reads_file, dir, k, initial_reads1, initial_reads2):
    basic.ensure_dir_existance(dir)
    basic.CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    contigs = ContigStorage().loadFromFasta(open(contigs_file, "r"), False)
#    contig = contigs[contig_name].asSegment().prefix(length=2000).asContig()
    contig = contigs[contig_name]
    reads = ContigStorage().loadFromFasta(open(reads_file, "r"), False)
    reads1 = ContigStorage()
    reads2 = ContigStorage()
    cnt = 0
    for read in reads.unique():
        cnt += 1
#        if cnt % 2 == 0:
        if read.id in initial_reads1:
            reads1.add(read)
        elif read.id in initial_reads2:
            reads2.add(read)
    polisher = Polisher(aligner, dd)
    contig1 = contig
    contig2 = contig
    scorer = Scorer()
    for i in range(3):
        diff = 0
        print "Iteration", i
        als1 = fixAlDir(aligner.overlapAlign(reads1.unique(), ContigStorage([contig])), contig)
        als2 = fixAlDir(aligner.overlapAlign(reads2.unique(), ContigStorage([contig])), contig)
        contig1 = Contig(polisher.polishSmallSegment(contig.asSegment(), als1).seg_from.Seq(), "1")
        contig2 = Contig(polisher.polishSmallSegment(contig.asSegment(), als2).seg_from.Seq(), "2")
        al = aligner.overlapAlign([contig1], ContigStorage([contig2])).next()
        als1 = fixAlDir(aligner.overlapAlign(reads.unique(), ContigStorage([contig1])), contig1)
        als1 = filter(lambda al: len(al.seg_to) > len(al.seg_to.contig) - 100, als1)
        als2 = fixAlDir(aligner.overlapAlign(reads.unique(), ContigStorage([contig2])), contig2)
        als2 = filter(lambda al: len(al.seg_to) > len(al.seg_to.contig) - 100, als2)
        als1 = sorted(als1, key = lambda al: al.seg_from.contig.id)
        als2 = sorted(als2, key = lambda al: al.seg_from.contig.id)
        reads1 = ContigStorage()
        reads2 = ContigStorage()
        dp = scorer.accurateScore(al.matchingSequence(), 10) #1 - al.percentIdentity()
        als_map = dict()
        for al in als1:
            als_map[al.seg_from.contig.id] = [al]
        for al in als2:
            if al.seg_from.contig.id in als_map:
                als_map[al.seg_from.contig.id].append(al)
        com_res = []
        diffs = []
        for tmp_als in als_map.values():
            if len(tmp_als) != 2:
                continue
            al1 = tmp_als[0]
            al2 = tmp_als[1]
            print al1, al2
            assert al1.seg_from.contig == al2.seg_from.contig
            pi1 = scorer.accurateScore(al1.matchingSequence(), 10) # al1.percentIdentity()
            pi2 = scorer.accurateScore(al2.matchingSequence(), 10) # al2.percentIdentity()
            com_res.append((al1, al2, pi1 - pi2))
            diffs.append(pi1 - pi2)
        diffs = sorted(diffs)
        th1 = diffs[len(diffs) / 4]
        th2 = diffs[len(diffs) * 3 / 4]
        print "Thresholds:", th1, th2
        for al1, al2, diff in com_res:
            if diff < th1:
                reads1.add(al1.seg_from.contig)
            elif diff > th2:
                reads2.add(al2.seg_from.contig)
#           if pi1 > pi2 + dp / 4:
#               reads1.add(al1.seg_from.contig)
#           elif pi2 > pi1 + dp / 4:
#               reads2.add(al2.seg_from.contig)
#           diff += abs(pi1 - pi2)
        print float(diff) / len(als1), len(reads1) / 2, len(reads2) / 2
    al = aligner.overlapAlign([contig1], ContigStorage([contig2])).next()
    print al
    print "\n".join(al.asMatchingStrings2())
    for read in reads1:
        if read.id in initial_reads1:
            sys.stdout.write(read.id + " ")
    print ""
    for read in reads2:
        if read.id in initial_reads2:
            sys.stdout.write(read.id + " ")
    print ""
    contig1 = prolong(aligner, polisher, contig1, reads1)
    contig2 = prolong(aligner, polisher, contig2, reads2)
    contig1.id = "1"
    contig2.id = "2"
    out = open(os.path.join(dir, "copies.fasta"), "w")
    SeqIO.write(contig1, out, "fasta")
    SeqIO.write(contig2, out, "fasta")
    out.close()
    out = open(os.path.join(dir, "reads1.fasta"), "w")
    for read in reads1.unique():
        SeqIO.write(read, out, "fasta")
    out.close()
    out = open(os.path.join(dir, "reads2.fasta"), "w")
    for read in reads2.unique():
        SeqIO.write(read, out, "fasta")
    out.close()
    print "Finished"




if __name__ == "__main__":
    contigs_file = sys.argv[1]
    contig_name = sys.argv[2]
    reads_file = sys.argv[3]
    dir = sys.argv[4]
    k = int(sys.argv[5])
    initial_reads1 = sys.argv[6].split(",")
    initial_reads2 = sys.argv[7].split(",")
    main(contigs_file, contig_name, reads_file, dir, k, initial_reads1, initial_reads2)

