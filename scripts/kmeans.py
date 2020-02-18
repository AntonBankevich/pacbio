import os
import shutil
import subprocess
import sys

from alignment.polishing import Polisher

sys.path.append("py")

from common.seq_records import NamedSequence
from common.SimpleGraph import SimpleGraph, Edge
from common.sequences import ContigStorage, Contig
from alignment.align_tools import DirDistributor, Aligner
from common import basic, params, SeqIO

def fixAlDir(als, contig):
    res = []
    for al in als:
        if al.seg_to.contig.id != contig.id:
            res.append(al.rc)
        else:
            res.append(al)
    return res

def main(contigs_file, contig_name, reads_file, dir, k, initial_reads):
    basic.ensure_dir_existance(dir)
    basic.CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    contigs = ContigStorage().loadFromFasta(open(contigs_file, "r"), False)
    contig = contigs[contig_name]
    reads = ContigStorage().loadFromFasta(open(reads_file, "r"), False)
    reads1 = ContigStorage()
    reads2 = ContigStorage()
    for read in reads:
        if read.id in initial_reads:
            reads1.add(read)
        else:
            reads2.add(read)
    polisher = Polisher(aligner, dd)
    contig1 = contig
    contig2 = contig
    diff = 0
    for i in range(10):
        print "Iteration", i
        als1 = fixAlDir(aligner.overlapAlign(reads1.unique(), ContigStorage([contig])), contig)
        als2 = fixAlDir(aligner.overlapAlign(reads2.unique(), ContigStorage([contig])), contig)
        contig1 = Contig(polisher.polishSmallSegment(contig.asSegment(), als1).seg_from.Seq(), "1")
        contig2 = Contig(polisher.polishSmallSegment(contig.asSegment(), als2).seg_from.Seq(), "2")
        als1 = fixAlDir(aligner.overlapAlign(reads.unique(), ContigStorage([contig1])), contig1)
        als2 = fixAlDir(aligner.overlapAlign(reads.unique(), ContigStorage([contig2])), contig2)
        als1 = sorted(als1, key = lambda al: al.seg_from.contig.id)
        als2 = sorted(als2, key = lambda al: al.seg_from.contig.id)
        reads1 = ContigStorage()
        reads2 = ContigStorage()
        for al1, al2 in zip(als1, als2):
            assert al1.seg_from.contig == al2.seg_from.contig
            pi1 = al1.percentIdentity()
            pi2 = al2.percentIdentity()
            if pi1 > pi2:
                reads1.add(al1.seg_from.contig)
            else:
                reads2.add(al2.seg_from.contig)
            diff += abs(pi1 - pi2)
        print diff, len(reads1), len(reads2)
    al = aligner.overlapAlign([contig1], ContigStorage([contig2])).next()
    print al
    print "\n".join(al.asMatchingStrings2())



if __name__ == "__main__":
    contigs_file = sys.argv[1]
    contig_name = sys.argv[2]
    reads_file = sys.argv[3]
    dir = sys.argv[4]
    k = int(sys.argv[5])
    initial_reads = sys.argv[6].split(",")
    main(contigs_file, contig_name, reads_file, dir, k, initial_reads)

