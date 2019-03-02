import os
import sys

sys.path.append("py")


from alignment.align_tools import Aligner, DirDistributor
from common import basic, SeqIO
from common.sam_parser import Samfile
from common.sequences import ContigCollection, ReadCollection, AlignedRead


def main(contig_file, reads_file, sam_file, dir, contig_id):
    # type: (str, str, str, str, str) -> None
    basic.ensure_dir_existance(dir)
    contigs = ContigCollection()
    contigs.loadFromFasta(open(contig_file, "r"))
    print "Contigs loaded"
    contig = contigs[contig_id]
    read_names = set()
    for rec in Samfile(open(sam_file, "r")):
        read_names.add(rec.query_name)
    reads = ReadCollection()
    cnt = 0
    for rec in SeqIO.parse_fasta(open(reads_file, "r")):
        if rec.id in read_names:
            rec.id = "Read" + str(cnt)
            reads.add(AlignedRead(rec))
            cnt += 1
    reads.print_fasta(open(os.path.join(dir, "reads.fasta"), "w"))
    print "Reads loaded", len(reads)
    reads.addAllRC()
    print "RC added", len(reads)

    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    aligner.alignReadCollection(reads, contigs)
    print "Reads aligned", len(reads)
    reads = reads.inter(contig.asSegment())
    print "Reads filtered", len(reads)
    sorted_reads = sorted(list(reads.reads.values()), key = lambda read: read.alignmentsTo(contig.asSegment()).next().seg_to.left)
    for read in sorted_reads:
        print read
        for al in read.alignmentsTo(contig.asSegment()):
            print "\n".join(al.asMatchingStrings())

# nohup python scripts/repeat_analysis.py results/celegans/contigs.fasta results/celegans/celegans.sam results/celegans/analysis contig_49
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])