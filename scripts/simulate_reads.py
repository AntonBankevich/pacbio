import os
import sys

sys.path.append("py")
from common import basic, SeqIO
from common.sequences import ContigCollection

def main(ref_file, contig_size, rlen, cov, dir):
    basic.ensure_dir_existance(dir)
    all_contigs = ContigCollection().loadFromFasta(open(ref_file, "r"), False)
    contig_file_name = os.path.join(dir, "contigs.fasta")
    contig_file = open(contig_file_name, "w")
    reads_file_name = os.path.join(dir, "reads.fasta")
    reads_file = open(reads_file_name, "w")
    for ref in all_contigs.unique():
        if len(ref) < contig_size:
            continue
        SeqIO.write(ref, contig_file, "fasta")
        for i in range(0, len(ref), max(1, rlen / cov)):
            read = ref.segment(i, min(i + rlen, len(ref))).asNamedSequence()
            SeqIO.write(read, reads_file, "fasta")
    reads_file.close()
    contig_file.close()
    print "Done"
    print contig_file_name
    print reads_file_name

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
