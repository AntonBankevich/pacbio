import os
import sys

sys.path.append("py")
from common import basic, SeqIO
from common.sequences import ContigCollection

def main(ref_file, contig_name, rlen, cov, dir):
    basic.ensure_dir_existance(dir)

    ref = ContigCollection().loadFromFasta(open(ref_file, "r"), False)[contig_name]
    contig_file = open(os.path.join(dir, contig_name + ".fasta"), "w")
    SeqIO.write(ref, contig_file, "fasta")
    contig_file.close()

    reads_file = open(os.path.join(dir, contig_name + "_reads.fasta"))
    for i in range(0, len(ref), max(1, rlen / cov)):
        read = ref.segment(i, min(i + rlen, len(ref))).asNamedSequence()
        SeqIO.write(reads_file, read, "fasta")
    reads_file.close()



if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])