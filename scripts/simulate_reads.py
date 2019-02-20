import os
import sys

sys.path.append("py")
from common import basic, SeqIO
from common.sequences import ContigCollection

def main(ref_file, contig_name, rlen, cov, dir):
    contig_names = contig_name.split(";")
    basic.ensure_dir_existance(dir)
    all_contigs = ContigCollection().loadFromFasta(open(ref_file, "r"), False)
    for contig_name in contig_names:
        if contig_name not in all_contigs:
            print "Fail", contig_name
            continue
        ref = all_contigs[contig_name]
        contig_file_name = os.path.join(dir, contig_name + ".fasta")
        contig_file = open(contig_file_name, "w")
        SeqIO.write(ref, contig_file, "fasta")
        contig_file.close()

        reads_file_name = os.path.join(dir, contig_name + "_reads.fasta")
        reads_file = open(reads_file_name, "w")
        for i in range(0, len(ref), max(1, rlen / cov)):
            read = ref.segment(i, min(i + rlen, len(ref))).asNamedSequence()
            SeqIO.write(reads_file, read, "fasta")
        reads_file.close()
        print "Done", contig_name
        print contig_file_name
        print reads_file_name

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])