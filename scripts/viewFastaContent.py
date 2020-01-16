import sys
sys.path.append("py")
from common import SeqIO

if __name__ == "__main__":
    contigs = list(SeqIO.parse_fasta(open(sys.argv[1], "r")))
    print "Total:", sum(map(len, contigs)), "nucleotides in", len(contigs), "contigs."
    for contig in contigs:
        print contig.id, len(contig)
