import sys

from common import SeqIO

rf = sys.argv[1]
outf = sys.argv[2]
total = int(sys.argv[3])
reads = list(SeqIO.parse_fasta(open(rf, "r")))
reads = sorted(reads, key = lambda read: -len(read))
out = open(outf, "w")
for read in reads:
    if total <= 0:
        break
    total -= len(read)
    SeqIO.write(read, out, "fasta")
out.close()
