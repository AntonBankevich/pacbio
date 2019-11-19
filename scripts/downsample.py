import sys
sys.path.append("py")

from common import SeqIO

rf = sys.argv[1]
outf = sys.argv[2]
total = int(sys.argv[3])
print "Reading reads"
reads = list(SeqIO.parse_fasta(open(rf, "r")))
print "Sorting reads"
reads = sorted(reads, key = lambda read: -len(read))
print "Printing reads"
out = open(outf, "w")
for read in reads:
    if total <= 0:
        break
    total -= len(read)
    SeqIO.write(read, out, "fasta")
out.close()
