import sys
sys.path.append("py")

from common import SeqIO
from common.sequences import ContigStorage


def printSegs(f, segs):
    c = ContigStorage().loadFromFasta(open(f, "r"), False)
    for seg in segs:
        if seg[2] == 0:
            seg[2] = len(c[seg[0]])
        SeqIO.write(c[seg[0]].segment(seg[1], seg[2]).asContig(), sys.stdout, "fasta")
        
if __name__ == "__main__":
    segs = []
    for i in range(len(sys.argv) / 3):
        segs.append([sys.argv[i * 3 + 2], int(sys.argv[i * 3 + 3]), int(sys.argv[i * 3 + 4])])
    printSegs(sys.argv[1], segs)
