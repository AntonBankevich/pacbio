import os
import sys
sys.path.append("py")
from common import SeqIO

def readsN50(dir):
    for fn in os.listdir(dir):
        tmp = []
        f = os.path.join(dir, fn, fn + ".fasta")
        for rec in SeqIO.parse_fasta(open(f, "r")):
            if len(tmp) >= 1000:
                break
            tmp.append(len(rec))
        print fn, sorted(tmp)[len(tmp) / 2]


if __name__ == "__main__":
    dir = sys.argv[1]
    readsN50(dir)