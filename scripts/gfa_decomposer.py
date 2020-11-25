import os
import sys
sys.path.append("py")

from common import basic, SeqIO
from common.SimpleGraph import SimpleGraph

if __name__ == "__main__":
    g = SimpleGraph()
    g.ReadGFA(sys.argv[1])
    fasta = open(sys.argv[2] + ".fasta", "w")
    dot = open(sys.argv[2] + ".dot", "w")
    if "merge" in sys.argv:
        g = g.Merge()
    g.Print(dot)
    g.PrintFasta(fasta)
    fasta.close()
    dot.close()
