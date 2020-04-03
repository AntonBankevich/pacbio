import sys
sys.path.append("py")
from common.sequences import ContigStorage

if __name__ == "__main__":
    fname = sys.argv[1]
    contigs = ContigStorage().loadFromFasta(open(fname, "r"), False)
    contigs.writeToFasta(sys.stdout)