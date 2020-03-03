import sys

sys.path.append("py")
from common.alignment_storage import ReadCollection


def cutReads(rf, cut_len):
    # type: (str, int) -> None
    rc = ReadCollection().loadFromFasta(open(rf, "r"), None, cut_len)
    rc.print_fasta(sys.stdout)

if __name__ == "__main__":
    rf = sys.argv[1]
    cut_len = int(sys.argv[2])
    cutReads(rf, cut_len)
