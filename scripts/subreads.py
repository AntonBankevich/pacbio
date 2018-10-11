import sys
sys.path.append("py")
from common.basic import *
from common import SeqIO


def extract(fname):
    res = dict()
    for read in SeqIO.parse_fasta(open(fname, "r")):
        res[read.id] = read.seq
    return res
