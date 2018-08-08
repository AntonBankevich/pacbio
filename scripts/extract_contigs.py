import sys
import os
sys.path.append("py")
import common.SeqIO as SeqIO
from common.basic import RC

names = sys.argv[1].split(";")
tmp = dict()
for s in names:
    s = s.split(",")
    tmp[s[0]] = s
names = tmp

for rec in SeqIO.parse_fasta(open(sys.argv[2], "r")):
    if rec.id in names:
        s = names[rec.id]
        if "RC" in s:
            rec.seq = RC(rec.seq)
        if "end" in s:
            rec.seq = rec.seq[-1000:]
        if "start" in s:
            rec.seq = rec.seq[:1000]
        SeqIO.write(rec, sys.stdout, "fasta")