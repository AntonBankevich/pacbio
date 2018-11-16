import sys
import os

import common.seq_records

sys.path.append("py")
import common.SeqIO as SeqIO
from common.basic import RC
# 1 - reads; 2 - dump
interest = ["5"]
repeat = None
cur = None
dump = open(sys.argv[2]).readlines()
d = dict()
for s in dump:
    s = s.strip()
    if s == "":
        continue
    if s.startswith("#"):
        s = s[1:].split()
        if s[0] =="Repeat":
            repeat = s[1]
        if s[0] in ["All", "Input", "Output"]:
            cur = s[1]
    else:
        if repeat not in interest:
            continue
        sign = s[0]
        s = s[1:]
        if s not in d:
            d[s] = []
        d[s].append((repeat, cur, sign))
for rec in SeqIO.parse_fasta(open(sys.argv[1])):
    id = rec.id.split()[0]
    if id in d:
        tmp = d[id]
        if ("reads", "-") in [(a[1], a[2]) for a in tmp]:
            rec.seq = RC(rec.seq)
        SeqIO.write(common.seq_records.SeqRecord(rec.seq, id + "_" + str(d[id])), sys.stdout, "fasta")
        sys.stderr.write(id + "_" + str(d[id]) + "\n")