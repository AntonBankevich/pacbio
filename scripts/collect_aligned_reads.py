import sys
import os
sys.path.append("py")

from common import sam_parser
from common.basic import RC
from flye_tools.alignment import make_alignment

import common.SeqIO as SeqIO

# # 1 - reads; 2 - dump
# dump = open(sys.argv[2]).readlines()
# d = dict()
# for s in dump:
#     s = s.strip()
#     s = s.split()
#     d[s[0]] = s[4]
# for rec in SeqIO.parse_fasta(open(sys.argv[1])):
#     id = rec.id.split()[0]
#     if id in d:
#         if d[id] == "-":
#             rec.seq = RC(rec.seq)
#         SeqIO.write(rec, sys.stdout, "fasta")
#         sys.stderr.write(id + "_" + str(d[id]) + "\n")

# 1 - contigs names 2 - contigs file 3- reads file 4 - dir
# 1 - contigs names 2 - contigs file 3- reads file 4 - dir

dir = os.path.join(sys.argv[4])
if not os.path.exists(dir):
    os.makedirs(dir)

names = sys.argv[1].split(";")
tmp = dict()
for s in names:
    if s.endswith("RC"):
        tmp[s[:-2]] = True
    else:
        tmp[s] = False
names = tmp

contigs_file = os.path.join(sys.argv[4], "contigs.fasta")
contigs_handler = open(contigs_file, "w")
for rec in SeqIO.parse_fasta(open(sys.argv[2], "r")):
    if rec.id in names:
        if names[rec.id]:
            rec.seq = RC(rec.seq)
            rec.id += "RC"
        SeqIO.write(rec, contigs_handler, "fasta")
contigs_handler.close()

alignment_dir = os.path.join(sys.argv[4], "alignment")
if not os.path.exists(alignment_dir):
    os.makedirs(alignment_dir)

alignment = os.path.join(sys.argv[4], "alignment.sam")
make_alignment(contigs_file, [sys.argv[3]], 8,
               alignment_dir, "pacbio", alignment)

aligned = set()
for rec in sam_parser.Samfile(open(alignment, "r")):
    if rec.is_unmapped:
        continue
    aligned.add(rec.query_name)

output = open(os.path.join(dir, "filtered_reads.fasta"), "w")
for rec in SeqIO.parse_fasta(open(sys.argv[3], "r")):
    if rec.id.split()[0] in aligned:
        SeqIO.write(rec, output, "fasta")
output.close()


