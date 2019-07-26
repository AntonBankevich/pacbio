import os
import shutil
import sys

from common import SeqIO

dir = sys.argv[1]
name = sys.argv[2]
f_name = os.path.join(dir, name + ".fasta")
tmp_name = os.path.join(dir, "tmp.fasta")
tmp = open(tmp_name, "w")
for rec in SeqIO.parse_fasta(open(f_name, "r")):
    rec.id = name + rec.id[rec.id.find("/"):]
    SeqIO.write(rec, tmp, "fasta")
shutil.move(tmp_name, f_name)