import json
import sys
import os

import pickle

sys.path.append("py")
from common import sam_parser, basic
from common.basic import RC
from flye_tools.alignment import make_alignment
from flye_tools import polysh_job

import common.SeqIO as SeqIO


def find_range(pairs, pos, val, range):
    rpos = pos
    while rpos + 1 < len(pairs) and pairs[rpos + 1][1] < val + range:
        rpos += 1
    lpos = pos
    while lpos - 1 >= 0 and pairs[lpos - 1][1] > val - range:
        lpos -= 1
    return pairs[lpos][0], pairs[rpos][0]

def find_val(pairs, pos, val):
    while pos + 1 < len(pairs) and pairs[pos + 1][1] < val:
        pos += 1
    return pos

sam_file = sys.argv[1]
contig = list(SeqIO.parse_fasta(open(sys.argv[2], "r")))[0].seq.upper()
pos = sorted(map(int, sys.argv[3].split(",")))
dir = sys.argv[4]
basic.ensure_dir_existance(dir)

result = []
dist = 10
for i in range(len(pos)):
    result.append([contig[pos[i] - dist: pos[i] + dist + 1]])
for rec in sam_parser.Samfile(open(sam_file)):
    if rec.is_unmapped:
        continue
    pairs = list(sam_parser.ParseCigar(rec.cigar, len(rec.seq), rec.pos))
    cpos = 0
    for i in range(len(pos)):
        cpos = find_val(pairs, cpos, pos[i])
        if cpos == 0 or cpos == len(pos) - 1:
            continue
        l, r = find_range(pairs, cpos, pos[i], dist)
        result[i].append(rec.seq[l: r+1].upper() + " " + rec.query_name)

for i in range(len(pos)):
    f = open(os.path.join(dir, str(pos[i]) + ".txt"), "w")
    f.write("\n".join(result[i]))
    f.close()
