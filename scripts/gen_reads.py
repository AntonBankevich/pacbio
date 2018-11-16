import sys
import os

import common.seq_records

sys.path.append("py")
import common.SeqIO as SeqIO

read_len = int(sys.argv[2])
cov = float(sys.argv[3])
for seq in SeqIO.parse_fasta(open(sys.argv[1], "r")):
    sys.stderr.write(seq.id + " " + str(len(seq)) + " " + str(int(len(seq) * cov / read_len)) + "\n")
    # if len(seq) > 100000000 or len(seq) < 10000000:
    #      continue
    cur = 100000
    for i in range(0, len(seq), int(read_len / cov)):
        if i > cur:
            sys.stderr.write(str(cur) + "\n")
            cur = cur * 3 / 2
        SeqIO.write(common.seq_records.SeqRecord(seq.seq[i:min(len(seq), i + read_len)], seq.id + "_" + str(i)), sys.stdout, "fasta")
    break

