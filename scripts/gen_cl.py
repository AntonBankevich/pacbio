import sys

import os

sys.path.append("py")

from common import basic

indir = sys.argv[1]
outdir = sys.argv[2]
reads_dir = sys.argv[3]
basic.ensure_dir_existance(outdir)

for dir in sorted(os.listdir(indir)):
    reads = os.path.join(indir, dir, dir + ".fasta")
    contigs = os.path.join(indir, dir, "2-repeat", "graph_final.fasta")
    dot = os.path.join(indir, dir, "2-repeat", "graph_final.dot")
    out = os.path.join(outdir, dir)
    print "python py/dag_resolve/main.py " + " ".join([os.path.join(indir, dir), "\"\"", out])