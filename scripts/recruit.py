import os
import shutil
import subprocess
import sys

sys.path.append("py")

from common.sequences import ContigStorage
from alignment.align_tools import DirDistributor, Aligner
from common import basic, params



def main(cf, rf, dir, k):
    params.k = k
    basic.ensure_dir_existance(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    reads = ContigStorage().loadFromFasta(open(rf, "r"), False)
    contigs = ContigStorage().loadFromFasta(open(cf, "r"), False)
    tl = sum(map(len, contigs.unique()))
    read_ids = set()
    aligner = Aligner(dd)
    for al in aligner.localAlign(reads, contigs):
        if len(al) > k:
            read_ids.add(al.seg_from.contig.id)
    relevant_reads = ContigStorage()
    for read in reads:
        if read.id in relevant_reads:
            relevant_reads.add(read)
    rrf = os.path.join(dir, reads)
    relevant_reads.writeToFasta(rrf)
    assembly_dir = dd.nextDir()
    subprocess.check_call(
        ["./bin/flye", "-o", assembly_dir, "-t", 16, "--pacbio-raw", rrf, "--genome-size", tl + 5000])
    shutil.copy(os.path.join(assembly_dir, "assembly.fasta"), os.path.join(dir, "new_seqs.fasta"))
    shutil.copy(os.path.join(assembly_dir, "assembly_graph.gv"), os.path.join(dir, "graph.gv"))
    subprocess.check_call([". dotconv", dir])

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))

