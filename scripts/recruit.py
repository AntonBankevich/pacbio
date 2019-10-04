import os
import shutil
import subprocess
import sys

sys.path.append("py")


from common.SimpleGraph import SimpleGraph
from common.sequences import ContigStorage, Contig
from alignment.align_tools import DirDistributor, Aligner
from common import basic, params, SeqIO



def main(flye_dir, rf, dir, edge_id, k):
    params.technology = "nano"
    params.k = k
    basic.ensure_dir_existance(dir)
    basic.CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    all = 0
    good = 0
    print "Reading graph"
    graph = SimpleGraph().ReadGFA(os.path.join(flye_dir, "assembly_graph.gfa"))
    for e in graph.e[edge_id].fin.out:
        print e.id
    for e in graph.e[edge_id].start.inc:
        print e.id
    return
    print "Reading reads"
    reads = ContigStorage()
    for read in SeqIO.parse_fasta(open(rf, "r")):
        all += 1
        if len(read) > k * 1.2:
            good += 1
            reads.add(Contig(read.seq, read.id))
        if all % 1000000 == 0:
            print all, good
    print "Reading contigs"
    contigs = ContigStorage().loadFromFasta(open(cf, "r"), False)
    tl = sum(map(len, contigs.unique()))
    read_ids = set()
    aligner = Aligner(dd)
    print "Aligning reads"
    for al in aligner.localAlign(reads, contigs):
        if len(al) > k:
            read_ids.add(al.seg_from.contig.id)
    relevant_reads = ContigStorage()
    print "Choosing relevant reads"
    for rid in read_ids:
        relevant_reads.add(reads[rid])
    rrf = os.path.join(dir, "reads.fasta")
    handler = open(rrf, "w")
    relevant_reads.writeToFasta(handler)
    handler.close()
    assembly_dir = dd.nextDir()
    print "Assembling relevant reads"
    subprocess.check_call(
        ["./bin/flye", "-o", assembly_dir, "-t", 16, "--pacbio-raw", rrf, "--genome-size", tl + 5000])
    shutil.copy(os.path.join(assembly_dir, "assembly.fasta"), os.path.join(dir, "new_seqs.fasta"))
    shutil.copy(os.path.join(assembly_dir, "assembly_graph.gv"), os.path.join(dir, "graph.gv"))
    subprocess.check_call([". dotconv", dir])

if __name__ == "__main__":
    flye_dir = sys.argv[1]
    rf = sys.argv[2]
    dir = sys.argv[3]
    edge_id = sys.argv[4]
    k = int(sys.argv[5])
    main(flye_dir, rf, dir, edge_id, k)

