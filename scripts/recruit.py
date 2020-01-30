import os
import shutil
import subprocess
import sys


sys.path.append("py")

from common.parsing import parseUPaths
from common.seq_records import NamedSequence
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
    print "Reading graph"
    graph = SimpleGraph().ReadGFA(os.path.join(flye_dir, "assembly_graph.gfa"))
    print "Parsing edge mapping"
    id_map = parseUPaths(flye_dir)
    edge_ids = edge_id.split(",")
    print "Extracting relevant graph component"
    res = open(os.path.join(dir, "contigs.fasta"), "w")
    unique = dict()
    for eid in edge_ids:
        for e in graph.v[graph.e[eid].start].inc:
            if basic.Normalize(e.id) in edge_ids:
                continue
            if len(e.seq) < 10000:
                if e.id.startswith("-"):
                    unique[e.id[1:]] = NamedSequence(basic.RC(e.seq), e.id[1:])
                else:
                    unique[e.id] = NamedSequence(e.seq, e.id)
            else:
                if e.id.startswith("-"):
                    unique[e.id[1:] + "l"] = NamedSequence(basic.RC(e.seq[:5000]), e.id[1:] + "l")
                else:
                    unique[e.id + "r"] = NamedSequence(e.seq[-5000:], e.id + "r")
        for e in graph.v[graph.e[eid].end].out:
            if basic.Normalize(e.id) in edge_ids:
                continue
            if len(e.seq) < 10000:
                if e.id.startswith("-"):
                    unique[e.id[1:]] = NamedSequence(basic.RC(e.seq), e.id[1:])
                else:
                    unique[e.id] = NamedSequence(e.seq, e.id)
            else:
                if e.id.startswith("-"):
                    unique[e.id[1:] + "r"] = NamedSequence(basic.RC(e.seq[-5000:]), e.id[1:] + "r")
                else:
                    unique[e.id + "l"] = NamedSequence(e.seq[:5000], e.id + "l")


    for c in unique.values():
        print c.id
        SeqIO.write(c, res, "fasta")
    res.close()
    old_ids = []
    for eid in edge_ids:
        for olde in id_map[eid[len("edge_"):]]:
            old_ids.append(basic.Normalize(olde))
    print "Finding reads that align to", edge_ids
    print "Old ids:", old_ids
    relevant_read_ids = set()
    for s in open(os.path.join(flye_dir, "20-repeat", "read_alignment_dump"), "r").readlines():
        s = s.split()
        if s[0] != "Aln":
            continue
        if s[6].split("_")[1] in old_ids:
            relevant_read_ids.add(s[2][1:])
            print s[2][1:], s[6].split("_")[1]
    print "Reading reads"
    res = open(os.path.join(dir, "reads.fasta"), "w")
    for read in SeqIO.parse_fasta(open(rf, "r")):
        if read.id in relevant_read_ids and len(read) > k * 1.2:
            SeqIO.write(read, res, "fasta")
    res.close()



if __name__ == "__main__":
    flye_dir = sys.argv[1]
    rf = sys.argv[2]
    dir = sys.argv[3]
    edge_id = sys.argv[4]
    k = int(sys.argv[5])
    main(flye_dir, rf, dir, edge_id, k)

