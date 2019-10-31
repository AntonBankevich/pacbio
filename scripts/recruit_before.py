import os
import shutil
import subprocess
import sys


sys.path.append("py")

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
    graph = SimpleGraph().ReadDot(os.path.join(flye_dir, "20-repeat", "graph_before_rr.gv"))
    graph.FillSeq(os.path.join(flye_dir, "20-repeat", "graph_before_rr.fasta"), True)
    print "Extracting relevant graph component"
    res = open(os.path.join(dir, "contigs.fasta"), "w")
    unique = dict()
    edge_ids = edge_id.split(",")
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
        for e in graph.v[graph.e[eid].fin].out:
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
    old_ids = list(edge_ids)
    print "Finding reads that align to", edge_ids
    relevant_read_ids = set()
    for s in open(os.path.join(flye_dir, "20-repeat", "read_alignment_dump"), "r").readlines():
        s = s.split()
        if s[0] != "Aln":
            continue
        if s[6].split("_")[1] in old_ids:
            relevant_read_ids.add(s[2][1:])
            print s[2][1:], s[6].split("_")[1]
    print "Reading reads"
    res_reads = []
    res = open(os.path.join(dir, "reads.fasta"), "w")
    for read in SeqIO.parse_fasta(open(rf, "r")):
        if read.id in relevant_read_ids and len(read) > k * 1.2:
            res_reads.append(read)
            SeqIO.write(read, res, "fasta")
    res.close()
    random_down = open(os.path.join(dir, "random_down.fasta"), "w")
    cnt = 0
    for read in res_reads:
        if cnt % 5 == 0:
            SeqIO.write(read, random_down, "fasta")
    random_down.close()
    res_reads = sorted(res_reads, key = lambda read: -len(read))
    largest_down = open(os.path.join(dir, "largest.fasta"), "w")
    all = sum(map(len, res_reads))
    cnt = 0
    for read in res_reads:
        if cnt > all / 5:
            break
        SeqIO.write(read, largest_down, "fasta")
    largest_down.close()



if __name__ == "__main__":
    flye_dir = sys.argv[1]
    rf = sys.argv[2]
    dir = sys.argv[3]
    edge_id = sys.argv[4]
    k = int(sys.argv[5])
    main(flye_dir, rf, dir, edge_id, k)
