import os
import shutil
import subprocess
import sys

from typing import List, Dict

sys.path.append("py")
from alignment.polishing import Polisher


from common.seq_records import NamedSequence
from common.SimpleGraph import SimpleGraph, Edge
from common.sequences import ContigStorage, Contig
from alignment.align_tools import DirDistributor, Aligner
from common import basic, params, SeqIO


def uniquePathForward(graph, edge, k):
    # type: (SimpleGraph, Edge, int) -> str
    print "Forward", edge.id, k
    if edge.len >= k:
        return edge.seq[:k]
    v = graph.v[edge.end]
    if len(v.out) == 1:
        return edge.seq + uniquePathForward(graph, v.out[0], k - edge.len)
    else:
        return edge.seq


def uniquePathBackward(graph, edge, k):
    # type: (SimpleGraph, Edge, int) -> str
    print "Backward", edge.id, k
    if edge.len >= k:
        return edge.seq[-k:]
    v = graph.v[edge.start]
    if len(v.inc) == 1:
        return uniquePathBackward(graph, v.inc[0], k - edge.len) + edge.seq
    else:
        return edge.seq

def parseReadDump(fn):
    for s in open(fn, "r").readlines():
        s = s.split()
        if s[0] != "Aln":
            continue
        yield s[2][1:], s[6].split("_")[1]

def fixAlDir(als, contig):
    res = []
    for al in als:
        if al.seg_to.contig.id != contig.id:
            res.append(al.rc)
        else:
            res.append(al)
    return res


def splitSeg(aligner, seg, mult, all_reads_list):
    all_reads = ContigStorage()
    base = seg.asContig()
    for al in fixAlDir(aligner.overlapAlign(all_reads_list, ContigStorage([base])), base):
        all_reads.add(al.seg_from.contig)
    all_reads_list = list(all_reads.unique())
    split_reads = []
    split_contigs = []
    for i in range(mult):
        split_reads.append([])
        split_contigs.append(base)
    cnt = 0
    for read in all_reads:
        split_reads[cnt % mult].append(read)
    polisher = Polisher(aligner, aligner.dir_distributor)
    diff = 0
    for i in range(10):
        print "Iteration", i
        split_contigs = []
        for reads in split_reads:
            tmp_als = fixAlDir(aligner.overlapAlign(reads, ContigStorage([base])), base)
            split_contigs.append(Contig(polisher.polishSmallSegment(base.asSegment(), tmp_als).seg_from.Seq(), str(len(split_contigs))))
        als = []
        for contig in split_contigs:
            als.append(fixAlDir(aligner.overlapAlign(all_reads_list, ContigStorage([contig])), contig))
            als[-1] = sorted(als[-1], key = lambda al: al.seg_from.contig.id)
        for i in range(mult):
            split_reads[i] = []
        for i in range(len(all_reads_list)):
            best_pi = 0
            best = 0
            for j in range(mult):
                assert als[j][i].seg_from.contig == als[0][i].seg_from.contig
                pi = als[j][i].percentIdentity()
                if best_pi > pi:
                    best = j
                    best_pi = pi
                split_reads[best].append(all_reads_list[i])
        print diff, " ".join(map(str, map(len, split_reads)))
    maxpi = 0
    print "pi matrix:"
    for i in range(mult):
        for j in range(mult):
            al = aligner.overlapAlign([split_contigs[i]], ContigStorage([split_contigs[j]])).next()
            sys.stdout.write(str(al.percentIdentity()) + " ")
            maxpi = max(maxpi, al.percentIdentity())
        print ""
    print "Maxpi:", maxpi
    if maxpi < 0.985:
        return zip(split_contigs, split_reads)
    else:
        return None


def splitRepeat(aligner, seq, mult, all_reads_list, min_contig_length):
    base = Contig(seq, "base")
    for i in range(len(seq) / min_contig_length):
        res = splitSeg(aligner, base.segment(i * min_contig_length, i * min_contig_length + min_contig_length), mult, all_reads_list)
        if res is not None:
             return res
    res = splitSeg(aligner, base.asSegment().suffix(length=min(min_contig_length, len(seq))), mult, all_reads_list)
    return res


def main(flye_dir, rf, dir, edge_id, to_resolve, min_contig_length):
    params.technology = "nano"
    basic.ensure_dir_existance(dir)
    basic.CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    print " ".join(sys.argv)
    print "Reading graph"
    graph = SimpleGraph().ReadDot(os.path.join(flye_dir, "20-repeat", "graph_before_rr.gv"))
    graph.FillSeq(os.path.join(flye_dir, "20-repeat", "graph_before_rr.fasta"), True)
    print "Extracting relevant graph component"
    edge_ids = edge_id.split(",")
    to_resolve = to_resolve.split(",")
    to_resolve = [(a, int(b)) for a, b in zip(to_resolve[0::2], to_resolve[1::2])]
    unique = uniqueNeighbours(edge_ids, graph, min_contig_length)

    if rf == "none":
        return
    print "Finding reads that align to", edge_ids
    reads_to_resolve = dict()  # type: Dict[str, List[str]]
    for eid, mult in to_resolve:
        reads_to_resolve[eid] = []
    for unique_edge, initial in unique:
        reads_to_resolve[initial] = []
    relevant_read_ids = set()
    for rid, eid in parseReadDump(os.path.join(flye_dir, "20-repeat", "read_alignment_dump")):
        if eid in edge_ids:
            relevant_read_ids.add(rid)
            print rid, eid
    for rid, eid in parseReadDump(os.path.join(flye_dir, "20-repeat", "read_alignment_dump")):
        if rid in relevant_read_ids and eid in reads_to_resolve:
            reads_to_resolve[eid].append(rid)
    for eid in reads_to_resolve:
        reads_to_resolve[eid] = list(set(reads_to_resolve))
    print "Reading reads"
    res_reads = ContigStorage()
    res = open(os.path.join(dir, "reads.fasta"), "w")
    for read in SeqIO.parse_by_name(rf):
        if read.id in relevant_read_ids:
            res_reads.add(read)
            SeqIO.write(read, res, "fasta")
    res.close()
    random_down = open(os.path.join(dir, "random_down.fasta"), "w")
    cnt = 0
    for read in res_reads:
        if cnt % 5 == 0:
            SeqIO.write(read, random_down, "fasta")
        cnt += 1
    random_down.close()
    res = open(os.path.join(dir, "contigs.fasta"), "w")
    lcf = open(os.path.join(dir, "contigs.lc"), "w")
    for eid, mult in to_resolve:
        repeat_reads = [res_reads[rid] for rid in reads_to_resolve[eid]]
        split_contigs = splitRepeat(aligner, graph.e[eid].seq, mult, repeat_reads, min_contig_length)
        print "Edge", eid, "was split into", mult, "copies"
        for contig, contig_reads in split_contigs:
            print contig.id
            SeqIO.write(contig, res, "fasta")
            lcf.write(contig.id + "\n")
            lcf.write(" ".join([r.id for r in contig_reads]) + "\n")
    res = open(os.path.join(dir, "contigs.fasta"), "w")
    for unique_edge, initial in unique:
        print unique_edge.id
        SeqIO.write(unique_edge, res, "fasta")
        lcf.write(unique_edge.id + "\n")
        lcf.write(" ".join(reads_to_resolve[initial]) + "\n")
    res.close()


def uniqueNeighbours(edge_ids, graph, min_contig_length):
    unique = []
    for eid in edge_ids:
        print "Finding neighbours of", eid
        for e in graph.v[graph.e[eid].start].inc:
            if basic.Normalize(e.id) in edge_ids:
                continue
            id = basic.Normalize(e.id)
            if len(e.seq) < min_contig_length + params.bad_end_length:
                seq = uniquePathForward(graph, e, min_contig_length + params.bad_end_length)
                id = id + "p"
                # seq = e.seq
            else:
                seq = e.seq[-min_contig_length - params.bad_end_length:]
                if e.id.startswith("-"):
                    id = id + "l"
                else:
                    id = id + "r"
            if e.id.startswith("-"):
                seq = basic.RC(seq)
            print "Right neighbour", eid, id
            unique.append((NamedSequence(seq, id), basic.Normalize(e.id)))
        for e in graph.v[graph.e[eid].end].out:
            if basic.Normalize(e.id) in edge_ids:
                continue
            id = basic.Normalize(e.id)
            if len(e.seq) < min_contig_length + params.bad_end_length:
                seq = uniquePathBackward(graph, e, min_contig_length + params.bad_end_length)
                id = id + "p"
                # seq = e.seq
            else:
                seq = e.seq[:min_contig_length + params.bad_end_length]
                if e.id.startswith("-"):
                    id = id + "r"
                else:
                    id = id + "l"
            if e.id.startswith("-"):
                seq = basic.RC(seq)
            print "Left neighbour", eid, id
            unique.append((NamedSequence(seq, id), basic.Normalize(e.id)))
    return unique


if __name__ == "__main__":
    flye_dir = sys.argv[1]
    rf = sys.argv[2]
    dir = sys.argv[3]
    edge_id = sys.argv[4]
    unique_edges = sys.argv[5]
    k = int(sys.argv[6])
    main(flye_dir, rf, dir, edge_id, unique_edges, k)

