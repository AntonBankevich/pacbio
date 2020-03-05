import os
import shutil
import subprocess
import sys

from sklearn.cluster import KMeans
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
    tmp = []
    for al in fixAlDir(aligner.overlapAlign(all_reads_list, ContigStorage([base])), base):
        if len(al.seg_to) < len(base) - 100:
            continue
        all_reads.add(al.seg_from.contig)
        tmp.append(al.seg_from.contig)
    all_reads_list = tmp
    split_reads = []
    split_contigs = []
    for i in range(mult):
        split_reads.append([])
        split_contigs.append(base)
    cnt = 0
    for read in all_reads_list:
        split_reads[cnt % mult].append(read)
    polisher = Polisher(aligner, aligner.dir_distributor)
    for i in range(10):
        print "Iteration", i
        split_contigs = []
        for reads in split_reads:
            tmp_als = fixAlDir(aligner.overlapAlign(reads, ContigStorage([base])), base)
            split_contigs.append(Contig(polisher.polishSmallSegment(base.asSegment(), tmp_als).seg_from.Seq(), str(len(split_contigs))))
        bestals = dict()
        for read in all_reads_list:
            bestals[read.id] = None
        for contig in split_contigs:
            for al in fixAlDir(aligner.overlapAlign(all_reads_list, ContigStorage([contig])), contig):
                if len(al.seg_to) < len(base) - 100:
                    continue
                if al.seg_from.contig.id not in bestals:
                    print bestals.keys()
                    print al
                if bestals[al.seg_from.contig.id] is None or al.percentIdentity() > bestals[al.seg_from.contig.id].percentIdentity():
                    bestals[al.seg_from.contig.id] = al
#            als.append(fixAlDir(aligner.overlapAlign(all_reads_list, ContigStorage([contig])), contig))
#            als[-1] = sorted(als[-1], key = lambda al: al.seg_from.contig.id)
        for i in range(mult):
            split_reads[i] = []
        for rid in bestals:
            al = bestals[rid]
            if al is None:
                print "Warning: no alignment for read", rid
            else:
                split_reads[int(al.seg_to.contig.id)].append(al.seg_from.contig)
        print " ".join(map(str, map(len, split_reads)))
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

w = 50

def toVector(al):
    res = []
    contig = al.seg_to
    m = al.matchingSequence(True)
    tmp = []
    for i in range(len(contig) / w + 1):
        tmp.append([])
    for a, b in m.matches:
        tmp[b / w].append((a, b))
    for i in range(len(contig) / w):
        if i + 1 < len(tmp) and len(tmp[i + 1]) > 0:
            tmp[i].append(tmp[i + 1][0])
    for i in range(len(contig) / w):
        seg = contig.segment(i * w, i * w + w)
        if al.seg_to.left >= seg.left and al.seg_from.left > params.bad_end_length:
            sys.stdout.write("B")
        elif al.seg_to.right <= seg.right and al.rc.seg_from.left > params.bad_end_length:
            sys.stdout.write("E")
        else:
            res.append(w - len(tmp[i]))
    return res

class ReadRecord:
    def __init__(self, al):
        self.read = al.seg_from
        self.al = al
        self.v = []

    def extend(self, v):
        self.v.extend(v)
        return self

def readsToVectors(aligner, reads_list, base):
    als = []
    rtv = dict()
    for al in fixAlDir(aligner.overlapAlign(reads_list, ContigStorage([base])), base):
        if len(al.seg_to) < len(base) - 30:
            continue
        else:
            als.append(al)
            rtv[al.seg_from.contig.id] = ReadRecord(al).extend(toVector(al))
    reads_list = [al.seg_from.contig for al in als]
    bases = [base]
    for base_al in als:
        rtr_als = []
        base_candidate = base_al.seg_from.asContig()
        for al in fixAlDir(aligner.overlapAlign(reads_list, ContigStorage([base_candidate])), base_candidate):
            if len(al.seg_to) < len(base_candidate) - 30:
                continue
            else:
                rtr_als.append(al)
        if len(rtr_als) == len(als):
            bases.append(base_candidate)
            for al in rtr_als:
                rtv[al.seg_from.contig.id].extend(toVector(al))
            if len(bases) > 10:
                break
    return rtv


def splitSegKmeans(aligner, seg, mult, all_reads_list):
    polisher = Polisher(aligner, aligner.dir_distributor)
    all_reads = ContigStorage()
    base = seg.asContig()
    tmp = []
    rtv = readsToVectors(aligner, all_reads_list, base)
    kmeans = KMeans(n_clusters=mult, precompute_distances=True)
    recs = list(rtv.values())
    result = kmeans.fit_predict(X=[rec.v for rec in recs])
    print result
    clusters = dict()
    for i, c in enumerate(result):
        if c not in clusters:
            clusters[c] = []
        clusters[c].append(recs[i].al)
    for c in clusters.values():
        print str(c), ":", len(c)
    split_contigs = []
    split_reads = []
    for c in clusters.values():
        split_contigs.append(Contig(polisher.polishSmallSegment(base.asSegment(), c).seg_from.Seq(), str(len(split_contigs))))
        split_reads.append([al.seg_from.contig for al in c])
    maxpi = 1
    for i in range(mult):
        for j in range(mult):
            if i == j:
                sys.stdout.write("1.0 ")
                continue
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
        res = splitSegKmeans(aligner, base.segment(i * min_contig_length, i * min_contig_length + min_contig_length), mult, all_reads_list)
        if res is not None:
             return res
    res = splitSegKmeans(aligner, base.asSegment().suffix(length=min(min_contig_length, len(seq))), mult, all_reads_list)
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
        reads_to_resolve[eid] = list(set(reads_to_resolve[eid]))
    print "Reading reads"
    res_reads = ContigStorage()
    res = open(os.path.join(dir, "reads.fasta"), "w")
    for read in SeqIO.parse_by_name(rf):
        if read.id in relevant_read_ids:
            res_reads.add(Contig(read.seq, read.id))
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
        print reads_to_resolve[eid]
        print map(str, repeat_reads)
        split_contigs = splitRepeat(aligner, graph.e[eid].seq, mult, repeat_reads, min_contig_length)
        if split_contigs is None:
            print "Failed to resove edge", eid, "Aborting"
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

