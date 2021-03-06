import sys
import os

from disjointig_resolve.initialization import constructDisjointigs

sys.path.append("py")
from typing import List, BinaryIO

from alignment.align_tools import DirDistributor, Aligner
from common import basic, params
from common.alignment_storage import AlignmentPiece
from common.sequences import ContigStorage, Contig

from common.basic import CreateLog


def alsToReads(als):
    # type: (List[AlignmentPiece]) -> ContigStorage
    readIds = set()
    res = ContigStorage()
    for al in als:
        if al.seg_from.contig.id in readIds:
            continue
        readIds.add(al.seg_from.contig.id)
        res.add(al.seg_from.contig)
    return res

def fakeGraph(contigs, handler):
    # type: (ContigStorage, BinaryIO) -> None
    handler.writelines(["digraph", "{"])
    for contig in contigs:
        handler.writelines(['"0" -> "1"[label = "id ' + contig.id + ' k 100x", color = "black"];'])
    handler.writelines(["}"])

def recruit(seqs, reads, k, dir):
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    params.k = k
    relevant_reads = ContigStorage()
    disjointigs = seqs
    for i in range(2):
        sys.stdout.info("Recruiting iteration", i)
        als = filter(lambda al: len(al) > k, aligner.localAlign(reads, disjointigs))
        print len(als), "alignments"
        relevant_reads = alsToReads(als)
        l = sum(map(len, seqs.unique()))
        disjointigs = constructDisjointigs(relevant_reads, l, dd.nextDir())
        print len(disjointigs), "disjointigs"
        print disjointigs
    disjointigs.writeToFasta(open(os.path.join(dir, "disjointigs.fasta"), "w"))
    relevant_reads.writeToFasta(open(os.path.join(dir, "reads.fasta"), "w"))
    sys.stdout.info("Aligning repeat sequences to disjointigs")
    als = list(aligner.localAlign(seqs, disjointigs))
    print "\n".join(map(str, als))
    starts = dict()
    for dis in disjointigs:
        starts[dis.id] = len(dis)
    for al in als:
        if len(al) > k:
            starts[al.seg_to.contig.id] = min(starts[al.seg_to.contig.id], al.seg_to.left)
            al = al.rc
            starts[al.seg_to.contig.id] = min(starts[al.seg_to.contig.id], al.seg_to.left)
    print "Starts:"
    for cid, val in starts.items():
        print cid, val
    contigs = ContigStorage()
    cnt = 1
    for dis in disjointigs:
        if starts[dis.id] > k and starts[dis.id] < len(dis):
            print cnt, dis.id, starts[dis.id]
            contigs.add(Contig(dis.prefix(starts[dis.id]).Seq(), str(cnt)))
            cnt += 1
    for dis in disjointigs.unique():
        if len(dis) > k and starts[dis.id] == len(dis):
            print cnt, dis.id
            contigs.add(Contig(dis.seq, str(cnt)))
            cnt += 1
    contigs.writeToFasta(open(os.path.join(dir, "contigs.fasta"), "w"))
    fakeGraph(contigs, open(os.path.join(dir, "graph.gv"), "w"))




def main(args):
    dir = args[4]
    basic.ensure_dir_existance(dir)
    CreateLog(dir)
    sys.stdout.info("Starting graph-free recruitment")
    print " ".join(args)
    sys.stdout.info("Loading repeat sequences")
    seqs = ContigStorage().loadFromFasta(open(args[1], "r"), False)
    sys.stdout.info("Loading reads")
    reads = ContigStorage().loadFromFasta(open(args[2], "r"), False)
    k = int(args[3])
    recruit(seqs, reads, k, dir)
    sys.stdout.info("Finised graph-free recruitment")


if __name__ == "__main__":
    main(sys.argv)

