import subprocess
import sys
import os
sys.path.append("py")
from typing import List, BinaryIO

from alignment.align_tools import DirDistributor, Aligner
from common import basic, params
from common.alignment_storage import AlignmentPiece
from common.sequences import ContigCollection, Contig

from common.basic import CreateLog

def constructDisjointigs(reads, total_length, dir):
    # type: (ContigCollection, int, str) -> ContigCollection
    reads_file = os.path.join(dir, "reads.fasta")
    disjointigs_file = os.path.join(dir, "disjointigs.fasta")
    log_file = os.path.join(dir, "log.txt")
    reads.writeToFasta(open(reads_file, "w"))
    subprocess.check_call(["./bin/flye-assemble", reads_file, disjointigs_file, str(total_length), "flye/config/bin_cfg/asm_raw_reads.cfg", "-v", "1500", "-t", "16", "-u", "-l", log_file])
    return ContigCollection().loadFromFasta(open(disjointigs_file, "r"), False)

def alsToReads(als):
    # type: (List[AlignmentPiece]) -> ContigCollection
    readIds = set()
    res = ContigCollection()
    for al in als:
        if al.seg_from.contig.id in readIds:
            continue
        readIds.add(al.seg_from.contig.id)
        res.add(al.seg_from.contig)
    return res

def fakeGraph(contigs, handler):
    # type: (ContigCollection, BinaryIO) -> None
    handler.writelines(["digraph", "{"])
    for contig in contigs:
        handler.writelines(['"0" -> "1"[label = "id ' + contig.id + ' k 100x", color = "black"];'])
    handler.writelines(["}"])

def recruit(seqs, reads, k, dir):
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    params.k = k
    relevant_reads = ContigCollection()
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
    als = aligner.localAlign(disjointigs, seqs)
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
    contigs = ContigCollection()
    cnt = 1
    for dis in disjointigs:
        if starts[dis.id] > k and starts[dis.id] < len(dis):
            print cnt, dis.id, starts[dis.id]
            contigs.add(Contig(dis.prefix(starts[dis.id]).Seq(), str(cnt)))
            cnt += 1
    for dis in disjointigs.unique():
        if len(dis) > k:
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
    seqs = ContigCollection().loadFromFasta(open(args[1], "r"), False)
    sys.stdout.info("Loading reads")
    reads = ContigCollection().loadFromFasta(open(args[2], "r"), False)
    k = int(args[3])
    recruit(seqs, reads, k, dir)
    sys.stdout.info("Finised graph-free recruitment")


if __name__ == "__main__":
    main(sys.argv)

