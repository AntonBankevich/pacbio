import os
import sys
sys.path.append("py")
from typing import List

from alignment.align_tools import Aligner, DirDistributor
from common import SeqIO
from common.alignment_storage import AlignmentPiece
from common.basic import CreateLog
from common.sequences import ContigStorage, Contig


def largestSubseq(als):
    # type: (List[AlignmentPiece]) -> List[AlignmentPiece]
    best = []#type: List[int]
    res = []#type: List[List[AlignmentPiece]]
    all_best = 0
    all_res = []
    for al in als:
        cbest = 0
        cres = []
        for val, r, in zip(best, res):
            if (len(r) == 0 or (r[-1].seg_from.right <= al.seg_from.left and r[-1].seg_to.right <= al.seg_to.left)) and val > cbest:
                cbest = val
                cres = list(r)
        cbest += len(al)
        cres.append(al)
        res.append(cres)
        best.append(cbest)
        if cbest > all_best:
            all_best = cbest
            all_res = cres
    print all_res
    return all_res



def iter_align(aligner, contig1, contig2):
    als = sorted(aligner.localAlign([contig1], ContigStorage([contig2])), key = lambda al: al.seg_from.left)
    als = [al for al in als if len(al) > 400 and al.seg_from.contig == contig1 and al.seg_to.contig == contig2]
    als = largestSubseq(als)
    res = []
    if len(als) > 0:
        for al1, al2 in zip(als[:-1], als[1:]):
            res.append(al1)
            if al1.seg_from.dist(al2.seg_from)> 400 and al1.seg_from.dist(al2.seg_from) > 400:
                seg1 = contig1.segment(al1.seg_from.right, al2.seg_from.left)
                seg2 = contig2.segment(al1.seg_to.right, al2.seg_to.left)
                tmp = iter_align(aligner, seg1.asContig(), seg2.asContig())
                for al in tmp:
                    res.append(al.queryAsSegment(seg1).targetAsSegment(seg2))
        res.append(als[-1])
    return res

def align(dir, contigs_file):
    CreateLog(dir)
    contigs = list(SeqIO.parse_fasta(open(contigs_file, "r")))
    assert len(contigs) == 2
    contigs = [Contig(contigs[0].seq, contigs[0].id), Contig(contigs[1].seq, contigs[1].id)]
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    als = iter_align(aligner, contigs[0], contigs[1])
    print "\n".join(map(str, als))


if __name__ == "__main__":
    dir = sys.argv[1]
    contigs_file = sys.argv[2]
    align(dir, contigs_file)