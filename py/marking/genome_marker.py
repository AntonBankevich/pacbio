from typing import Dict, List

from alignment.align_tools import Aligner
from common import basic
from common.sequences import ContigStorage, Segment
from disjointig_resolve.smart_storage import SegmentStorage


def markUniqueInGenome(aligner, genome, k):
    # type: (Aligner, ContigStorage, int) -> Dict[str, SegmentStorage]
    badsegs = dict() #type: Dict[str, SegmentStorage]
    res = dict()
    for contig in genome:
        badsegs[contig.id] = SegmentStorage()
    for al in aligner.localAlign(genome, genome):
        if len(al.seg_to) >= k:
            badsegs[al.seg_to.contig.id].add(al.seg_to)
    for contig in genome.unique():
        ss = badsegs[contig.id]
        ss.mergeSegments(k)
        res[contig.id] = ss.reverse(contig, k - 1)
    return res

def checkUnique(aligner, seqs, genome, k):
    # type: (Aligner, ContigStorage, ContigStorage, int) -> None
    ethalon = markUniqueInGenome(aligner, genome, k)
    for ss in ethalon.values():
        print ss
    als = dict()
    for al in aligner.overlapAlign(seqs, genome):
        if not basic.isCanonocal(al.seg_from.contig.id):
            al = al.rc
        if al.seg_from.left < 10 and al.rc.seg_from.left < 10:
            if al.seg_from.contig.id not in als or al.percentIdentity() > als[al.seg_from.contig.id].percentIdentity():
                als[al.seg_from.contig.id] = al
    for contig in seqs.unique():
        if contig.id not in als:
            print "No alignment to genome for candidate:", contig
        else:
            al = als[contig.id]
            if ethalon[al.seg_to.contig.id].interSize(al.seg_to) == len(al.seg_to):
                print "Confirmed candidate:", contig, al
            else:
                print "Failed candidate:", contig, al