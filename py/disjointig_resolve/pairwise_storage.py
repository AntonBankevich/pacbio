import sys
from typing import Dict, List, Generator

from alignment.align_tools import Aligner
from common import params, basic
from common.alignment_storage import AlignmentPiece, ReadCollection
from common.sequences import Segment, Contig, ContigStorage
from disjointig_resolve.accurate_line import NewLine


class PseudoAlignment:
    def __init__(self, seg_from, seg_to):
        # type: (Segment, Segment) -> None
        self.seg_from = seg_from
        self.seg_to = seg_to

class PairwiseStorage:
    def __init__(self):
        self.storage = dict() #type: Dict[str, List[PseudoAlignment]]

    def addAlignment(self, seg_from, seg_to):
        # type: (Segment, Segment) -> None
        self.innerAdd(seg_from, seg_to)
        self.innerAdd(seg_from.RC(), seg_to.RC())
        self.innerAdd(seg_to, seg_from)
        self.innerAdd(seg_to.RC(), seg_from.RC())

    def innerAdd(self, seg_from, seg_to):
        if seg_to.contig.id not in self.storage:
            self.storage[seg_to.contig.id] = [PseudoAlignment(seg_from, seg_to)]
        else:
            self.storage[seg_to.contig.id].append(PseudoAlignment(seg_from, seg_to))

    def getAlignments(self, contigid, min_overlap):
        # type: (str, int) -> Generator[Contig]
        if contigid not in self.storage:
            return
        for al in self.storage[contigid]:
            yield al.seg_from.contig

class PairwiseReadRecruiter:
    def __init__(self, aligner, reads):
        # type: (Aligner, ReadCollection) -> None
        self.aligner = aligner
        self.reads = reads
        self.als = PairwiseStorage()
        for al in self.aligner.pairwiseAlign(reads):
            if al.__len__() > params.k:
                self.als.addAlignment(al.seg_from, al.seg_to)

    def getRelevantAlignments(self, seg, min_overlap):
        # type: (Segment, int) -> Generator[AlignmentPiece]
        sys.stdout.trace("Requesting read alignments for", seg, " using palignments")
        line = seg.contig #type: NewLine
        reads = ContigStorage()
        for base_read_al in line.read_alignments.allInter(seg, min_overlap):
            sys.stdout.write("Base " + str(base_read_al))
            for read in self.als.getAlignments(base_read_al.seg_from.contig.id, params.k):
                reads.add(read)
                sys.stdout.write(" " + read.id)
            sys.stdout.write("\n")
        for al in self.aligner.localAlign(reads, ContigStorage([seg.contig])):
            if al.seg_to.interSize(seg) > min_overlap and al.__len__() > params.k:
                print "Relevant read", al
                yield al
        sys.stdout.trace("Request for read alignments for", seg, "finished")



