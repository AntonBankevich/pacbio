import sys
from typing import Dict, List, Generator

from alignment.align_tools import Aligner
from common import params, basic
from common.alignment_storage import AlignmentPiece, ReadCollection
from common.sequences import Segment, Contig, ContigStorage
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_storage import NewLineStorage


class PseudoAlignment:
    def __init__(self, seg_from, seg_to):
        # type: (Segment, Segment) -> None
        self.seg_from = seg_from
        self.seg_to = seg_to

    def __str__(self):
        return str(self.seg_from) + "->" + str(self.seg_to)

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

    def dump(self, f):
        for key in self.storage:
            for al in self.storage[key]:
                f.write(str(al) + "\n")

class PairwiseReadRecruiter:
    def __init__(self, aligner, reads, lines):
        # type: (Aligner, ReadCollection, NewLineStorage) -> None
        sys.stdout.info("Preparing read overlaps for", len(reads), "reads")
        self.aligner = aligner
        self.reads = reads
        read_set = set()
        for line in lines:
            sys.stdout.trace(line.completely_resolved)
            for al in line.read_alignments:
                seg = al.seg_to.expand(params.k)
                sys.stdout.trace(al, line.completely_resolved.interSize(seg))
                if line.completely_resolved.interSize(seg) == len(seg) and \
                        al.seg_from.left < params.bad_end_length and \
                        al.rc.seg_from.left < params.bad_end_length and \
                        al.seg_to.left > params.k and \
                        al.rc.seg_to.left > params.k:
                    read_set.add(al.seg_from.contig.id)
        sys.stdout.info(len(read_set) / 2, "reads filtered out")
        reads = [read for read in reads if read.id not in read_set]
        sys.stdout.info(len(reads) / 2, "reads left for pairwise alignment")
        self.als = PairwiseStorage()
        for al in self.aligner.pairwiseAlign(reads):
            if al.__len__() > params.k:
                self.als.addAlignment(al.seg_from, al.seg_to)
        sys.stdout.info("Finished read overlap collection")

    def getRelevantAlignments(self, seg, min_overlap):
        # type: (Segment, int) -> Generator[AlignmentPiece]
        sys.stdout.trace("Requesting read alignments for", seg, " using palignments")
        line = seg.contig #type: NewLine
        reads = ContigStorage()
        print "Using reads ", line.read_alignments.allInter(seg, min_overlap)
        for base_read_al in line.read_alignments.allInter(seg, min_overlap):
            for read in self.als.getAlignments(base_read_al.seg_from.contig.id, params.k):
                reads.add(read)
        cnt = 0
        for al in self.aligner.localAlign(reads, ContigStorage([seg.contig])):
            if al.seg_to.interSize(seg) > min_overlap and al.__len__() > params.k:
                yield al
                cnt += 1
        sys.stdout.trace("Request for read alignments for", seg, "yielded", cnt, "alignments")



