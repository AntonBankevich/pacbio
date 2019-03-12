import itertools

from typing import Optional, Tuple, List, Dict

from common.sequences import Segment
from common.alignment_storage import AlignmentPiece
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.disjointigs import DisjointigCollection

k = 1000

class LineCorrector:
    def __init__(self):
        self.alignments = [] # type: List[AlignmentPiece]
        self.completely_restored = [] # type: List[Segment]
        self.correct_segments = [] # type: List[Segment]

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.alignments.append(al)

    def addRestored(self, seg):
        # type: (Segment) -> None
        self.completely_restored.append(seg)

    def addCorrect(self, seg):
        # type: (Segment) -> None
        self.correct_segments.append(seg)

    def dumpRestored(self):
        lines = dict() # type: Dict[str, NewLine]
        for seg in self.completely_restored:
            line = seg.contig # type: NewLine
            lines[line.id] = line
            line.completely_resolved.add(seg)
        for line in lines.values():
            line.completely_resolved.mergeSegments(k)
        self.completely_restored = []

    def dumpCorrect(self):
        lines = dict() # type: Dict[str, NewLine]
        for seg in self.correct_segments:
            line = seg.contig # type: NewLine
            lines[line.id] = line
            line.correct_segments.add(seg)
        for line in lines.values():
            line.correct_segments.mergeSegments()
        self.completely_restored = []

    def dumpAlignments(self):
        self.alignments = sorted(self.alignments, key = lambda al: al.seg_to.contig.id)
        for line, it in itertools.groupby(self.alignments, key = lambda al: al.seg_to.contig):
            line.correctSequence(list(it))
        self.alignments = []

    def dump(self):
        self.dumpRestored()
        self.dumpCorrect()
        self.dumpAlignments()


class LineExtender:
    def __init__(self, disjointigs):
        # type: (DisjointigCollection) -> None
        self.disjointigs = disjointigs

    def tryExtend(self, line):
        # type: (NewLine) -> bool
        line.completely_resolved.mergeSegments(k)
        corrector = LineCorrector()
        changed = False
        for seg, bound in zip(line.completely_resolved, [seg.left + k for seg in line.completely_resolved[1:]].extend([len(line)])):
             changed |= self.attemptCleanResolution(seg, bound, corrector)
        corrector.dump()
        return changed


    def attemptCleanResolution(self, resolved, bound, corrector):
        # type: (Segment, int, LineCorrector) -> bool
        # IMPLEMENT!
        pass