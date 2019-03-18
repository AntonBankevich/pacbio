import itertools

from typing import Optional, Tuple, List, Dict, Iterable

from common import basic
from common.line_align import Scorer
from common.sequences import Segment, ReadCollection
from common.alignment_storage import AlignmentPiece, AlignedRead
from disjointig_resolve.accurate_line import NewLine, LinePosition
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.smart_storage import SegmentStorage

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
        self.scorer = Scorer()

    def tryExtend(self, line):
        # type: (NewLine) -> int
        line.completely_resolved.mergeSegments(k)
        bound = LinePosition(line, line.left())
        new_recruits = 0
        while True:
            bounds = [seg.left + k for seg in line.completely_resolved.items][1:]
            bounds.append(len(line))
            seg_to_resolve = None
            for seg, new_bound in zip(line.completely_resolved, bounds): #type: Segment, int
                if seg.right > bound:
                    seg_to_resolve = seg
                    bound = LinePosition(line, new_bound)
                    break
            if seg_to_resolve is None:
                break
            while True:
                result = self.attemptCleanResolution(seg_to_resolve)
                to_correct = [seg for seg, num in result if num > 0]
                to_correct = sorted(to_correct, key = lambda seg: (basic.Normalize(seg.contig.id), seg.left))
                for line, it in itertools.groupby(to_correct, key = lambda seg: seg.contig): # type: NewLine, Iterable[Segment]
                    to_polysh = []
                    for seg in it:
                        if seg.contig != line:
                           to_polysh.append(seg.RC())
                        else:
                            to_polysh.append(seg)
                    line.polyshSegments(to_polysh)
                total = sum([num for seg, num in result])
                new_recruits += total
                if total == 0:
                    break
        return new_recruits


    def attemptCleanResolution(self, resolved):
        # type: (Segment) -> List[Tuple[Segment, int]]
        # Find all lines that align to at leasr k nucls of resolved segment. Since this segment is resolve we get all
        line_alignments = self.alignedLines(resolved) # type: List[AlignmentPiece]
        #Find all reads that align to at least k nucls of resolved segment or corresponding segments on other lines
        reads = self.relevantReadAlignments(resolved, line_alignments) # type: List[AlignmentPiece]
        #Find all reads that align to at least k nucls of resolved segment or corresponding segments on other lines
        recruited_reads = set()
        for al in line_alignments:
            seg = al.seg_from
            line = seg.contig # type: NewLine
            for al in line.getReads(seg):
                recruited_reads.add(al.seg_from.contig.id)
        reads = filter(lambda al: al.seg_from.contig.id not in recruited_reads, reads)
        # For each read we find all its alignmens that can compete with alignments to this line
        read_alignments = self.generateReadToLineAlignments(reads, line_alignments) # type: List[List[AlignmentPiece]]
        new_recruits = 0
        for als in read_alignments:
            winner = self.tournament(als) #type: AlignmentPiece
            if winner is not None:
                line = winner.seg_to.contig # type: NewLine
                line.addReadAlignment(winner)
                if line == resolved.contig:
                    new_recruits += 1
        return new_recruits

    def compare(self, c1, c2):
        # type: (AlignmentPiece, AlignmentPiece) -> Tuple[Optional[int], Optional[int], Optional[int]]

        return None, None, None


    def fight(self, c1, c2):
        # type: (AlignmentPiece, AlignmentPiece) -> Optional[AlignmentPiece]
        assert c1.seg_from.contig == c2.seg_from.contig
        s1, s2, s12 = self.compare(c1, c2)
        winner = None
        if s12 is None:
            if s1 is None:
                winner = c2
            else:
                winner = c1
        else:
            if s12 < 25 or (s12 < 100 and abs(s1 - s2) < s12 * 0.8) or abs(s1 - s2) < s12 * 0.65:
                winner = None
            elif s1 > s2:
                winner = c2
            else:
                winner = c1
        if winner is None:
            print "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "No winner"
        else:
            print "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "Winner:", c2
        return winner

    def tournament(self, candidates):
        # type: (list[AlignmentPiece]) -> Optional[AlignmentPiece]
        best = None
        for candidate in candidates:
            if best is None:
                best = candidate
            else:
                best = self.fight(candidate, best)
        if best is None:
            return None
        if len(candidates) > 2:
            for candidate in candidates:
                if candidate == best:
                    continue
                fight_results = self.fight(candidate, best)
                if fight_results is None or fight_results != best:
                    return None
        return best
