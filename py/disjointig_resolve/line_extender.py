import itertools

from typing import Optional, Tuple, List, Dict, Iterable

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import basic, params
from common.line_align import Scorer
from common.sequences import Segment, ReadCollection
from common.alignment_storage import AlignmentPiece, AlignedRead
from disjointig_resolve.accurate_line import NewLine, LinePosition
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage

k = 1000

# class LineCorrector:
#     def __init__(self):
#         self.alignments = [] # type: List[AlignmentPiece]
#         self.completely_restored = [] # type: List[Segment]
#         self.correct_segments = [] # type: List[Segment]
#
#     def addAlignment(self, al):
#         # type: (AlignmentPiece) -> None
#         self.alignments.append(al)
#
#     def addRestored(self, seg):
#         # type: (Segment) -> None
#         self.completely_restored.append(seg)
#
#     def addCorrect(self, seg):
#         # type: (Segment) -> None
#         self.correct_segments.append(seg)
#
#     def dumpRestored(self):
#         lines = dict() # type: Dict[str, NewLine]
#         for seg in self.completely_restored:
#             line = seg.contig # type: NewLine
#             lines[line.id] = line
#             line.completely_resolved.add(seg)
#         for line in lines.values():
#             line.completely_resolved.mergeSegments(k)
#         self.completely_restored = []
#
#     def dumpCorrect(self):
#         lines = dict() # type: Dict[str, NewLine]
#         for seg in self.correct_segments:
#             line = seg.contig # type: NewLine
#             lines[line.id] = line
#             line.correct_segments.add(seg)
#         for line in lines.values():
#             line.correct_segments.mergeSegments()
#         self.completely_restored = []
#
#     def dumpAlignments(self):
#         self.alignments = sorted(self.alignments, key = lambda al: al.seg_to.contig.id)
#         for line, it in itertools.groupby(self.alignments, key = lambda al: al.seg_to.contig):
#             line.correctSequence(list(it))
#         self.alignments = []
#
#     def dump(self):
#         self.dumpRestored()
#         self.dumpCorrect()
#         self.dumpAlignments()


# LineExtender contains algorithms for resolving repeats by extending lines.
# The results of its work is reflected in lines themselves as they grow and knot to each other
class LineExtender:
    def __init__(self, aligner, knotter, disjointigs, dot_plot):
        # type: (Aligner, Knotter, DisjointigCollection, LineDotPlot) -> None
        self.aligner = aligner
        self.polisher = Polisher(aligner, aligner.dir_distributor)
        self.knotter = knotter
        self.disjointigs = disjointigs
        self.dot_plot = dot_plot
        self.scorer = Scorer()

    def tryExtend(self, line):
        # type: (NewLine) -> int
        line.completely_resolved.mergeSegments(k)
        bound = LinePosition(line, line.left())
        new_recruits = 0
        while True:
            seg_to_resolve = line.completely_resolved.find(bound.suffix(), k)
            if seg_to_resolve is None:
                break
            result = self.attemptCleanResolution(seg_to_resolve)
            total = sum([num for seg, num in result])
            new_recruits += total
            if total == 0:
                bound = seg_to_resolve.right - k + 1
                continue
            to_correct = [seg for seg, num in result if num > 0]
            self.updateAllStructures(to_correct)
            if seg_to_resolve.right > len(line) - 2000:
                self.attemptExtend(line)
                if knotter.tryKnot(line) is not None:
                    return new_recruits
        return new_recruits

    def updateAllStructures(self, to_correct):
        to_correct = sorted(to_correct, key=lambda seg: (basic.Normalize(seg.contig.id), seg.left))
        for line, it in itertools.groupby(to_correct,
                                               key=lambda seg: seg.contig):  # type: NewLine, Iterable[Segment]
            to_polysh = []
            for seg in it:
                if seg.contig != line:
                    to_polysh.append(seg.RC())
                else:
                    to_polysh.append(seg)
            self.polyshSegments(line, to_polysh)
            self.updateCorrectSegments(line)
            self.updateCompletelyResolved(line)  # IMPLEMENT

    def attemptCleanResolution(self, resolved):
        # type: (Segment) -> List[Tuple[Segment, int]]
        # Find all lines that align to at least k nucls of resolved segment. Since this segment is resolve we get all
        resolved = resolved.suffix(length = min(len(resolved), k))
        line_alignments = filter(lambda al: len(al.seg_to) >= k, self.dot_plot.getAlignmentsTo(resolved)) # type: List[AlignmentPiece]
        line_alignments = [al.reduce(target=resolved) for al in line_alignments]
        read_alignments = [] # type: List[Tuple[AlignmentPiece, Segment]]
        correct_segments = []
        for ltl in line_alignments:
            line = ltl.seg_from.contig # type: NewLine
            correct_segments.append(line.correct_segments.find(ltl.seg_from))
            assert correct_segments[-1] is not None and correct_segments[-1].contains(ltl.seg_from)
            read_alignments.extend(zip(line.getPotentialAlignmentsTo(ltl.seg_from), itertools.cycle([correct_segments[-1]])))
        read_alignments = sorted(read_alignments, key=lambda al: al[0].seg_from.contig.name)
        # removing all reads that are already sorted to one of the contigs
        alignments_by_read = itertools.groupby(lambda al: al.seg_from.contig.name, read_alignments)
        new_recruits = 0
        # TODO: parallel
        for name, it in alignments_by_read:
            als = list(it) # type: List[Tuple[AlignmentPiece, Segment]]
            read = als[0][0].seg_from.contig # type: AlignedRead
            skip = False
            for al1 in als:
                for al2 in read.alignments:
                    if al1[0].seg_to.inter(al2.seg_to):
                        skip = True
                        break
                if skip:
                    break
            if skip:
                continue
            winner = self.tournament(als) #type: AlignmentPiece
            if winner is not None:
                line = winner.seg_to.contig # type: NewLine
                line.addReadAlignment(winner)
                if line == resolved.contig:
                    new_recruits += 1
        return new_recruits

    def fight(self, c1, c2):
        # type: (Tuple[AlignmentPiece, Segment], Tuple[AlignmentPiece, Segment]) -> Optional[AlignmentPiece]
        assert c1[0].seg_from.contig == c2[0].seg_from.contig
        s1, s2, s12 = self.scorer.scoreInCorrectSegments(c1[0], c1[1], c2[0], c2[1])
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
        # type: (list[Tuple[AlignmentPiece, Segment]]) -> Optional[AlignmentPiece]
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

    def attemptExtend(self, line):
        # type: (NewLine) -> None
        new_contig, relevant_als = self.polisher.polishEnd(list(line.read_alignments.allInter(line.suffix(1000))))
        if len(new_contig) == len(line):
            return
        assert line.seq == new_contig.prefix(len=len(line)).Seq()
        line.extendRight(new_contig.suffix(pos = len(line)).Seq(), relevant_als)

    def polyshSegments(self, line, to_polysh):
        # type: (NewLine, List[Segment]) -> None
        segs = SegmentStorage()
        corrections = AlignmentStorage()
        line.addListener(segs)
        segs.addAll(to_polysh)
        segs.sort()
        for i in range(len(segs)):
            corrections.add(self.polisher.polishSegment(segs[i], list(line.read_alignments.allInter(segs[i]))))
        line.correctSequence(corrections.items)
        line.removeListener(segs)

    def updateCorrectSegments(self, line):
        # type: (NewLine) -> None
        line.updateCorrectSegments(line.asSegment())