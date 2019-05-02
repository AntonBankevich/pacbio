import itertools

from typing import Optional, Tuple, List, Dict, Iterable, Set, Generator, Union, Iterator

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import basic, params
from common.line_align import Scorer
from common.sequences import Segment
from common.alignment_storage import AlignmentPiece, AlignedRead, ReadCollection
from disjointig_resolve.accurate_line import NewLine, LinePosition
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.knotter import LineKnotter
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage

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
        # type: (Aligner, LineKnotter, DisjointigCollection, LineDotPlot) -> None
        self.aligner = aligner
        self.polisher = Polisher(aligner, aligner.dir_distributor)
        self.knotter = knotter
        self.disjointigs = disjointigs
        self.dot_plot = dot_plot
        self.scorer = Scorer()

    def tryExtend(self, line):
        # type: (NewLine) -> int
        line.completely_resolved.mergeSegments(params.k)
        bound = LinePosition(line, line.left())
        new_recruits = 0
        while True:
            seg_to_resolve = line.completely_resolved.find(bound.suffix(), params.k)
            if seg_to_resolve is None:
                break
            result = self.attemptCleanResolution(seg_to_resolve)
            total = sum([len(arr) for seg, arr in result])
            new_recruits += total
            if total == 0:
                bound = LinePosition(line, seg_to_resolve.right - params.k + 1)
                continue
            self.updateAllStructures([seg for seg, arr in result])
            if seg_to_resolve.right > len(line) - 1000:
                if self.attemptExtend(line):
                    self.updateAllStructures([seg_to_resolve])
                    if self.knotter.tryKnotRight(line) is not None:
                        return new_recruits
        return new_recruits

    # input: a collection of segments that had reads recruited to.
    def updateAllStructures(self, interesting_segments):
        # Correct read sequences, update correct segment storages. Return segments that were corrected.
        corrected = self.correctSequences(interesting_segments)
        # Collect all relevant contig segments, collect all reads that align to relevant segments. Mark resolved bound for each read.
        records = self.collectRecords(corrected)
        # Update resolced segmetns on all relevant contig positions
        self.updateResolved(records)

    def updateResolved(self, records):
        # type: (List[Record]) -> None
        ok = True
        while ok:
            ok = False
            for rec in records:
                if self.attemptProlongResolved(rec):
                    ok = True
                    print "prolonged"
            if not ok:
                for rec in records:
                    if self.attemptJump(rec):
                        ok = True
                        print "jumped"
        for rec in records:
            line = rec.resolved.contig  # type: NewLine
            line.completely_resolved.add(rec.resolved)
            for seg in rec.old_resolved:
                line.completely_resolved.add(seg)
            line.completely_resolved.mergeSegments(params.k)

    def collectRecords(self, corrected):
        # type: (List[Segment]) -> List[Record]
        read_bounds = dict()
        records = dict()  # type: Dict[Segment, LineExtender.Record]
        good_reads = set()
        for seg in corrected:
            for al in self.dot_plot.allInter(seg):
                line = al.seg_from.contig  # type:NewLine
                for seg_correct in line.correct_segments.allInter(al.seg_from):
                    for seg_resolved in line.completely_resolved.allInter(seg_correct):
                        next = line.completely_resolved.find(line.asSegment().suffix(pos=seg_resolved.right), 1)
                        if next is None:
                            next_start = len(line)
                        else:
                            next_start = next.left
                        records[seg_resolved] = self.createRecord(seg_resolved, next_start, seg_correct, good_reads, read_bounds)
        records = list(records.values())  # type: List[LineExtender.Record]
        return records

    def correctSequences(self, interesting_segments):
        to_correct = []
        for seg in interesting_segments:
            line = seg.contig # type: NewLine
            correct = line.correct_segments.find(seg)
            next = line.correct_segments.find(line.suffix(correct.right), 1)
            if next is None:
                right = len(line)
            else:
                right = min(len(line), next.left + params.k / 2)
            to_correct.append(line.segment(correct.left, right))
        to_correct = sorted(to_correct, key=lambda seg: (basic.Normalize(seg.contig.id), seg.left))
        corrected = []
        for line_id, it in itertools.groupby(to_correct,
                                          key=lambda seg: basic.Normalize(seg.contig.id)):  # type: NewLine, Iterable[Segment]
            line = None # type: NewLine
            to_polysh = []
            for seg in it:
                if seg.contig.id != line_id:
                    to_polysh.append(seg.RC())
                    line = seg.contig.rc
                else:
                    to_polysh.append(seg)
                    line = seg.contig
            new_segments = self.polyshSegments(line, to_polysh)
            corrected.extend(new_segments)
            self.updateCorrectSegments(line)
        return corrected

    def attemptCleanResolution(self, resolved):
        # type: (Segment) -> List[Tuple[Segment, List[AlignmentPiece]]]
        # Find all lines that align to at least k nucls of resolved segment. Since this segment is resolve we get all
        print resolved
        resolved = resolved.suffix(length = min(len(resolved), params.k * 2))
        line_alignments = filter(lambda al: len(al.seg_to) >= params.k and resolved.interSize(al.seg_to) > params.k / 2,
                                 self.dot_plot.allInter(resolved)) # type: List[AlignmentPiece]
        line_alignments = [al.reduce(target=resolved) for al in line_alignments]
        read_alignments = [] # type: List[Tuple[AlignmentPiece, Segment]]
        correct_segments = []
        for ltl in line_alignments:
            line = ltl.seg_from.contig # type: NewLine
            correct_segments.append(line.correct_segments.find(ltl.seg_from))
            assert correct_segments[-1] is not None and correct_segments[-1].contains(ltl.seg_from)
            read_alignments.extend(zip(line.getRelevantAlignmentsFor(ltl.seg_from), itertools.cycle([correct_segments[-1]])))
        print list(self.dot_plot.allInter(resolved))
        print line_alignments
        print read_alignments
        read_alignments = sorted(read_alignments, key=lambda al: al[0].seg_from.contig.id)
        # removing all reads that are already sorted to one of the contigs
        alignments_by_read = itertools.groupby(read_alignments, lambda al: al[0].seg_from.contig.id)
        new_recruits = []
        # TODO: parallel
        for name, it in alignments_by_read:
            als = list(it) # type: List[Tuple[AlignmentPiece, Segment]]
            print als
            ok = False
            for al in als:
                if al[0].seg_to.interSize(resolved) >= params.k:
                    ok = True
                    break
            print ok
            if not ok:
                continue
            read = als[0][0].seg_from.contig # type: AlignedRead
            skip = False
            for al1 in als:
                for al2 in read.alignments:
                    if al1[0].seg_to.inter(al2.seg_to):
                        skip = True
                        break
                if skip:
                    break
            print skip
            if skip:
                continue
            winner, seg = self.tournament(als) #type: AlignmentPiece, Segment
            if winner is not None:
                line = winner.seg_to.contig # type: NewLine
                line.addReadAlignment(winner)
                new_recruits.append((seg, winner))
        new_recruits = sorted(new_recruits, key = lambda rec: (rec[0].contig.id, rec[0].left, rec[0].right))
        return [(seg, [al for seg, al in it]) for seg, it in itertools.groupby(new_recruits, key = lambda rec: rec[0])]

    def fight(self, c1, c2):
        # type: (Tuple[AlignmentPiece, Segment], Tuple[AlignmentPiece, Segment]) -> Optional[Tuple[AlignmentPiece, Segment]]
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
            print c1, c2, s1, s2, s12
            print "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "Winner:", winner
        return winner

    def tournament(self, candidates):
        # type: (list[Tuple[AlignmentPiece, Segment]]) -> Tuple[Optional[AlignmentPiece], Optional[Segment]]]
        best = None
        for candidate in candidates:
            if best is None:
                best = candidate
            else:
                best = self.fight(candidate, best)
        if best is None:
            return None, None
        if len(candidates) > 2:
            for candidate in candidates:
                if candidate == best:
                    continue
                fight_results = self.fight(candidate, best)
                if fight_results is None or fight_results != best:
                    return None, None
        return best

    def attemptExtend(self, line):
        # type: (NewLine) -> bool
        new_contig, relevant_als = self.polisher.polishEnd(list(line.read_alignments.allInter(line.suffix(1000))))
        if len(new_contig) == len(line):
            return False
        assert line.seq == new_contig.prefix(len=len(line)).Seq()
        line.extendRight(new_contig.suffix(pos = len(line)).Seq(), relevant_als)
        return True

    def polyshSegments(self, line, to_polysh):
        # type: (NewLine, List[Segment]) -> List[Segment]
        segs = SegmentStorage()
        corrections = AlignmentStorage()
        line.addListener(segs)
        segs.addAll(to_polysh)
        segs.sort()
        for seg in segs:
            corrections.add(self.polisher.polishSegment(seg, list(line.read_alignments.allInter(seg))))
        line.correctSequence(corrections.items)
        line.removeListener(segs)
        return list(segs)

    def updateCorrectSegments(self, line):
        # type: (NewLine) -> None
        line.updateCorrectSegments(line.asSegment())

    class Record:
        def __init__(self, resolved, next, correct, good_reads, read_bounds):
            # type: (Segment, int, Segment, Set[str], Dict[str, int]) -> None
            self.line = resolved.contig # type: NewLine
            self.resolved = resolved
            self.old_resolved = []
            self.next_resolved_start = next
            self.correct = correct
            self.good_reads = good_reads
            self.read_bounds = read_bounds
            self.reads = [] # type: List[AlignmentPiece]
            self.sorted = True
            self.potential_good = []

        def setResolved(self, seg):
            # type: (Segment) -> None
            if seg.interSize(self.resolved) >= params.k - 1:
                self.resolved = self.resolved.cup(seg)
            else:
                self.old_resolved.append(self.resolved)
                self.resolved = seg
            self.updateGood()

        def add(self, al):
            # type: (AlignmentPiece) -> None
            if al.seg_from.left < params.k / 2:
                self.potential_good.append(al)
            else:
                self.reads.append(al)
            read = al.seg_from.contig # type: AlignedRead
            if read.id not in self.read_bounds:
                self.read_bounds[read.id] = len(read)
            if al.rc.seg_to.left < 50:
                self.read_bounds[read.id] = min(self.read_bounds[read.id], al.seg_from.right)
            self.sorted = False

        def addAll(self, als):
            # type: (Iterator[AlignmentPiece]) -> None
            for al in als:
                self.add(al)

        def sort(self):
            if not self.sorted:
                self.reads = sorted(self.reads, key = lambda al: -al.seg_to.left)
                self.potential_good = sorted(self.potential_good, key = lambda al: -al.seg_to.left)
                self.sorted = True

        def get(self, num = None, right = None, min_inter = 0):
            # type: (int, Segment, int) -> List[AlignmentPiece]
            self.sort()
            if num is None:
                num = len(self.reads)
            if right is None:
                right = self.resolved.right
            popped = []
            res = []
            while len(res) < num and len(self.reads) > 0 and self.reads[-1].seg_to.left < right:
                al = self.reads.pop()
                if al.seg_from.contig.id not in self.good_reads or al.seg_from.right > self.read_bounds[al.seg_from.contig.id]:
                    popped.append(al)
                    if len(al.seg_to) >= min_inter:
                        res.append(al)
            self.reads.extend(popped)
            return res

        def __iter__(self):
            for al in self.reads[::-1]:
                if al.seg_from.contig.id not in self.good_reads:
                    yield al

        def updateGood(self):
            self.sort()
            while len(self.reads) > 0 and self.reads[-1].seg_to.left <= self.resolved.right - params.k:
                self.good_reads.add(self.reads.pop().seg_from.contig.id)
            while len(self.potential_good) > 0 and self.potential_good[-1].seg_to.left <= self.resolved.right - params.k:
                self.good_reads.add(self.potential_good.pop().seg_from.contig.id)

        Remove all reads that align before current resolved segment

        def pop(self):
            return self.reads.pop()

        def __str__(self):
            return str([self.resolved, self.correct, self.next_resolved_start, self.reads])

    def createRecord(self, resolved, next_start, correct, good_reads, read_bounds):
        # type: (Segment, int, Segment, Set[str], Dict[str, int]) -> Record
        line = resolved.contig # type: NewLine
        focus = line.segment(resolved.right - params.k, correct.right)
        # for al in self.dot_plot.allInter(focus):
        #     line1 = al.seg_from.contig # type: NewLine
        #     ok = False
        #     for seg in line1.correct_segments.allInter(al.seg_from): # type: Segment
        #         if seg.left <= al.seg_from.left:
        #             ok = True
        #             if seg.right < al.seg_from.right - 20:
        #                 next_start = min(next_start, al.matchingSequence(True).mapDown(seg.right))
        #     if not ok:
        #         next_start = min(next_start, al.seg_to.left + params.k / 2)
        res = self.Record(resolved, next_start, correct, good_reads, read_bounds)
        als = line.getRelevantAlignmentsFor(focus)
        # als = filter(lambda al: al.seg_from.left > params.k / 2 + 20, als)
        res.addAll(als)
        res.updateGood()
        return res

    def findResolvedBound(self, rec, inter_size):
        # type: (Record, int) -> int
        bad_reads = []
        for read in rec:
            if len(read.seg_to) >= inter_size:
                bad_reads.append(read)
            if len(bad_reads) >= 3:
                if bad_reads[-1].seg_to.left - bad_reads[0].seg_to.left > 50:
                    bad_reads = bad_reads[1:]
                else:
                    break
        if len(bad_reads) < 3:
            return len(rec.line)
        else:
            return bad_reads[0].seg_to.left

    def attemptProlongResolved(self, rec):
        # type: (Record) -> bool
        bound = min(rec.correct.right, rec.next_resolved_start + params.k, self.findResolvedBound(rec, params.k) + params.k / 2)
        if bound > rec.resolved.right:
            candidate = self.segmentsWithGoodCopies(rec.line.segment(rec.resolved.right - params.k, bound), params.k)[0]
            if candidate.right > rec.resolved.right:
                rec.resolved = rec.resolved.contig.segment(rec.resolved.left, candidate.right)
                rec.updateGood()
                return True
        return False

    # TODO: filter chimeric reads
    def attemptJump(self, rec):
        # type: (Record) -> bool
        bound = min(rec.correct.right, rec.next_resolved_start + params.k, self.findResolvedBound(rec, params.l) + params.k / 2)
        bad_segments = SegmentStorage()
        for al in rec:
            if al.seg_to.left > bound:
                break
            if al.seg_from.left > params.k / 2 and al.rc.seg_from.left > params.k / 2:
                bad_segments.add(al.seg_to)
        bad_segments.mergeSegments()
        good_segments = bad_segments.reverse().reduce(rec.line.segment(rec.resolved.right - params.k, bound))
        for seg in good_segments:
            for seg1 in self.segmentsWithGoodCopies(seg.expand(params.k / 2), params.k):
                if len(seg) >= params.k and seg.right > rec.resolved.right:
                    rec.setResolved(seg1)
                    return True
        return False

    def segmentsWithGoodCopies(self, seg, inter_size):
        # type: (Segment, int) -> List[Segment]
        als = self.dot_plot.allInter(seg)
        segs = SegmentStorage()
        for al in als:
            if len(al.seg_to) >= inter_size:
                line = al.seg_from.contig # type: NewLine
                correct = line.correct_segments.reverse().reduce(al.seg_from)
                matching = al.matchingSequence()
                for seg1 in correct:
                    segs.add(matching.mapSegDown(seg.contig, seg1))
        segs.mergeSegments()
        return list(segs.reverse().reduce(seg).filterBySize(min=inter_size))

