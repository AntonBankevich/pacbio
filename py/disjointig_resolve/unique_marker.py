import sys

from typing import Generator, Iterable, List, Tuple

from alignment.align_tools import Aligner
from common import params
from common.alignment_storage import AlignedRead, AlignmentPiece
from common.sequences import Segment
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.disjointigs import DisjointigCollection, Disjointig
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage


class UniqueMarker:
    def __init__(self, aligner):
        # type: (Aligner) -> None
        self.aligner = aligner

    # Mark unique regions on a contig as correct
    def findUniqueInDisjointigs(self, disjointigs):
        # type: (DisjointigCollection) -> Generator[Segment]
        # TODO: find unique sequences
        pass

    def link(self, line, arr, dist):
        # type: (NewLine, Iterable[int], int) -> List[Tuple[Segment, int]]
        arr = sorted(arr)
        if len(arr) == 0:
            return []
        res = []
        left = arr[0]
        prev = arr[0]
        cur = 1
        for pos in arr[1:]:
            if pos > prev + dist:
                res.append((line.segment(left, prev), cur))
                left = pos
                prev = pos
                cur = 1
            else:
                prev = pos
                cur += 1
        res.append((line.segment(left, prev), cur))
        return res

    def markUniqueInLine(self, line):
        # type: (NewLine) -> None
        sys.stdout.info("Finding unique in", line)
        alignments = list(line.read_alignments) # type: List[AlignmentPiece]
        # for i in range(len(alignments)):
        #     for j in range(len(alignments)):
        #         if i != j and alignments[j] is not None and alignments[j].contains(alignments[i]):
        #             print "Removed alignment", alignments[i], alignments[j]
        #             print "\n".join(alignments[j].asMatchingStrings())
        #             alignments[i] = None
        #             break
        # alignments = filter(lambda al: al is not None, alignments)
        alignments = sorted(alignments, key=lambda al:al.seg_to.left)
        sys.stdout.trace("Sorting finished")
        inc = self.link(line, [al.seg_to.left for al in alignments if al.seg_from.left > 1000 and al.seg_to.left > 50], 20)
        # for al in alignments:
        #     if al.seg_from.left > 1000 and al.seg_to.left > 50:
        #         print "Breakpoint read inc:", al
        inc.append((line.segment(len(line) - 1, len(line)), params.min_k_mer_cov))
        alignments = sorted(alignments, key=lambda al:al.seg_to.right)
        out = self.link(line, [al.seg_to.right for al in alignments if al.rc.seg_from.left > 1000 and al.rc.seg_to.left > 50 ], 20)
        sys.stdout.trace("Linking finished")
        # for al in alignments:
        #     if al.rc.seg_from.left > 1000 and al.rc.seg_to.left > 50:
        #         print "Breakpoint read out:", al
        out.insert(0, (line.segment(0, 1), params.min_k_mer_cov))
        print "inc:", inc
        print "out:", out
        events = []
        for seg, val in inc:
            if val >= params.min_k_mer_cov:
                events.append((seg.right, -1))
        for seg, val in out:
            if val >= params.min_k_mer_cov:
                events.append((seg.left, 1))
        events= sorted(events)
        sys.stdout.trace("Events collected and sorted", len(events))
        print events
        segs = SegmentStorage()
        for e1, e2 in zip(events[:-1], events[1:]):
            seg = line.segment(e1[0], e2[0])
            if e1[1] == 1 and e2[1] == -1:
                if len(seg) > params.max_allowed_unaligned:
                    seg = seg.expand(params.k / 2).expandToSize(params.k + 50)
                    if len(seg) >= params.k:
                        segs.add(seg)
            elif len(seg) > 50000:
                segs.add(seg.shrink(3000))
        sys.stdout.trace("Unique segmetns selected")
        # inc = SegmentStorage().addAll([seg for seg, cov in inc if cov >= params.min_k_mer_cov]).reverse(line)
        # out = SegmentStorage().addAll([seg for seg, cov in out if cov >= params.min_k_mer_cov]).reverse(line)
        # print "inc:", inc
        # print "out:", out
        # segs1 = inc.orderedCap(out)
        # print "segs1:", segs1
        # # segs = segs1.cap(segs2).expand(params.k / 2).filterBySize(min = params.k)
        # segs = segs1.filterBySize(min = params.max_allowed_unaligned).expand(params.k / 2).expandTo(params.k + 50).filterBySize(min=params.k)
        # # segs = segs1.expand(params.k / 2).filterBySize(min = params.k)
        line.cleanReadAlignments()
        line.read_alignments.clean()
        all = 0
        inter = 0
        contradicting = 0
        print "Unique segments:", segs
        if len(segs) == 0:
            print "WARNING: line with no resolved segments. Removing", line
            return
        for al in alignments:
            all += 1
            if segs.inter(al.seg_to, params.k):
                inter += 1
                if not al.contradictingRTC(tail_size=params.bad_end_length):
                    line.addReadAlignment(al)
                    # print "Added read alignment", al, str(al.seg_from.contig.alignments)
                else:
                    contradicting += 1
                    print "Contradicting read alignment", al, str(al.seg_from.contig.alignments)
            # else:
            #     print "Ambiguous read alignment", al, str(al.seg_from.contig.alignments)
        sys.stdout.trace("Read recruitment results:", all, inter, float(contradicting) / inter)
        line.updateCorrectSegments(line.asSegment())
        segs = segs.cap(line.correct_segments, params.k)
        line.completely_resolved.addAll(segs)
        sys.stdout.trace("The end")

    def markAllUnique(self, lines, reads):
        # type: (NewLineStorage, Iterable[AlignedRead]) -> None
        sys.stdout.info("Aligning reads to contigs")
        for al in self.aligner.localAlign(reads, lines):
            if len(al.seg_to) >= params.k:
                line = al.seg_to.contig # type: NewLine
                line.addReadAlignment(al)
        sys.stdout.info("Removing bad regions from lines")
        for line in list(lines.unique()):
            self.splitBad(line, lines)
        sys.stdout.info("Marking unique regions in lines")
        for line in lines.unique():
            self.markUniqueInLine(line)
        for line in list(lines.unique()):  # type:NewLine
            if len(line.completely_resolved) == 0:
                lines.remove(line)
            # else:
            #     line.initial.clean()
            #     for seg in line.completely_resolved:
            #         line.initial.add(AlignmentPiece.Identical(seg.asContig().asSegment(), seg))

    def medianCoverage(self, covs, line):
        clen = 0
        median_cov = 0
        for seg, cov in sorted(covs, key=lambda cov: cov[1]):
            clen += len(seg)
            if clen * 2 > len(line):
                median_cov = cov
                break
        return median_cov

    def splitBad(self, line, lines):
        # type: (NewLine, NewLineStorage) -> None
        segs = list(line.read_alignments.filterByCoverage(mi=params.reliable_coverage, k=params.k)) # type: List[Segment]
        assert len(segs) > 0, "No part of a unique edge is covered by reads"
        if len(segs) == 1 and len(segs[0]) > len(line) - 10:
            sys.stdout.warning("Whole line", line.id, "is covered by reads")
            lines.removeLine(line)
            return
        segs = filter(lambda seg: len(seg) >= params.k, segs)
        print "Line", line.id, "has poorly covered regions. Splitting into", len(segs), "parts"
        print segs
        next_left = segs[-1].left
        line.cutRight(segs[-1].right)
        for seg in list(segs)[-2:0:-1]:
            if next_left < seg.right:
                line, new_line = lines.splitLine(line.segment(next_left, seg.right))
            else:
                line, new_line = lines.splitLine(line.segment(next_left, next_left))
                line.cutRight(seg.right)
            next_left = seg.left
        line.rc.cutRight(len(segs[0]))


