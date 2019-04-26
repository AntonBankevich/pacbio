from typing import Generator, Iterable, List, Tuple

from common import params
from common.sequences import Segment
from disjointig_resolve.accurate_line import NewLine, NewLineStorage
from disjointig_resolve.disjointigs import DisjointigCollection, Disjointig
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage


class UniqueMarker:
    def __init__(self):
        # type: () -> None
        pass

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

    def markUniqueInLine(self, line, dot_plot):
        # type: (NewLine, LineDotPlot) -> None
        alignments = filter(lambda al: len(al.seg_to) > params.k, line.getRelevantAlignmentsFor(line.asSegment()))
        inc = self.link(line, [al.seg_to.left for al in alignments if al.seg_from.left > 1000 and al.seg_to.left > 50], 20)
        inc.append((line.segment(len(line) - 1, len(line)), 5))
        out = self.link(line, [al.seg_to.right for al in alignments if al.rc.seg_from.left > 1000 and al.rc.seg_to.left > 50 ], 20)
        out.insert(0, (line.segment(0, 1), 5))
        # print inc
        # print [al for al in alignments if al.rc.seg_from.left < 1000 and al.seg_to.left > 50]
        # print out
        inc = SegmentStorage().addAll([seg for seg, cov in inc if cov >= 5]).reverse()
        out = SegmentStorage().addAll([seg for seg, cov in out if cov >= 5]).reverse()
        segs1 = inc.orderedCap(out)
        # print segs1
        line_als = AlignmentStorage().addAll(dot_plot.allInter(line.asSegment()))
        segs2 = line_als.filterByCoverage(1, 2)
        # print segs2
        segs = segs1.cap(segs2).expand(params.k / 2).filterBySize(min = params.k)
        # print segs
        for al in alignments:
            if segs.inter(al.seg_to, params.k):
                line.addReadAlignment(al)
        line.updateCorrectSegments(line.asSegment())
        segs = segs.cap(line.correct_segments, params.k)
        # print segs
        line.completely_resolved.addAll(segs)

    def markAllUnique(self, lines, dot_plot):
        # type: (NewLineStorage, LineDotPlot) -> None
        for line in lines.unique():
            self.markUniqueInLine(line, dot_plot)

    def medianCoverage(self, covs, line):
        clen = 0
        median_cov = 0
        for seg, cov in sorted(covs, key=lambda cov: cov[1]):
            clen += len(seg)
            if clen * 2 > len(line):
                median_cov = cov
                break
        return median_cov

