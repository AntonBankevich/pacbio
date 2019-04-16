from typing import Generator, Iterable, List, Tuple

from common import params
from common.sequences import Segment
from disjointig_resolve.accurate_line import NewLine, CoverageCalculator
from disjointig_resolve.disjointigs import DisjointigCollection
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
        alignments = filter(lambda al: len(al.seg_to) > params.k, line.getPotentialAlignmentsTo(line.asSegment()))
        inc = self.link(line, [al.seg_to.left for al in alignments if al.seg_from.left > 1000], 20)
        out = self.link(line, [al.seg_to.right for al in alignments if len(al.seg_from.contig) - al.seg_from.right > 1000], 20)
        inc = SegmentStorage().addAll([seg for seg, cov in inc if cov >= 5]).reverse()
        out = SegmentStorage().addAll([seg for seg, cov in out if cov >= 5]).reverse()
        segs1 = inc.orderedCap(out)
        line_als = AlignmentStorage().addAll(dot_plot.allInter(line.asSegment()))
        segs2 = line_als.filterByCoverage(line_als, 1, 2)
        segs = segs1.cap(segs2).expand(params.k / 2).filterBySize(min = params.k)
        for al in alignments:
            if segs.inter(al.seg_to):
                line.addReadAlignment(al)
        line.updateCorrectSegments(line.asSegment())
        segs = segs.cap(line.correct_segments, params.k)
        line.completely_resolved.addAll(segs)

    def medianCoverage(self, covs, line):
        clen = 0
        median_cov = 0
        for seg, cov in sorted(covs, key=lambda cov: cov[1]):
            clen += len(seg)
            if clen * 2 > len(line):
                median_cov = cov
                break
        return median_cov

