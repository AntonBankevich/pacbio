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

    def markUniqueInLine(self, line, alter_cov):
        # type: (NewLine, float) -> None
        sys.stdout.info("Finding unique in", line)
        alignments = list(line.read_alignments) # type: List[AlignmentPiece]
        alignments = sorted(alignments, key=lambda al:al.seg_to.left)
        sys.stdout.trace("Sorting finished")
        inc = self.link(line, [al.seg_to.left for al in alignments if
                               al.seg_from.left > 1000 and al.seg_to.left > 50 and al.percentIdentity() > 0.85], 20)
        inc.append((line.segment(len(line) - 1, len(line)), alter_cov))
        alignments = sorted(alignments, key=lambda al:al.seg_to.right)
        out = self.link(line, [al.seg_to.right for al in alignments if al.rc.seg_from.left > 1000 and al.rc.seg_to.left > 50 ], 20)
        sys.stdout.trace("Linking finished")
        out.insert(0, (line.segment(0, 1), alter_cov))
        sys.stdout.trace( "inc:", inc)
        sys.stdout.trace( "out:", out)
        events = []
        for seg, val in inc:
            if val >= alter_cov:
                events.append((seg.left, -1))
        for seg, val in out:
            if val >= alter_cov:
                events.append((seg.right, 1))
        events= sorted(events)
        sys.stdout.trace("Events collected and sorted", len(events))
        events = [(pos, dir) for pos, dir in events if (dir == -1 or pos < len(line) - 200) and (dir == 1 or pos > - 200)]
        sys.stdout.trace( events)
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
        sys.stdout.trace("Unique segments selected")
        line.cleanReadAlignments()
        all = 0
        inter = 0
        contradicting = 0
        bad_quality = 0
        sys.stdout.trace( "Unique segments:", segs)
        if len(segs) == 0:
            sys.stdout.trace( "WARNING: line with no resolved segments. Removing", line)
            return
        for al in alignments:
            all += 1
            if segs.inter(al.seg_to, params.k):
                inter += 1
                if al.contradictingRTC(tail_size=params.bad_end_length):
                    contradicting += 1
                    sys.stdout.trace( "Contradicting read alignment", al, str(al.seg_from.contig.alignments))
                elif al.percentIdentity() < 0.85:
                    bad_quality += 1
                    sys.stdout.trace( "Read with bad alignment quality:", al)
                else:
                    sys.stdout.trace( "Recruited read:", al)
                    line.addReadAlignment(al)
            else:
                sys.stdout.trace("Alignment out of resolved:", al)

        sys.stdout.trace("Read recruitment results. All:", all, "In resolved regions:", inter,
                         "Contradicting:", float(contradicting) / inter, "Bad quality", float(bad_quality) / inter)
        line.updateCorrectSegments(line.asSegment())
        segs = segs.cap(line.correct_segments, params.k)
        line.completely_resolved.addAll(segs)
        sys.stdout.trace("The end")

    def markAllUnique(self, lines, reads):
        # type: (NewLineStorage, Iterable[AlignedRead]) -> None
        sys.stdout.info("Aligning reads to contigs")
        total_al_size = 0
        for al in self.aligner.localAlign(reads, lines):
            if len(al.seg_to) >= params.k:
                line = al.seg_to.contig # type: NewLine
                line.addReadAlignment(al)
                if (len(al.seg_to) > params.k + params.bad_end_length):
                    total_al_size += len(al.seg_to) - params.k - params.bad_end_length
        total_line_size = sum([len(line) for line in lines.unique()])
        alter_cov = max(float(total_al_size) / total_line_size / 2, params.min_k_mer_cov)
        sys.stdout.info("Average coverage for alternative k-mers calculated as ", alter_cov)
        sys.stdout.info("Removing bad regions from lines")
        self.splitBad(lines)
        sys.stdout.info("Marking unique regions in lines")
        for line in lines.unique():
            self.markUniqueInLine(line, alter_cov)
        for line in list(lines.unique()):  # type:NewLine
            if len(line.completely_resolved) == 0:
                lines.remove(line)
            # else:
            #     line.initial.clean()
            #     for seg in line.completely_resolved:
            #         line.initial.add(AlignmentPiece.Identical(seg.asContig().asSegment(), seg))

    def medianCoverage(self, covs):
        total = sum(len(seg) for seg, cov in covs)
        clen = 0
        median_cov = 0
        for seg, cov in sorted(covs, key=lambda cov: cov[1]):
            clen += len(seg)
            if clen * 2 >= total:
                return cov

    def splitBad(self, lines):
        # type: (NewLineStorage) -> None
        all_covs = []
        for line in lines:
            for rec in  line.read_alignments.calculateCoverage(params.k):
                all_covs.append(rec)
        median = self.medianCoverage(all_covs)
        sys.stdout.info("Median coverage determined as", median)
        lids = [line.id for line in lines.unique()]
        for line_id in lids:
            line = lines[line_id]
            s = AlignmentStorage()
            s.addAll(al for al in line.read_alignments if not al.contradictingRTC())
            segs = SegmentStorage().addAll(s.filterByCoverage(mi=params.reliable_coverage, ma=median * 7 /4, k=params.k))
            segs.mergeSegments(max(params.k - params.bad_end_length * 2, params.k / 2))
            if len(segs) == 0:
                sys.stdout.warn("No part of a unique edge is covered by reads", line.id)
                lines.removeLine(line)
                continue
            if len(segs) == 1 and len(segs[0]) > len(line) - 10:
                sys.stdout.info("Whole line", line.id, "is covered by reads")
                continue
            sys.stdout.info( "Line", line.id, "has poorly covered regions. Splitting into", len(segs), "parts")
            sys.stdout.trace(segs)
            next_left = segs[-1].left
            line.cutRight(segs[-1].right)
            for seg in list(segs)[-2::-1]:
                if next_left < seg.right:
                    line, new_line = lines.splitLine(line.segment(next_left, seg.right))
                else:
                    line, new_line = lines.splitLine(line.segment(next_left, next_left))
                    line.cutRight(seg.right)
                next_left = seg.left
            line.rc.cutRight(len(segs[0]))

