import itertools

from typing import List

from alignment.align_tools import Aligner
from common import basic
from common.alignment_storage import ReadCollection
from common.sequences import Segment, ContigStorage
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_storage import NewLineStorage


class CoverageAnalyser:
    def __init__(self, aligner, reads):
        # type: (Aligner, ReadCollection) -> None
        self.aligner = aligner
        self.reads = reads

    def analyseSegments(self, segs):
        # type: (List[Segment]) -> List[List[int]]
        contigs = ContigStorage()
        contigs.addAll([seg.asContig() for seg in segs if len(seg) > 5000])
        res = [] # type: List[Segment]
        for al in self.aligner.overlapAlign(self.reads, contigs):
            if basic.isCanonocal(al.seg_to.contig.id):
                res.append(al.seg_to)
            else:
                res.append(al.seg_to.RC())
        res = sorted(res, key=lambda seg: (seg.contig.id, seg.left))
        covs = [[0] * 20 for i in range(100)]
        for contig, it in itertools.groupby(res, key = lambda seg: seg.contig):
            segs = list(it)
            shrink = contig.asSegment().shrink(1000)
            for i in range(len(covs)):
                k = 500 + i * 100
                tmp = []
                for seg in segs:
                    if not seg.interSize(shrink) >= k:
                        seg = seg.cap(shrink)
                        tmp.append((seg.left, -1))
                        tmp.append((seg.right - k + 1, 1))
                tmp = sorted(tmp)
                cur_cov = 0
                prev = 1000
                for pos, diff in tmp:
                    assert pos >= prev
                    covs[i][min(cur_cov, len(covs[i]) - 1)] += pos - prev
                    cur_cov -= diff
                assert cur_cov == 0
                if prev < len(contig) - 1000 - k + 1:
                    covs[0] += len(contig) - 1000 - k + 1 - prev
        return covs

    def analyseLines(self, lines):
        # type: (NewLineStorage) -> List[List[int]]
        segs = []
        for line in lines.unique(): #type: NewLine
            segs.extend(line.completely_resolved)
        return self.analyseSegments(segs)


    def printAnalysis(self, covs):
        for i in range(len(covs)):
            print 500 + i * 100, ":", map(lambda cov: cov / sum(covs), covs[i])
