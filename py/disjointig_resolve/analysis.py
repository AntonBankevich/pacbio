import itertools

from typing import List

from alignment.align_tools import Aligner
from common import basic, params
from common.alignment_storage import ReadCollection
from common.sequences import Segment, ContigStorage
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_storage import NewLineStorage


class CoverageAnalyser:
    def __init__(self, aligner, reads):
        # type: (Aligner, ReadCollection) -> None
        self.aligner = aligner
        self.reads = reads
        self.recs = None #type: List[CoverageAnalyser.CoverageRecord]

    class CoverageRecord:
        def __init__(self, k, covs):
            self.k = k
            self.covs = covs
            for i in range(0, len(covs) - 1):
                self.covs[i + 1] += self.covs[i]

        def __str__(self):
            if len(self.covs) > 0 and self.covs[-1] > 0:
                return "k:" + str(self.k) + ": " + " ".join(map(lambda cov: "%0.3f" % (float(cov) / self.covs[-1]), self.covs))
            else:
                return "None"

    def analyseSegments(self, segs):
        # type: (List[Segment]) -> None
        contigs = ContigStorage()
        contigs.addAll([seg.asContig() for seg in segs if len(seg) > 5000])
        res = [] # type: List[Segment]
        for al in self.aligner.overlapAlign(self.reads, contigs):
            if basic.isCanonocal(al.seg_to.contig.id):
                res.append(al.seg_to)
            else:
                res.append(al.seg_to.RC())
        res = sorted(res, key=lambda seg: (seg.contig.id, seg.left))
        covs = [[0] * params.maxCoverageThreshold for i in range(100)]
        for contig, it in itertools.groupby(res, key = lambda seg: seg.contig):
            segs = list(it)
            shrink = contig.asSegment().shrink(1000)
            for i in range(len(covs)):
                k = 500 + i * 100
                tmp = []
                for seg in segs:
                    if seg.interSize(shrink) >= k:
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
                    prev = pos
                assert cur_cov == 0
                if prev < len(contig) - 1000 - k + 1:
                    covs[i][0] += len(contig) - 1000 - k + 1 - prev
        self.recs = [CoverageAnalyser.CoverageRecord(500 + i * 100, covs[i]) for i in range(len(covs))]

    def analyseLines(self, lines):
        # type: (NewLineStorage) -> None
        segs = []
        for line in lines.unique(): #type: NewLine
            segs.extend(line.completely_resolved)
        self.analyseSegments(segs)


    def printAnalysis(self):
        for cov in self.recs:
            print cov

    def chooseK(self, threshold):
        assert threshold < params.maxCoverageThreshold
        res = None
        for rec in self.recs:
            if rec.covs[threshold] < params.uncoveredFractionForK:
                res = rec.k
            else:
                break
        return res

