import itertools

from typing import Iterator, Tuple, Optional, List

from alignment.polishing import Polisher
from common.alignment_storage import AlignedRead, AlignmentPiece
from common.line_align import Scorer
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_storage import NewLineStorage


class LineKnotter:
    def __init__(self, storage, polisher):
        # type: (NewLineStorage, Polisher) -> None
        self.storage = storage
        self.polisher = polisher

    class Record:
        def __init__(self, al1, al2):
            # type: (AlignmentPiece, AlignmentPiece) -> None
            self.read = al1.seg_from.contig # type: AlignedRead
            line = al1.seg_to.contig # type: NewLine
            self.gap = len(self.read) - al1.seg_from.right - (len(line) - al1.seg_to.right) - (len(self.read) - al2.seg_from.left) - al2.seg_to.left
            self.al1 = al1
            self.al2 = al2
            self.other = al2.seg_to.contig # type: NewLine

        def __str__(self):
            return str([self.other, self.gap, self.al1, self.al2])

    # Find connection of line to any other line using reads. Line is supposed to contain or precede the other line.
    def tryKnotRight(self, line):
        # type: (NewLine) -> Optional[NewLine]
        read_alignments = line.read_alignments.allInter(line.asSegment().suffix(length=1000))
        candidates = [] # type: List[LineKnotter.Record]
        for al1 in read_alignments:
            read = al1.seg_from.contig # type: AlignedRead
            for al2 in read.alignments:
                if al1.canMergeTo(al2):
                    continue
                new_rec = self.Record(al1, al2)
                if len(line) + new_rec.gap > 0:
                    candidates.append(new_rec)
        #             (al2.seg_from.contig, gap, read, al1, al2)
        candidates = sorted(candidates, key = lambda rec: rec.other.id)
        final_candidates = []
        for other_line, iter in itertools.groupby(candidates, lambda rec: rec.other): # type: NewLine, Iterator[LineKnotter.Record]
            recs = list(iter)
            if (recs[-1].gap - recs[0].gap) > min(100, abs(recs[-1].gap) / 8):
                print "\n".join(map(str, candidates))
                assert False, "Ambiguous knotting to the same line" + str(recs[0])
            avg = sum([rec.gap for rec in recs]) / len(recs)
            final_candidates.append((avg, recs[0].other, recs))
        if len(final_candidates) == 0:
            return None
        final_candidates = sorted(final_candidates)
        final = final_candidates[0]
        for candidate in final_candidates[1:]:
            if final[0] + len(final[1]) > candidate[0]:
                print "\n".join(map(str, candidates))
                assert False, "Contradicting candidates" + str(final[0]) + " " + str(final[1]) + " " + str(candidate[0]) + " " + str(candidate[1])
        if final[0] > 0 or len(final[2]) <= 1:
            print "Positive gap. Can not merge."
            return None
        else:
            print "Merging", line, "with", final[1], "with gap", final[0]
            line_alignment = final[2][0].al1.composeTargetDifference(final[2][0].al2)
            pref = line_alignment.seg_from.left
            suff = len(line_alignment.seg_to.contig) - line_alignment.seg_to.right
            line_alignment = Scorer().polyshAlignment(line_alignment)
            new_line = self.storage.mergeLines(line_alignment, 500)
            seg = new_line.segment(pref, len(new_line) - suff)
            correction = self.polisher.polishSegment(seg, list(new_line.read_alignments.allInter(seg)))
            new_line.correctSequence([correction])
            new_line.updateCorrectSegments(new_line.segment(pref, len(new_line) - suff).expand(100))
            return new_line




