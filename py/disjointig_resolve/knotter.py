import itertools

from typing import Iterator, Tuple

from alignment.polishing import Polisher
from common.alignment_storage import AlignedRead, AlignmentPiece
from common.line_align import Scorer
from disjointig_resolve.accurate_line import NewLine, NewLineStorage


class LineKnotter:
    def __init__(self, storage, polisher):
        # type: (NewLineStorage, Polisher) -> None
        self.storage = storage
        self.polisher = polisher

    # Find connection of line to any other line using reads. Line is supposed to contain or precede the other line.
    def tryKnotRight(self, line):
        # type: (NewLine) -> bool
        read_alignments = line.read_alignments.allInter(line.asSegment().suffix(length=1000))
        candidates = []
        for al1 in read_alignments:
            read = al1.seg_from.contig # type: AlignedRead
            for al2 in read.alignments:
                if al1.canMergeTo(al2):
                    continue
                gap = len(read) - al1.seg_from.right - (len(line) - al1.seg_to.right) - (len(read) - al2.seg_from.left) - al2.seg_to.left
                if len(line) + gap > 0:
                    candidates.append((al2.seg_from.contig, gap, read, al1, al2))
        candidates = sorted(candidates)
        final_candidates = []
        for other_line, iter in itertools.groupby(lambda rec: rec[0], candidates): # type: NewLine, Iterator[Tuple[NewLine, int, AlignedRead, AlignmentPiece, AlignmentPiece]]
            recs = list(iter)
            if (recs[-1][1] - recs[0][1]) > min(100, abs(recs[-1][1]) / 8):
                print "\n".join(map(str, candidates))
                assert False, "Ambiguous knotting to the same line" + str(recs[0][0])
            avg = sum([rec[1] for rec in recs]) / len(recs)
            final_candidates.append((avg, avg + len(recs[0][0]), recs[0][0], [(rec[2], rec[3], rec[4]) for rec in recs]))
        if len(final_candidates) == 0:
            return False
        final_candidates = sorted(final_candidates)
        final = final_candidates[0]
        for candidate in final_candidates[1:]:
            if final[1] > candidate[0]:
                print "\n".join(map(str, candidates))
                assert False, "Contradicting candidates" + str(candidates[0][0]) + " " + str(candidates[0][1]) + str(candidate[0]) + " " + str(candidate[1])
        if final[0] > 0 or len(final[3]) <= 1:
            print "Positive gap. Can not merge."
            return False
        else:
            print "Merging", line, "with", final[2], "with gap", final[0]
            line_alignment = final[3][0][1].composeTargetDifference(final[3][0][2])
            pref = line_alignment.seg_from.left
            suff = len(line_alignment.seg_to.contig) - line_alignment.seg_to.right
            line_alignment = Scorer().polyshAlignment(line_alignment)
            new_line = self.storage.mergeLines(line_alignment, 500)
            seg = new_line.segment(pref, len(new_line) - suff)
            correction = self.polisher.polishSegment(seg, list(new_line.read_alignments.allInter(seg)))
            new_line.correctSequence([correction])
            new_line.updateCorrectSegments(new_line.segment(pref, len(new_line) - suff).expand(100))
            return True




