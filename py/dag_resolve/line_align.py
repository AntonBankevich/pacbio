from typing import Optional, Tuple

from dag_resolve.sequences import AlignmentPiece
from dag_resolve import params


class Scorer:
    def __init__(self):
        pass

    def countHomo(self, seq, pos, step):
        cpos = pos + step
        while 0 <= cpos < len(seq):
            if seq[cpos] != seq[pos]:
                break
            cpos += step
        return (cpos - pos) * step

    def score(self, alignment):
        # type: (MatchingSequence, int, int) -> int
        matches = alignment.matches
        res = 0
        for match1, match2 in zip(matches[:-1], matches[1:]):
            if match2[0] - match1[0] == 1 and match2[1] - match1[1] == 1:
                continue
            l = [match2[0] - match1[0] - 1, match2[1] - match1[1] - 1]
            if l[1] > l[0]:
                homo = 0
                homo += self.countHomo(alignment.seq_from, match1[1], 1) - 1
                homo += self.countHomo(alignment.seq_from, match2[1], -1) - 1
                homo = min(homo, (match2[1] - match1[1]) - (match2[0] - match1[0]))
                res += params.Scores.sub_score * (match2[0] - match1[0] - 1) + homo * params.Scores.homo_score + ((match2[1] - match1[1]) - (match2[0] - match1[0]) - homo) * params.Scores.ins_score
        return res

    def score3(self, piece1, piece2, alignments):
        # type: (AlignmentPiece, AlignmentPiece, list[AlignmentPiece]) -> Tuple[Optional[int], Optional[int], Optional[int]]
        pid1 = piece1.percentIdentity()
        pid2 = piece2.percentIdentity()
        if pid1 < 0.8 and pid2 < 0.8:
            return None, None, None
        if pid1 < 0.8:
            return None, self.score(piece2.MatchingSequence()), None
        if pid2 < 0.8:
            return self.score(piece1.MatchingSequence()), None, None
        matches1 = piece1.MatchingSequence()
        matches2 = piece2.MatchingSequence()
        composite = matches1.compose(matches2)
        matches1 = matches1.reduce(composite.matches[0][0], composite.matches[-1][0])
        matches2 = matches2.reduce(composite.matches[0][1], composite.matches[-1][1])
        composite = composite.combine(map(lambda al: al.MatchingSequence(), alignments))
        return self.score(matches1), self.score(matches2), self.score(composite)



