from typing import Optional, Tuple

from dag_resolve.sequences import AlignmentPiece, MatchingSequence
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
        # type: (MatchingSequence) -> int
        matches = alignment.matches
        res = 0
        for match1, match2 in zip(matches[:-1], matches[1:]):
            if match2[0] - match1[0] == 1 and match2[1] - match1[1] == 1:
                continue
            l = [match2[0] - match1[0] - 1, match2[1] - match1[1] - 1]
            homo = 0
            if l[1] > l[0]:
                homo += self.countHomo(alignment.seq_to, match1[1], 1) - 1
                homo += self.countHomo(alignment.seq_to, match2[1], -1) - 1
                homo = min(homo, (match2[1] - match1[1]) - (match2[0] - match1[0]))
                res += params.Scores.sub_score * l[0] + homo * params.Scores.homo_score + (l[1] - l[0] - homo) * params.Scores.del_score
            else:
                homo += self.countHomo(alignment.seq_from, match1[0], 1) - 1
                homo += self.countHomo(alignment.seq_from, match2[0], -1) - 1
                homo = min(homo, l[0] - l[1])
                res += params.Scores.sub_score * l[1] + homo * params.Scores.homo_score + (l[0] - l[1] - homo) * params.Scores.del_score
        return res

    def accurateScore(self, alignment): #This score is nonsymmetric!!! Insertions and deletions have different cost
        # type: (MatchingSequence) -> int
        # print "Accurate scoring:", alignment[0], alignment[-1]
        cur = 0
        prev = Storage(alignment[0][1], alignment[1][1] + params.alignment_correction_radius, params.Scores.inf)
        for j in range(alignment[0][1], alignment[1][1] + params.alignment_correction_radius + 1):
            prev.set(j, (j - alignment[0][1]) * params.Scores.del_score)
        for i in range(alignment[0][0] + 1, alignment[-1][0] + 1):
            j_min = max(alignment[cur][1] - params.alignment_correction_radius, alignment[0][1])
            if alignment[cur + 1][0] == i and cur + 2 < len(alignment):
                cur += 1
                if alignment[cur + 1][1] - alignment[cur][1] > 10 and alignment[cur + 1][0] - alignment[cur][0] > 10:
                    print "Long gap:", alignment[cur], alignment[cur + 1]
            j_max = min(alignment[cur + 1][1] + params.alignment_correction_radius, alignment[-1][1])
            ne = Storage(j_min, j_max + 1, params.Scores.inf)
            # assert j_max - j_min < 100
            for j in range(j_min, j_max + 1):
                res = params.Scores.inf
                if alignment.seq_from[i] == alignment.seq_to[j]:
                    res = min(res, prev.get(j - 1))
                    if i > 0 and j > 0 and alignment.seq_from[i - 1] == alignment.seq_from[i] and alignment.seq_to[j - 1] == alignment.seq_to[j]:
                        if i > 1 and alignment.seq_from[i - 2] == alignment.seq_from[i - 1]:
                            res = min(res, prev.get(j) + params.Scores.homo_score)
                        if j > 1 and alignment.seq_to[j - 2] == alignment.seq_to[j - 1]:
                            res = min(res, ne.get(j - 1) + params.Scores.homo_score)
                else:
                    res = min(res, prev.get(j - 1) + params.Scores.sub_score)
                res = min(res, prev.get(j) + params.Scores.ins_score)
                res = min(res, ne.get(j - 1) + params.Scores.del_score)
                ne.set(j, res)
            prev = ne
            # if alignment[cur][0] == i and cur + 1 < len(alignment):
            #     print alignment[cur], ne.get(alignment[cur][1])
            #     if cur > 0:
            #         print alignment.seq_from[alignment[cur - 1][0] : alignment[cur][0]], alignment.seq_to[alignment[cur - 1][1]: alignment[cur][1]]
        return prev.get(alignment[-1][1])


    def score3(self, piece1, piece2, alignments):
        # type: (AlignmentPiece, AlignmentPiece, list[AlignmentPiece]) -> Tuple[Optional[int], Optional[int], Optional[int]]
        pid1 = piece1.percentIdentity()
        pid2 = piece2.percentIdentity()
        # print "Percent identities:", pid1, pid2
        if pid1 < params.min_expected_Pacbio_PI and pid2 < params.min_expected_Pacbio_PI:
            return None, None, None
        if pid1 < params.min_expected_Pacbio_PI:
            return None, self.score(piece2.matchingSequence()), None
        if pid2 < params.min_expected_Pacbio_PI:
            return self.score(piece1.matchingSequence()), None, None
        matches1 = piece1.matchingSequence()
        matches2 = piece2.matchingSequence()
        # print "Matches:"
        # print matches1.matches
        # print matches2.matches
        composite = matches1.compose(matches2)
        # print "Composite:", composite.matches
        matches1 = matches1.reduceReference(composite.matches[0][0], composite.matches[-1][0] + 1)
        matches2 = matches2.reduceReference(composite.matches[0][1], composite.matches[-1][1] + 1)
        # print "New matches:"
        # print matches1.matches
        # print matches2.matches
        # print "Line alignment:"
        # print map(str, alignments)
        tmp_composite = composite.combine(map(lambda al: al.matchingSequence(), alignments)) # refactor
        # print "Old scores:", self.score(matches1), self.score(matches2), self.score(tmp_composite)
        accurate1 = self.accurateScore(matches1)
        accurate2 = self.accurateScore(matches2)
        if accurate1 > accurate2:
            composite = composite.reverse()
        accurate12 = self.accurateScore(composite)
        if not (abs(accurate1 - accurate2) <= accurate12 <= accurate1 + accurate2):
            print "Triangle inequality failed: " + \
                  str(accurate1) + " " + str(accurate2) + " " + \
                  str(abs(accurate1 - accurate2)) + "<=" + str(accurate12) + "<=" +str(accurate1 + accurate2)
        return self.accurateScore(matches1), self.accurateScore(matches2), self.accurateScore(composite)

class Storage:
    def __init__(self, left, right, default = None):
        self.left = left
        self.right = right
        self.default = default
        self.vals = []
        for i in range(left, right + 1):
            self.vals.append(default)

    def set(self, a, val):
        assert self.left <= a <= self.right
        self.vals[a - self.left] = val

    def get(self, a):
        if self.left <= a <= self.right:
            return self.vals[a - self.left]
        else:
            return self.default

class RectStorage:
    def __init__(self, left, right, default = None):
        self.left = left
        self.right = right
        self.default = default
        self.vals = []
        for i in range(left, right + 1):
            self.vals.append(dict())

    def set(self, a, b, val):
        assert self.left <= self.right <= b
        self.vals[a - self.left][b] = val

    def get(self, a, b):
        assert self.left <= self.right <= b
        if b in self.vals[a - self.left]:
            return self.vals[a - self.left][b]
        return self.default



