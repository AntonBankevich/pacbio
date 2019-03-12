from typing import Optional, Tuple, List

if __name__ == "__main__":
    import sys
    sys.path.append("py")
from common.alignment_storage import AlignmentPiece, MatchingSequence
from dag_resolve import params


class Scorer:
    def __init__(self, scores = None):
        if scores is None:
            scores = params.Scores()
        self.scores = scores

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
                res += self.scores.sub_score * l[0] + homo * self.scores.homo_score + (l[1] - l[0] - homo) * self.scores.del_score
            else:
                homo += self.countHomo(alignment.seq_from, match1[0], 1) - 1
                homo += self.countHomo(alignment.seq_from, match2[0], -1) - 1
                homo = min(homo, l[0] - l[1])
                res += self.scores.sub_score * l[1] + homo * self.scores.homo_score + (l[0] - l[1] - homo) * self.scores.del_score
        return res

    def accurateScore(self, alignment): #This score is nonsymmetric!!! Insertions and deletions have different cost
        # type: (MatchingSequence) -> int
        # print "Accurate scoring:", alignment[0], alignment[-1]
        cur = 0
        prev = Storage(alignment[0][1], alignment[1][1] + params.alignment_correction_radius, self.scores.inf)
        for j in range(alignment[0][1], alignment[1][1] + params.alignment_correction_radius + 1):
            prev.set(j, (j - alignment[0][1]) * self.scores.del_score)
        for i in range(alignment[0][0] + 1, alignment[-1][0] + 1):
            j_min = max(alignment[cur][1] - params.alignment_correction_radius, alignment[0][1])
            if alignment[cur + 1][0] == i and cur + 2 < len(alignment):
                cur += 1
                if alignment[cur + 1][1] - alignment[cur][1] > 10 and alignment[cur + 1][0] - alignment[cur][0] > 10:
                    print "Long gap:", alignment[cur], alignment[cur + 1]
            j_max = min(alignment[cur + 1][1] + params.alignment_correction_radius, alignment[-1][1])
            ne = Storage(j_min, j_max + 1, self.scores.inf)
            for j in range(j_min, j_max + 1):
                res = self.scores.inf
                if alignment.seq_from[i] == alignment.seq_to[j]:
                    res = min(res, prev.get(j - 1))
                    if i > 0 and j > 0 and alignment.seq_from[i - 1] == alignment.seq_from[i] and alignment.seq_to[j - 1] == alignment.seq_to[j]:
                        if i > 1 and alignment.seq_from[i - 2] == alignment.seq_from[i - 1]:
                            res = min(res, prev.get(j) + self.scores.homo_score)
                        if j > 1 and alignment.seq_to[j - 2] == alignment.seq_to[j - 1]:
                            res = min(res, ne.get(j - 1) + self.scores.homo_score)
                else:
                    res = min(res, prev.get(j - 1) + self.scores.sub_score)
                res = min(res, prev.get(j) + self.scores.ins_score)
                res = min(res, ne.get(j - 1) + self.scores.del_score)
                ne.set(j, res)
            prev = ne
        return prev.get(alignment[-1][1])

    def adjustMin(self, old_val, old_shift, new_val, new_shift):
        # type: (int, Tuple[int, int], int, Tuple[int, int]) -> Tuple[int, Tuple[int, int]]
        if new_val < old_val:
            return new_val, new_shift
        else:
            return old_val, old_shift

    def polyshAlignment(self, alignment):
        # type: (MatchingSequence) -> MatchingSequence
        # TEST
        cur = 0
        storage = RectStorage(alignment[0][0], alignment[-1][0])
        prev = Storage(alignment[0][1], alignment[1][1] + params.alignment_correction_radius, self.scores.inf)
        storage.set(alignment[0][0], Storage(alignment[0][1], alignment[1][1] + params.alignment_correction_radius))
        for j in range(alignment[0][1], alignment[1][1] + params.alignment_correction_radius + 1):
            prev.set(j, (j - alignment[0][1]) * self.scores.del_score)
        for i in range(alignment[0][0] + 1, alignment[-1][0] + 1):
            j_min = max(alignment[cur][1] - params.alignment_correction_radius, alignment[0][1])
            if alignment[cur + 1][0] == i and cur + 2 < len(alignment):
                cur += 1
                if alignment[cur + 1][1] - alignment[cur][1] > 10 and alignment[cur + 1][0] - alignment[cur][0] > 10:
                    print "Long gap:", alignment[cur], alignment[cur + 1]
            j_max = min(alignment[cur + 1][1] + params.alignment_correction_radius, alignment[-1][1])
            ne = Storage(j_min, j_max + 1, self.scores.inf)
            storage.set(i, Storage(j_min, j_max + 1, None))
            for j in range(j_min, j_max + 1):
                res = self.scores.inf
                res_shift = (0, 0)
                if alignment.seq_from[i] == alignment.seq_to[j]:
                    res, res_shift = self.adjustMin(res, res_shift, prev.get(j - 1), (-1, -1))
                    if i > 0 and j > 0 and alignment.seq_from[i - 1] == alignment.seq_from[i] and alignment.seq_to[j - 1] == alignment.seq_to[j]:
                        if i > 1 and alignment.seq_from[i - 2] == alignment.seq_from[i - 1]:
                            res, res_shift = self.adjustMin(res, res_shift, prev.get(j) + self.scores.homo_score, (-1, 0))
                        if j > 1 and alignment.seq_to[j - 2] == alignment.seq_to[j - 1]:
                            res, res_shift = self.adjustMin(res, res_shift, ne.get(j - 1) + self.scores.homo_score, (0, -1))
                else:
                    res, res_shift = self.adjustMin(res, res_shift, prev.get(j - 1) + self.scores.sub_score, (-1, -1))
                res, res_shift = self.adjustMin(res, res_shift, prev.get(j) + self.scores.ins_score, (-1, 0))
                res, res_shift = self.adjustMin(res, res_shift, ne.get(j - 1) + self.scores.del_score, (0, -1))
                ne.set(j, res)
                storage.get(i).set((i + res_shift[0], j + res_shift[1]))
            prev = ne
        cur_i = alignment[-1][0]
        cur_j = alignment[-1][1]
        matches = [(cur_i, cur_j)]
        while cur_i != alignment[0][0] or cur_j != alignment[0][1]:
            cur_i, cur_j = storage.get(cur_i).get(cur_j)
            matches.append((cur_i, cur_j))
        return MatchingSequence(alignment.seq_from, alignment.seq_to, matches[::-1])

    def score3(self, piece1, piece2):
        # type: (AlignmentPiece, AlignmentPiece) -> Tuple[Optional[int], Optional[int], Optional[int]]
        pid1 = piece1.percentIdentity()
        pid2 = piece2.percentIdentity()
        # print "Percent identities:", pid1, pid2
        if pid1 < params.min_expected_Pacbio_PI and pid2 < params.min_expected_Pacbio_PI:
            return None, None, None
        # we only consider contradictions with the right part of the line since reads were aligned to line center suffixes
        contra1 = piece1.contradicting(piece1.seg_to.contig.centerPos.suffix())
        contra2 = piece2.contradicting(piece2.seg_to.contig.centerPos.suffix())
        if pid1 < params.min_allowed_Pacbio_PI or (contra1 and piece2.seg_from.right > piece1.seg_from.right + 500):
            return None, self.score(piece2.matchingSequence()), None
        if pid2 < params.min_allowed_Pacbio_PI or (contra2 and piece1.seg_from.right > piece2.seg_from.right + 500):
            return self.score(piece1.matchingSequence()), None, None
        matches1 = piece1.matchingSequence()
        matches2 = piece2.matchingSequence()
        composite = matches1.composeDifference(matches2)
        matches1 = matches1.reduceTarget(composite.matches[0][0], composite.matches[-1][0] + 1)
        matches2 = matches2.reduceTarget(composite.matches[0][1], composite.matches[-1][1] + 1)
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

class RectStorage(Storage):
    def __init__(self, left, right, default = None):
        Storage.__init__(left, right, default)
        self.vals = self.vals # type: List[Storage]

if __name__ == "__main__":
    tmp = MatchingSequence(sys.argv[1], sys.argv[2], [(0, 0), (len(sys.argv[1]) - 1, len(sys.argv[2]) - 1)])
    print Scorer().accurateScore(tmp)

