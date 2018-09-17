from typing import Optional

from alignment.align_tools import Aligner, AlignedSequences, ReadToAlignedSequences
from dag_resolve import params
from dag_resolve.line_tools import LineTail
from dag_resolve.sequences import AlignedRead, Contig, ReadCollection, ContigCollection

class Scorer:
    def __init__(self):
        pass

    def countHomo(self, seq, pos, step):
        cpos = pos + step
        while cpos >= 0 and cpos < len(seq):
            if seq[cpos] != seq[pos]:
                break
            cpos += step
        return (cpos - pos) * step

    def score(self, alignment, first, last):
        # type: (AlignedSequences, int, int) -> int
        assert alignment[first] != -1
        assert alignment[last] != -1
        matches = [(i, alignment[i]) for i in range(first, last + 1) if alignment[i] != 1]
        res = 0
        for match in matches:
            if not alignment.checkPosEquals(match[0]):
                res += params.Scores.sub_score
        for match1, match2 in zip(matches[:-1], matches[1:]):
            if not alignment.checkNoIndelAfter(match1[0]):
                continue
            l = [match2[0] - match1[0] - 1, match2[1] - match1[1] - 1]
            seqs = [alignment.seq_to[match1[0]: match2[0] + 1], alignment.seq_to[match1[1]: match2[1] + 1]]
            if l[1] > l[0]:
                homo = 0
                if alignment.checkPosEquals(match1[0]):
                    homo += self.countHomo(alignment.seq_from, match1[1], 1) - 1
                if alignment.checkPosEquals(match2[0]):
                    homo += self.countHomo(alignment.seq_from, match2[1], -1) - 1
                homo = min(max, (match2[1] - match1[1]) - (match2[0] - match1[0]))
                res += params.Scores.sub_score * (match2[0] - match1[0] - 1) + homo * params.Scores.homo_score + ((match2[1] - match1[1]) - (match2[0] - match1[0]) - homo) * params.Scores.ins_score


        return res


class LineAlignRecruter:
    def __init__(self, tails, consensus, aligner):
        # type: (list[LineTail], Contig, Aligner) -> LineAlignRecruter
        self.tails = tails
        self.tail_contigs = ContigCollection([Contig(tail.tail_consensus.seq, i) for i, tail in enumerate(tails)]) # type: ContigCollection
        self.consensus = consensus
        self.aligner = aligner

    def prepareReads(self, reads):
        # type: (ReadCollection) -> ReadCollection
        res = ReadCollection(self.tail_contigs)
        for i, tail in enumerate(self.tails):
            contig = self.tail_contigs[i]
            res.loadFromSam(self.aligner.align(reads, [contig]))
        return res


    def Recruit(self, read):
        # type: (AlignedRead) -> Optional[LineTail]
        alignment_list = []
        for i, tail in enumerate(self.tails):
            contig = self.tail_contigs[i]
            new_al = ReadToAlignedSequences(read, contig)
            alignment_list.append(new_al)
        return None