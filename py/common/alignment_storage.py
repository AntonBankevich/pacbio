import itertools

from typing import Generator, Tuple, Optional, Any, List, Dict, Callable, Iterator, Iterable, BinaryIO

import common.seq_records
from common import sam_parser, params, easy_cigar, basic, SeqIO
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, ContigCollection, Contig, UniqueList


# TODO merge with MatchingSequence. Make cigar a GlobalAlignment class with effective operations.
class AlignmentPiece:
    def __init__(self, seg_from, seg_to, cigar, rc = None):
        # type: (Segment, Segment, str, Optional[AlignmentPiece]) -> None
        self.seg_from = seg_from #type: Segment
        self.seg_to = seg_to #type: Segment
        assert cigar != "X"
        assert cigar.find("H") == -1 and cigar.find("S") == -1
        if cigar == "=":
            cigar = str(len(seg_from)) + "M"
        self.cigar = cigar
        assert len(seg_from) > 0 and len(seg_to) > 0
        assert seg_from.contig[seg_from.left] == seg_to.contig[seg_to.left]
        assert seg_from.contig[seg_from.right - 1] == seg_to.contig[seg_to.right - 1]
        if params.assert_pi:
            pi = self.matchingPercentIdentity()
            if pi < 0.05:
                print seg_from.Seq()
                print seg_to.Seq()
            assert pi >= 0.5, str(self)
        if rc is None:
            self.rc = AlignmentPiece(seg_from.RC(), seg_to.RC(), easy_cigar.RCCigar(self.cigar), self)
        else:
            self.rc = rc

    @staticmethod
    def FromSamRecord(seq_from, seq_to, rec):
        # type: (Contig, Contig, sam_parser.SAMEntryInfo) -> AlignmentPiece
        cigar_list = list(easy_cigar.CigarToList(rec.cigar))
        ls = 0
        rs = 0
        if cigar_list[0][1] in "HS":
            ls = cigar_list[0][0]
        if cigar_list[-1][1] in "HS":
            rs = cigar_list[-1][0]
        new_cigar = []
        for num, s in cigar_list:
            if s not in "HS":
                new_cigar.append(num)
                new_cigar.append(s)
        new_cigar = "".join(map(str, new_cigar))
        seg_from = Segment(seq_from, ls, len(seq_from) - rs)
        seg_to = Segment(seq_to, rec.pos - 1, rec.pos - 1 + rec.alen)
        if rec.rc:
            seg_from = Segment(seq_from.rc, ls, len(seq_from) - rs).RC()
            seg_to = seg_to.RC()
            new_cigar = easy_cigar.RCCigar(new_cigar)
        piece = AlignmentPiece(seg_from, seg_to, new_cigar)
        # seg_from.contig.alignments.append(piece)
        # seg_from.contig.rc.alignments.append(piece.rc)
        return piece

    @staticmethod
    def Identical(seg_from, other = None):
        # type: (Segment, Optional[Segment]) -> AlignmentPiece
        if other is None:
            other = seg_from
        return AlignmentPiece(seg_from, other, str(len(seg_from)) + "M")

    @staticmethod
    def GlueOverlappingAlignments(als):# Glue a list of alignments such that consecutive target segments overlap. Can return None
        # type: (List[AlignmentPiece]) -> AlignmentPiece
        contig = als[0].seg_to.contig
        als1 = [al.matchingSequence() for al in als]
        truncated = [als[0].seg_to]
        last = als[0].matchingSequence()
        for al in als[1:]:
            next = al.matchingSequence()
            shared = list(last.common_to(next))
            if len(shared) == 0:
                return None
            pos1, pos2 = shared[len(shared) / 2]
            truncated[-1] = truncated[-1].prefix(pos=last.matches[pos1][1])
            truncated.append(al.seg_to.suffix(pos=next.matches[pos2][1]))
            last = next
        als = [al.reduce(target=seg) for al, seg in zip(als, truncated)]
        return AlignmentPiece.GlueFittingAlignments(als)

    @staticmethod
    def MergeOverlappingAlignments(als):# Glue a list of alignments such that consecutive query and target segments overlap. Can return None
        # type: (List[AlignmentPiece]) -> AlignmentPiece
        als1 = [al.matchingSequence() for al in als]
        truncated = [als[0].seg_to]
        last = als[0].matchingSequence()
        for al in als[1:]:
            next = al.matchingSequence()
            shared = list(last.common_to(next))
            if len(shared) == 0:
                return None
            pos1, pos2 = shared[len(shared) / 2]
            truncated[-1] = truncated[-1].prefix(pos=last.matches[pos1][1])
            truncated.append(al.seg_to.suffix(pos=next.matches[pos2][1]))
            last = next
        als = [al.reduce(target=seg) for al, seg in zip(als, truncated)]
        return AlignmentPiece.MergeFittingAlignments(als)

    @staticmethod
    def GlueFittingAlignments(als): # Glue a list of alignments with seg_to exactly fitting together.
        # type: (List[AlignmentPiece]) -> AlignmentPiece
        contig = als[0].seg_to.contig
        new_seq = "".join((al.seg_from.Seq() for al in als))
        new_contig = Contig(new_seq, "glued")
        new_cigar = "".join(al.cigar for al in als)
        return AlignmentPiece(new_contig.asSegment(), contig.segment(als[0].seg_to.left, als[-1].seg_to.right),
                              new_cigar)

    @staticmethod
    def MergeFittingAlignments(als): # Glue a list of alignments with seg_from and seg_to exactly fitting together
        # type: (List[AlignmentPiece]) -> AlignmentPiece
        contig_from = als[0].seg_from.contig # type: Contig
        contig_to = als[0].seg_to.contig
        new_cigar = "".join(al.cigar for al in als)
        return AlignmentPiece(contig_from.segment(als[0].seg_from.left, als[-1].seg_from.right),
                              contig_to.segment(als[0].seg_to.left, als[-1].seg_to.right),
                              new_cigar)

    def isIdentical(self):
        return self.seg_from == self.seg_to

    def __str__(self):
        # type: () -> str
        if len(self) < 20000:
            pid = self.percentIdentity()
        else:
            pid = 1
        if pid > 0.99:
            spid = "%0.3f" % pid
        else:
            spid = "%0.2f" % pid
        suffix = ""
        if self.contradicting():
            suffix = "!!!"
        return "(" + str(self.seg_from) + "->" + str(self.seg_to) + ":" + spid + suffix + ")"

    def __repr__(self):
        return self.__str__()

    def changeQueryContig(self, read):
        return AlignmentPiece(self.seg_from.changeContig(read), self.seg_to, self.cigar)

    def changeTargetContig(self, contig):
        return AlignmentPiece(self.seg_from, self.seg_to.changeContig(contig), self.cigar)

    def changeQuerySegment(self, seg):
        # type: (Segment) -> AlignmentPiece
        return AlignmentPiece(seg, self.seg_to, self.cigar)

    def changeTargetSegment(self, seg):
        # type: (Segment) -> AlignmentPiece
        return AlignmentPiece(self.seg_from, seg, self.cigar)

    def precedes(self, other, delta=0):
        # type: (AlignmentPiece, int) -> bool
        return self.seg_from.precedes(other.seg_from, delta) and self.seg_to.precedes(other.seg_to, delta)

    def __le__(self, other):
        # type: (AlignmentPiece) -> bool
        return self.seg_from <= other.seg_from and self.seg_to <= other.seg_to

    def canMergeTo(self, other):
        # type: (AlignmentPiece) -> bool
        return self <= other or other <= self or self.contains(other) or other.contains(self)

    def connects(self, other, delta=0):
        # type: (AlignmentPiece, int) -> bool
        return self.seg_from.connects(other.seg_from, delta) and self.seg_to.connects(other.seg_to, delta)

    def contains(self, other):
        # type: (AlignmentPiece) -> bool
        return self.seg_from.contains(other.seg_from) and self.seg_to.contains(other.seg_to)

    # TODO: do not use this
    def mergeDistant(self, other):
        ins = ""
        if not self.connects(other):
            d = (other.seg_from.left - self.seg_from.right, other.seg_to.left - self.seg_to.right)
            assert d[0] < 10 and d[1] < 10
            if min(d[0], d[1]) > 0:
                ins += str(min(d[0], d[1])) + "M"
            if d[0] > d[1]:
                ins += str(d[0] - d[1]) + "I"
            elif d[1] > d[0]:
                ins += str(d[1] - d[0]) + "D"
        return AlignmentPiece(self.seg_from.merge(other.seg_from), self.seg_to.merge(other.seg_to),
                              self.cigar + ins + other.cigar)

    def __eq__(self, other):
        return self.seg_from == other.seg_from and self.seg_to == other.seg_to and self.cigar == other.cigar

    def __ne__(self, other):
        return self.seg_from != other.seg_from or self.seg_to != other.seg_to or self.cigar != other.cigar

    def __len__(self):
        return len(self.seg_from)

    # def RC(self):
    #     return AlignmentPiece(self.seg_from.RC(), self.seg_to.RC(), sam_parser.RCCigar(self.cigar))

    def matchingPositions(self, equalOnly=False):
        # type: (bool) -> Generator[Tuple[int, int]]
        assert self.cigar != "="
        cur_query = self.seg_from.left
        cur_tar = self.seg_to.left
        for n, c in easy_cigar.CigarToList(self.cigar):
            if c == 'M':
                for i in range(n):
                    if not equalOnly or self.seg_from.contig[cur_query] == self.seg_to.contig[cur_tar]:
                        yield (cur_query, cur_tar)
                    cur_tar += 1
                    cur_query += 1
            elif c == "D":
                cur_tar += n
            elif c == "I":
                cur_query += n

    def asMatchingStrings(self):
        pos_pairs = list(self.matchingPositions())
        l1 = []
        l2 = []
        for p1, p2 in zip(pos_pairs[:-1], pos_pairs[1:]):
            l1.append(self.seg_from.contig[p1[0]])
            l2.append(self.seg_to.contig[p1[1]])
            lens = [p2[0] - p1[0] - 1, p2[1] - p1[1] - 1]
            for i in range(min(lens[0], lens[1])):
                l1.append(self.seg_from.contig[p1[0] + i + 1])
                l2.append(self.seg_to.contig[p1[1] + i + 1])
            for i in range(min(lens[0], lens[1]), lens[0]):
                l1.append(self.seg_from.contig[p1[0] + i + 1])
                l2.append("-")
            for i in range(min(lens[0], lens[1]), lens[1]):
                l1.append("-")
                l2.append(self.seg_to.contig[p1[1] + i + 1])
        l1.append(self.seg_from.contig[pos_pairs[-1][0]])
        l2.append(self.seg_to.contig[pos_pairs[-1][1]])
        return "".join(l1), "".join(l2)

    def percentIdentity(self):
        res = len(list(self.matchingPositions(True)))
        return float(res) / max(len(self.seg_from), len(self.seg_to))

    def matchingPercentIdentity(self):
        res = len(list(self.matchingPositions(True)))
        all = len(list(self.matchingPositions(False)))
        return float(res) / all

    def matchingSequence(self, equalOnly=True):
        return MatchingSequence(self.seg_from.contig.seq, self.seg_to.contig.seq,
                                list(self.matchingPositions(equalOnly)))

    def contradicting(self, seg=None, tail_size = 500):
        # type: (Segment, int) -> bool
        if seg is None:
            seg = self.seg_to.contig.asSegment()
        if not self.seg_to.inter(seg):
            return False
        if self.seg_from.left >= tail_size and self.seg_to.left >= seg.left + tail_size:
            return True
        if self.seg_from.right <= len(self.seg_from.contig) - tail_size and self.seg_to.right <= seg.right - tail_size:
            return True
        return False

    def contradictingRTC(self,
                         seg=None):  # contradiction between read and consensus sequence. Stricter consensus condition
        # type: (Segment) -> bool
        if seg is None:
            seg = self.seg_to.contig.asSegment()
        if not self.seg_to.inter(seg):
            return False
        if self.seg_from.left >= 500 and self.seg_to.left >= seg.left + 50:
            return True
        if self.seg_from.right <= len(self.seg_from.contig) - 500 and self.seg_to.right <= seg.right - 50:
            return True
        return False

    def targetAsSegment(self, seg):
        # type: (Segment) -> AlignmentPiece
        seg_to = self.seg_to.contigAsSegment(seg)
        return AlignmentPiece(self.seg_from, seg_to, self.cigar)

    def queryAsSegment(self, seg):
        # type: (Segment) -> AlignmentPiece
        seg_from = self.seg_from.contigAsSegment(seg)
        return AlignmentPiece(seg_from, self.seg_to, self.cigar)

    def reduce(self, segment=None, query=None, target=None):
        # type: (Optional[Segment], Optional[Segment], Optional[Segment]) -> AlignmentPiece
        # TODO alignment reduce without additional translation. Same for alignment glue and composition
        if query is not None:
            return self.matchingSequence().reduceQuery(query.left, query.right).asAlignmentPiece(self.seg_from.contig,
                                                                                                 self.seg_to.contig)
        else:
            return self.matchingSequence().reduceTarget(target.left, target.right).asAlignmentPiece(
                self.seg_from.contig, self.seg_to.contig)

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.seg_from.save(handler)
        self.seg_to.save(handler)
        handler.newLine()
        handler.writeTokenLine(self.cigar)

    @staticmethod
    def load(handler, collection_from, collection_to):
        # type: (TokenReader, Any, Any) -> AlignmentPiece
        return AlignmentPiece(Segment.load(handler, collection_from), Segment.load(handler, collection_to), handler.readToken())

    # composes alignments A->B and B->C into alignment A->C
    def compose(self, other):
        # type: (AlignmentPiece) -> AlignmentPiece
        return self.matchingSequence(False).compose(other.matchingSequence(False)).asAlignmentPiece(self.seg_from.contig, other.seg_to.contig)

    # composes alignments A->B and A->C into alignment B->C
    def composeTargetDifference(self, other):
        # type: (AlignmentPiece) -> AlignmentPiece
        return self.matchingSequence(False).composeDifference(other.matchingSequence(False)).\
            asAlignmentPiece(self.seg_to.contig, other.seg_to.contig)

    # composes alignments A->B and C->B into alignment A->C
    def composeQueryDifference(self, other):
        # type: (AlignmentPiece) -> AlignmentPiece
        return self.matchingSequence(False).compose(other.matchingSequence(False).reverse()).\
            asAlignmentPiece(self.seg_from.contig, other.seg_from.contig)

    def reverse(self):
        # type: () -> AlignmentPiece
        return AlignmentPiece(self.seg_to, self.seg_from, easy_cigar.ReverseCigar(self.cigar))


class MatchingSequence:
    def __init__(self, seq_from, seq_to, matchingPositions):
        # type: (str, str, list[Tuple[int, int]]) -> MatchingSequence
        self.seq_from = seq_from
        self.seq_to = seq_to
        self.matches = matchingPositions

    def __str__(self):
        return str(self.matches)

    def common_from(self, other):
        # type: (MatchingSequence) -> Generator[Tuple[int, int]]
        assert self.seq_from == other.seq_from
        cur_self = 0
        cur_other = 0
        while cur_self < len(self) and cur_other < len(other):
            if self.matches[cur_self][0] < other.matches[cur_other][0]:
                cur_self += 1
            elif self.matches[cur_self][0] > other.matches[cur_other][0]:
                cur_other += 1
            else:
                yield (cur_self, cur_other)
                cur_self += 1
                cur_other += 1

    def common_to(self, other):
        # type: (MatchingSequence) -> Generator[Tuple[int, int]]
        assert self.seq_from == other.seq_from
        cur_self = 0
        cur_other = 0
        while cur_self < len(self) and cur_other < len(other):
            if self.matches[cur_self][1] < other.matches[cur_other][1]:
                cur_self += 1
            elif self.matches[cur_self][1] > other.matches[cur_other][1]:
                cur_other += 1
            else:
                yield (cur_self, cur_other)
                cur_self += 1
                cur_other += 1

    def reduceTarget(self, left, right):
        return MatchingSequence(self.seq_from, self.seq_to,
                                filter(lambda match: left <= match[1] < right, self.matches))

    def reduceQuery(self, left, right):
        return MatchingSequence(self.seq_from, self.seq_to,
                                filter(lambda match: left <= match[0] < right, self.matches))

    def mapDown(self, pos, roundDown = True):
        # type: (int, bool) -> Optional[int]
        if roundDown:
            for pos1, pos2 in self.matches[::-1]:
                if pos1 <= pos:
                    return pos2
            return None
        else:
            for pos1, pos2 in self.matches:
                if pos1 >= pos:
                    return pos2
            return None

    def mapUp(self, pos, roundDown = True):
        # type: (int, bool) -> Optional[int]
        if roundDown:
            for pos1, pos2 in self.matches[::-1]:
                if pos2 <= pos:
                    return pos1
            return None
        else:
            for pos1, pos2 in self.matches:
                if pos2 >= pos:
                    return pos1
            return None

    def asAlignmentPiece(self, contig_from, contig_to):
        # type: (NamedSequence, NamedSequence) -> AlignmentPiece
        return AlignmentPiece(self.SegFrom(contig_from), self.SegTo(contig_to), self.cigar())

    def cigar(self):
        # type: () -> str
        curm = 1
        res = []
        for m1, m2 in zip(self.matches[:-1], self.matches[1:]):
            d = (m2[0] - m1[0] - 1, m2[1] - m1[1] - 1)
            if d[0] == d[1]:
                curm += d[0] + 1
            else:
                if curm > 0:
                    res.append(str(curm))
                    res.append("M")
                res.append(str(abs(d[0] - d[1])))
                if d[0] > d[1]:
                    res.append("I")
                else:
                    res.append("D")
                curm = min(d[0], d[1]) + 1
        res.append(str(curm))
        res.append("M")
        return "".join(res)

    def inter(self, other):
        assert self.seq_from == other.seq_from and self.seq_to == other.seq_to
        cur_self = 0
        cur_other = 0
        matches = []
        while cur_self < len(self) and cur_other < len(other):
            if self.matches[cur_self] < other.matches[cur_other]:
                cur_self += 1
            elif self.matches[cur_self] < other.matches[cur_other]:
                cur_other += 1
            else:
                matches.append(self.matches[cur_self])
                cur_self += 1
                cur_other += 1
        return MatchingSequence(self.seq_from, self.seq_to, matches)

    def composeDifference(self, other):
        # type: (MatchingSequence) -> MatchingSequence
        matchings = [(self.matches[pos_self][1], other.matches[pos_other][1]) for pos_self, pos_other in
                     self.common_from(other)]
        return MatchingSequence(self.seq_to, other.seq_to, matchings)

    def compose(self, other):
        # type: (MatchingSequence) -> MatchingSequence
        return self.reverse().composeDifference(other)

    def concat(self, others):
        # type: (List[MatchingSequence]) -> MatchingSequence
        new_matches = list(itertools.chain(*[other.matches for other in others]))
        return MatchingSequence(self.seq_from, self.seq_to, self.matches + new_matches)

    def combine(self, others):
        # type: (List[MatchingSequence]) -> MatchingSequence
        if len(others) == 0:
            return self
        others = filter(lambda other: len(self.inter(other)) > 20, others)
        # print "Combining"
        # print self.matches
        matchings = sorted(list(itertools.chain(*others)))
        # print matchings
        res = [self.matches[0]]
        for matching in matchings:
            if matching[0] < self.matches[-1][0] and matching[1] < self.matches[-1][1] and matching[0] > res[-1][0] and \
                    matching[1] > res[-1][1]:
                res.append(matching)
        if len(self.matches) > 1:
            res.append(self.matches[-1])
        # print "Combining result:"
        # print res
        return MatchingSequence(self.seq_from, self.seq_to, res)

    def reverse(self):
        return MatchingSequence(self.seq_to, self.seq_from, map(lambda pair: (pair[1], pair[0]), self.matches))

    def SegFrom(self, contig):
        return Segment(contig, self.matches[0][0], self.matches[-1][0] + 1)

    def SegTo(self, contig):
        return Segment(contig, self.matches[0][1], self.matches[-1][1] + 1)

    def __getitem__(self, item):
        return self.matches[item]

    def __len__(self):
        return len(self.matches)


class AlignedRead(Contig):
    def __init__(self, rec, rc=None):
        # type: (NamedSequence, Optional[AlignedRead]) -> None
        if rec.id.startswith("contig_"):
            rec.id = rec.id[len("contig_"):]
        self.alignments = []  # type: list[AlignmentPiece]
        if rc is None:
            rc = AlignedRead(rec.RC(), self)
        Contig.__init__(self, rec.seq, rec.id, rc)
        self.rc = rc  # type: AlignedRead

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)
        handler.writeIntLine(len(self.alignments))
        for al in self.alignments:
            al.save(handler)

    @staticmethod
    def load(handler, collection):
        # type: (TokenReader, Any) -> AlignedRead
        id = handler.readToken()
        seq = handler.readToken()
        res = AlignedRead(NamedSequence(seq, id))
        n = handler.readInt()
        for i in range(n):
            res.alignments.append(AlignmentPiece.load(handler, res, collection))
        return res

    def alignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        self.sort()
        for al in self.alignments:
            if al.seg_to.inter(seg):
                yield al

    def AddSamAlignment(self, rec, contig):
        # type: (sam_parser.SAMEntryInfo, Contig) -> AlignmentPiece
        return AlignmentPiece.FromSamRecord(self, contig, rec)

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.alignments.append(al)
        self.rc.alignments.append(al.rc)

    def sort(self):
        self.alignments = sorted(self.alignments,
                                 key=lambda alignment: (alignment.seg_to.contig.id, alignment.seg_from.left))

    def __str__(self):
        self.sort()
        return str(self.id) + "(" + str(len(self.seq)) + ")" + "[" + ".".join(map(str, self.alignments)) + "]"

    def removeContig(self, contig):
        self.alignments = filter(lambda alignment: alignment.seg_to.contig.id != contig.id, self.alignments)
        self.rc.alignments = filter(lambda alignment: alignment.seg_to.contig.id != contig.rc.id, self.rc.alignments)

    def clean(self):
        self.alignments = []

    def contigAlignment(self, contig):
        # type: (Contig) -> Tuple[int,int]
        left = None
        right = None
        for alignment in self.alignments:
            if alignment.seg_to.contig.id != contig.id:
                continue
            if left is None or left > alignment.seg_to.left:
                left = alignment.seg_to.left
            if right is None or right < alignment.seg_to.right:
                right = alignment.seg_to.right
        return left, right

    def inter(self, other):
        # type: (Segment) -> bool
        for ap in self.alignments:
            if ap.seg_to.inter(other):
                return True
        return False

    def noncontradicting(self, seg):
        # type: (Segment) -> bool
        for al in self.alignments:
            if al.contradicting(seg):
                return False
        return True

    def changeTargets(self, contigs):
        # type: (ContigCollection) -> None
        new_alignments = []
        for al in self.alignments:
            if al.seg_to.contig in contigs:
                new_alignments.append(al.changeTargetContig(contigs[al.seg_to.contig.id]))
        self.alignments = new_alignments

    def contains(self, other):
        # type: (Segment) -> bool
        for ap in self.alignments:
            if ap.seg_to.contains(other):
                return True
        return False

    def contains_start(self, contig):
        # type: (Contig) -> bool
        return self.inter(Segment(contig, 0, 200))

    def contains_end(self, contig):
        # type: (Contig) -> bool
        return self.inter(Segment(contig, len(contig) - 200, len(contig)))

    def contigsAsSegments(self, seg_dict):
        # type: (Dict[int, Segment]) -> None
        # print "ContigsAsSegments", self, [str(i) + " " + str(seg_dict[i]) for i in seg_dict]
        self.rc.alignments = []
        alignments = []
        for al in self.alignments:
            if al.seg_to.contig.id in seg_dict:
                seg = seg_dict[al.seg_to.contig.id]
            elif al.seg_to.contig.rc.id in seg_dict:  # Fix it!!! This should never happen
                seg = seg_dict[al.seg_to.contig.rc.id].RC()
            else:
                assert False
            alignments.append(al.targetAsSegment(seg))
            self.rc.alignments.append(al.rc)
        self.alignments = alignments

    def mergeAlignments(self, other):
        # type: (AlignedRead) -> AlignedRead
        for al in other.alignments:
            new_al = al.changeQueryContig(self)
            self.alignments.append(new_al)
            self.rc.alignments.append(new_al.rc)
        return self

    def invalidate(self, seg):
        # type: (Segment) -> None
        self.alignments = filter(lambda al: not al.seg_to.inter(seg), self.alignments)
        self.rc.alignments = filter(lambda al: not al.seg_to.inter(seg.RC()), self.rc.alignments)

    def suffix(self, pos):
        # type: (int) -> Segment
        if pos < 0:
            pos = self.__len__() + pos
        if pos < 0:
            pos = 0
        if pos > len(self):
            pos = len(self)
        return Segment(self, pos, self.__len__())

    def prefix(self, len):
        # type: (int) -> Segment
        len = min(len, self.__len__())
        return Segment(self, 0, len)


class ReadCollection:
    def __init__(self, reads=None):
        # type: (Optional[Iterable[AlignedRead]]) -> None
        self.reads = dict()  # type: Dict[str, AlignedRead]
        if reads is not None:
            for read in reads:
                self.add(read)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine("ReadCollection")
        handler.writeTokens(self.reads.keys())
        handler.writeIntLine(len(list(UniqueList(self.reads.values()))))
        for read in UniqueList(self.reads.values()):
            read.save(handler)

    def load(self, handler, contigs):
        # type: (TokenReader, ContigCollection) -> None
        message = handler.readToken()
        assert message == "ReadCollection", message
        keys = set(handler.readTokens())
        n = handler.readInt()
        for i in range(n):
            read = AlignedRead.load(handler, contigs)
            self.add(read)
            if read.rc.id in keys:
                self.add(read.rc)

    def extend(self, other_collection):
        # type: (Iterable[AlignedRead]) -> ReadCollection
        for read in other_collection:
            self.add(read)
        return self

    def extendClean(self, other_collection):
        # type: (Iterable[NamedSequence]) -> ReadCollection
        for read in other_collection:
            if read.id not in self.reads:
                if basic.Reverse(read.id) in self.reads:
                    self.add(self.reads[basic.Reverse(read.id)].rc)
                else:
                    self.addNewRead(read)
        return self

    def addNewRead(self, rec):
        # type: (NamedSequence) -> AlignedRead
        new_id = str(rec.id).split()[0]
        if new_id in self.reads:
            return self.reads[rec.id]
        if basic.Reverse(new_id) in self.reads:
            self.add(self.reads[basic.Reverse(new_id)].rc)
            return self.reads[rec.id]
        read = AlignedRead(rec)
        self.add(read)
        return read

    def addNewAlignment(self, rec, contig):
        # type: (sam_parser.SAMEntryInfo, Contig) -> bool
        if rec.is_unmapped:
            return False
        rname = rec.query_name.split()[0]
        if rname.startswith("contig_"):
            rname = rname[len("contig_"):]
        if rname in self.reads:
            self.reads[rname].AddSamAlignment(rec, contig)
            return True
        elif basic.Reverse(rname) in self.reads:
            self.reads[basic.Reverse(rname)].rc.AddSamAlignment(rec, contig)
            return True
        return False

    def fillFromSam(self, sam, contigs):
        # type: (sam_parser.Samfile, ContigCollection) -> ReadCollection
        for rec in sam:
            self.addNewAlignment(rec, contigs[rec.tname])
        return self

    def add(self, read):
        # type: (AlignedRead) -> AlignedRead
        assert read not in self.reads or read == self.reads[read.id], str(read) + " " + str(self.reads[read.id])
        self.reads[read.id] = read
        return read

    def addAllRC(self):
        # type: () -> ReadCollection
        for read in self.reads.values():
            self.add(read.rc)
        return self

    def filter(self, condition):
        # type: (callable(AlignedRead)) -> ReadCollection
        res = ReadCollection()
        for read in self.reads.values():
            if condition(read):
                res.add(read)
        return res

    def copy(self):
        # type: () -> ReadCollection
        return self.filter(lambda read: True)

    def remove(self, read):
        # type: (AlignedRead) -> None
        if read in self.reads:
            del self.reads[read.id]

    def minus(self, other):
        # type: (ReadCollection) -> ReadCollection
        return self.filter(lambda read: read not in other)

    def minusBoth(self, other):
        # type: (ReadCollection) -> ReadCollection
        return self.filter(lambda read: read not in other and read.rc not in other)

    def minusAll(self, others):
        # type: (list[ReadCollection]) -> ReadCollection
        tmp = ReadCollection()
        for other in others:
            tmp.extend(other)
        return self.minus(tmp)

    def cap(self, other):
        # type: (ReadCollection) -> ReadCollection
        return self.filter(lambda read: read in other)

    def inter(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.inter(segment))

    def noncontradicting(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.noncontradicting(segment))

    def contain(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.contain(segment))

    def __iter__(self):
        # type: () -> Iterator[AlignedRead]
        return self.reads.values().__iter__()

    def __getitem__(self, read_id):
        # type: (str) -> AlignedRead
        return self.reads[read_id.split()[0]]

    def __contains__(self, read):
        # type: (AlignedRead) -> bool
        return read.id in self.reads

    def __len__(self):
        return len(self.reads)

    def print_fasta(self, hander):
        # type: (BinaryIO) -> None
        for read in self:
            SeqIO.write(read, hander, "fasta")

    def print_alignments(self, handler):
        # type: (file) -> None
        for read in self:
            handler.write(read.__str__() + "\n")

    def asSeqRecords(self):
        # type: () -> Generator[SeqIO.SeqRecord]
        for read in self.reads.values():
            yield common.seq_records.SeqRecord(read.seq, read.id)

    def loadFromFasta(self, handler):
        # type: (BinaryIO) -> ReadCollection
        for rec in SeqIO.parse_fasta(handler):
            new_read = self.add(AlignedRead(rec))
        return self

    def nontontradictingCopy(self, contig):
        res = ReadCollection()
        for read in self.inter(contig.asSegment()):
            add = True
            for al in read.alignments:
                if al.contradicting(contig.asSegment()):
                    add = False
            if add:
                assert read.rc not in res
                new_read = res.addNewRead(read)
                for al in read.alignments:
                    if al.seg_to.contig == contig:
                        new_read.addAlignment(al.changeQueryContig(new_read))
        return res

    def contigsAsSegments(self, seg_dict):
        # type: (Dict[str, Segment]) -> ReadCollection
        for read in UniqueList(self.reads.values()):
            read.contigsAsSegments(seg_dict)
        return self

    def changeTargets(self, contigs):
        # type: (ContigCollection) -> ReadCollection
        for read in self.reads.values():
            read.changeTargets(contigs)

    def cleanCopy(self, contigs, filter=lambda read: True):
        # type: (ContigCollection, Callable[[AlignedRead], bool]) -> ReadCollection
        res = ReadCollection()
        for read in self.reads.values():
            if not filter(read):
                continue
            if read.rc in res.reads:
                res.add(res.reads[read.rc.id].rc)
            else:
                res.addNewRead(read)
        return res

    def RC(self):
        res = ReadCollection()
        for read in self.reads.values():
            res.add(read.rc)

    def mergeAlignments(self, other):
        # type: (ReadCollection) -> ReadCollection
        for read in UniqueList(self):
            if read in other:
                read.mergeAlignments(other.reads[read.id])
            elif read.rc in other:
                read.rc.mergeAlignments(other.reads[read.rc.id])
        return self