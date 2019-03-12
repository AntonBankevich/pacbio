import itertools

from typing import Generator, Tuple, Optional, Any, List, Dict

from common import sam_parser
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, Contig, ContigCollection
from dag_resolve import params


class AlignmentPiece:
    def __init__(self, seg_from, seg_to, cigar, rc = None):
        # type: (Segment, Segment, str, Optional[AlignmentPiece]) -> None
        self.seg_from = seg_from #type: Segment
        self.seg_to = seg_to #type: Segment
        assert cigar != "X"
        assert cigar.find("H") == -1 and cigar.find("S") == -1
        self.cigar = cigar
        if params.assert_pi:
            pi = self.matchingPercentIdentity()
            if pi < 0.05:
                print seg_from.Seq()
                print seg_to.Seq()
            assert pi > 0.5, str(self)
        if rc is None:
            self.rc = AlignmentPiece(seg_from.RC(), seg_to.RC(), sam_parser.RCCigar(self.cigar), self)
        else:
            self.rc = rc

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

    def changeQuery(self, read):
        return AlignmentPiece(self.seg_from.changeContig(read), self.seg_to, self.cigar)

    def changeTarget(self, contig):
        return AlignmentPiece(self.seg_from, self.seg_to.changeContig(contig), self.cigar)

    def precedes(self, other, delta=0):
        # type: (AlignmentPiece, int) -> bool
        return self.seg_from.precedes(other.seg_from, delta) and self.seg_to.precedes(other.seg_to, delta)

    def connects(self, other, delta=0):
        # type: (AlignmentPiece, int) -> bool
        return self.seg_from.connects(other.seg_from, delta) and self.seg_to.connects(other.seg_to, delta)

    def contains(self, other):
        # type: (AlignmentPiece) -> bool
        return self.seg_from.contains(other.seg_from) and self.seg_to.contains(other.seg_to)

    def merge(self, other):
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
        if self.cigar == "=":
            for i in range(len(self.seg_from)):
                yield (self.seg_from.left + i, self.seg_to.left + i)
            return
        cur_query = self.seg_from.left
        cur_tar = self.seg_to.left
        for n, c in sam_parser.pattern.findall(self.cigar):
            if n:
                n = int(n)
            else:
                n = 1
            if c == 'M':
                for i in range(n):
                    if not equalOnly or self.seg_from.contig[cur_query] == self.seg_to.contig[cur_tar]:
                        yield (cur_query, cur_tar)
                    cur_tar += 1
                    cur_query += 1
            elif c in 'DPN':
                cur_tar += n
            elif c in "I":
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

    def contradicting(self, seg=None):
        # type: (Segment) -> bool
        if seg is None:
            seg = self.seg_to.contig.asSegment()
        if not self.seg_to.inter(seg):
            return False
        if self.seg_from.left >= 500 and self.seg_to.left >= seg.left + 500:
            return True
        if self.seg_from.right <= len(self.seg_from.contig) - 500 and self.seg_to.right <= seg.right - 500:
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
        # type: (Segment) -> None
        # print "ContigAsSegment", self, seg
        self.seg_to = self.seg_to.contigAsSegment(seg)

    def queryAsSegment(self, seg):
        # type: (Segment) -> None
        self.seg_from = self.seg_from.contigAsSegment(seg)

    def reduce(self, segment=None, query=None, target=None):
        # type: (Optional[Segment], Optional[Segment], Optional[Segment]) -> AlignmentPiece
        if segment is not None:
            if segment.contig == self.seg_from.contig:
                query = segment
            else:
                target = segment
        assert (query is None) != (target is None), str(query) + " " + str(target)
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

    def compose(self, other):
        # type: (AlignmentPiece) -> AlignmentPiece
        return self.matchingSequence(False).compose(other.matchingSequence(False)).asAlignmentPiece(self.seg_from.contig, other.seg_to.contig)


class MatchingSequence:
    def __init__(self, seq_from, seq_to, matchingPositions):
        # type: (str, str, list[Tuple[int, int]]) -> MatchingSequence
        self.seq_from = seq_from
        self.seq_to = seq_to
        self.matches = matchingPositions

    def common(self, other):
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

    def reduceTarget(self, left, right):
        return MatchingSequence(self.seq_from, self.seq_to,
                                filter(lambda match: left <= match[1] < right, self.matches))

    def reduceQuery(self, left, right):
        return MatchingSequence(self.seq_from, self.seq_to,
                                filter(lambda match: left <= match[0] < right, self.matches))

    def asAlignmentPiece(self, contig_from, contig_to):
        # type: (NamedSequence, NamedSequence) -> AlignmentPiece
        return AlignmentPiece(self.SegFrom(contig_from), self.SegTo(contig_to), self.cigar())

    def cigar(self):
        # type: () -> str
        curm = 0
        res = []
        for m1, m2 in zip(self.matches[:-1], self.matches[1:]):
            d = (m2[0] - m1[0] - 1, m2[1] - m1[1] - 1)
            if d[0] == d[1]:
                curm += d[0] + 1
            else:
                curm += min(d[0], d[1])
                res.append(str(curm))
                res.append("M")
                res.append(str(abs(d[0] - d[1])))
                if d[0] > d[1]:
                    res.append("I")
                else:
                    res.append("D")
                curm = 1
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
                     self.common(other)]
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
        return Segment(contig, self.matches[0][0], self.matches[-1][0])

    def SegTo(self, contig):
        return Segment(contig, self.matches[0][1], self.matches[-1][1])

    def __getitem__(self, item):
        return self.matches[item]

    def __len__(self):
        return len(self.matches)


class AlignedRead(NamedSequence):
    def __init__(self, rec, rc=None):
        # type: (NamedSequence, Optional[AlignedRead]) -> None
        if rec.id.startswith("contig_"):
            rec.id = rec.id[len("contig_"):]
        NamedSequence.__init__(self, rec.seq, rec.id)
        self.alignments = []  # type: list[AlignmentPiece]
        if rc is None:
            rc = AlignedRead(rec.RC(), self)
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
            al.seg_from.save(handler)
            al.seg_to.save(handler)
            handler.writeTokenLine(al.cigar)

    @staticmethod
    def load(handler, collections):
        # type: (TokenReader, List) -> AlignedRead
        id = handler.readToken()
        seq = handler.readToken()
        res = AlignedRead(NamedSequence(seq, id))
        n = handler.readInt()
        for i in range(n):
            seg_from = Segment.load(handler, res)
            seg_to = Segment.load(handler, collections)
            cigar = handler.readToken()
            res.alignments.append(AlignmentPiece(seg_from, seg_to, cigar))
        return res

    def alignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        self.sort()
        for al in self.alignments:
            if al.seg_to.inter(seg):
                yield al

    def AddSamAlignment(self, rec, contig):
        # type: (sam_parser.SAMEntryInfo, Contig) -> AlignmentPiece
        cigar_list = list(sam_parser.CigarToList(rec.cigar))
        ls = 0
        rs = 0
        if cigar_list[0][0] in "HS":
            ls = cigar_list[0][1]
        if cigar_list[-1][0] in "HS":
            rs = cigar_list[-1][1]
        new_cigar = []
        for s, num in cigar_list:
            if s not in "HS":
                new_cigar.append(num)
                new_cigar.append(s)
        new_cigar = "".join(map(str, new_cigar))
        seg_from = Segment(self, ls, len(self.seq) - rs)
        seg_to = Segment(contig, rec.pos - 1, rec.pos - 1 + rec.alen)
        if rec.rc:
            seg_from = Segment(self.rc, ls, len(self.seq) - rs).RC()
            seg_to = seg_to.RC()
            new_cigar = sam_parser.RCCigar(new_cigar)
        piece = AlignmentPiece(seg_from, seg_to, new_cigar)
        seg_from.contig.alignments.append(piece)
        seg_from.contig.rc.alignments.append(piece.rc)
        # self.rc.alignments.append(piece.RC())
        # if piece.percentIdentity() < 0.3:
        #     print piece
        #     print rec.pos, rec.rc, rec.query_name, rec.tname, rec.alen
        #     print rec.cigar
        #     print piece.seg_from.Seq()
        #     print piece.seg_to.Seq()
        #     print "\n".join(map(str, piece.matchingPositions()))
        #     assert False
        return piece

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
                new_alignments.append(al.changeTarget(contigs[al.seg_to.contig.id]))
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
        for al in self.alignments:
            if al.seg_to.contig.id in seg_dict:
                seg = seg_dict[al.seg_to.contig.id]
            elif al.seg_to.contig.rc.id in seg_dict:  # Fix it!!! This should never happen
                seg = seg_dict[al.seg_to.contig.rc.id].RC()
            else:
                assert False
            al.targetAsSegment(seg)
            self.rc.alignments.append(al.rc)

    def mergeAlignments(self, other):
        # type: (AlignedRead) -> AlignedRead
        for al in other.alignments:
            new_al = al.changeQuery(self)
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


