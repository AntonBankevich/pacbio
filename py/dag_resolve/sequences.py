import itertools
import sys

from common.SeqIO import NamedSequence
from dag_resolve import params

sys.path.append("py")
from common import sam_parser, SeqIO, basic
from typing import Generator, Iterator, Dict, Tuple, Optional, Union, Callable, Iterable, Any


class ContigCollection():
    def __init__(self, contigs_list = None):
        # type: (Optional[list[Contig]]) -> ContigCollection
        self.contigs = dict() # type: Dict[str, Contig]
        if contigs_list is not None:
            for contig in contigs_list:
                self.add(contig)

    def add(self, contig):
        # type: (Contig) -> None
        self.contigs[contig.id] = contig

    def filter(self, condition):
        # type: (callable(Contig)) -> ContigCollection
        res = ContigCollection()
        for contig in self.contigs.values():
            if condition(contig):
                res.add(contig)
        return res

    def incoming(self):
        # type: () -> ContigCollection
        return self.filter(lambda contig: "in" in contig.info.misc)

    def outgoing(self):
        # type: () -> ContigCollection
        return self.filter(lambda contig: "out" in contig.info.misc)

    def print_names(self, handler):
        # type: (file) -> None
        for contig in self:
            handler.write(str(contig.id) + " " + " ".join(contig.info.misc) + "\n")

    def __iter__(self):
        return self.contigs.values().__iter__()

    def __getitem__(self, contig_id):
        # type: (int) -> Contig
        return self.contigs[contig_id]

    def loadFromFasta(self, handler):
        # type: (file) -> ContigCollection
        for rec in SeqIO.parse_fasta(handler):
            self.add(Contig(rec.seq, basic.parseNegativeNumber(rec.id)))
        return self

    def print_fasta(self, handler):
        # type: (file) -> None
        for contig in self.contigs.values():
            contig.print_fasta(handler)

    def __len__(self):
        # type: () -> int
        return len(self.contigs)

    def __contains__(self, item):
        return item.id in self.contigs

class TmpInfo:
    def __init__(self, l):
        self.misc = l

class Contig(NamedSequence):
    def __init__(self, seq, id, info = None, rc = None):
        # type: (str, Union[str,int], Optional[TmpInfo], Optional[Contig]) -> Contig
        if info is None:
            info = []
        if isinstance(info, list):
            info = TmpInfo(info)
        self.info = info
        if rc is None:
            rc = Contig(basic.RC(seq), basic.Reverse(id), info, self)
        self.rc = rc
        NamedSequence.__init__(self, seq, id)

    def suffix(self, pos):
        # type: (int) -> Segment
        if pos < 0:
            pos = self.__len__() + pos
        return Segment(self, pos, self.__len__())

    def prefix(self, len):
        # type: (int) -> Segment
        len = min(len, self.__len__())
        return Segment(self, 0, len)

    def asSegment(self):
        # type: () -> Segment
        return Segment(self, 0, len(self))

    def print_fasta(self, handler):
        # type: (file) -> None
        SeqIO.write(self, handler, "fasta")

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def __eq__(self, other):
        return self.id == other.id

    def __ne__(self, other):
        return self.id != other.id

    def __str__(self):
        return str(self.id) + "(" + str(self.__len__()) + ")"


class Segment:
    def __init__(self, contig, left = None, right = None):
        # type: (Union[Contig, AlignedRead, NamedSequence], Optional[int], Optional[int]) -> Segment
        if left is None:
            left = 0
        if right is None:
            right = len(contig)
        assert 0 <= left <= right <= len(contig), str([0, left, right, len(contig), contig.id])
        self.contig = contig
        self.left = left
        self.right = right

    def RC(self):
        l = len(self.contig)
        return Segment(self.contig.rc, l - self.right, l - self.left)

    def precedes(self, other, delta = 0):
        # type: (Segment, int) -> bool
        return self.contig == other.contig and self.right <= other.left + delta

    def connects(self, other, delta = 0):
        # type: (Segment, int) -> bool
        return self.contig == other.contig and self.right - delta <= other.left <= self.right + delta

    def inter(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and not (self.right < other.left or self.left > other.right)

    def contains(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and self.left <= other.left and self.right >= other.right

    def subcontig(self):
        # type: () -> Contig
        return Contig(self.contig.seq[self.left:self.right], "(" + str(self.contig.id) + ")[" + str(self.left) + "," + str(self.right) + "]", self.contig.info)

    def merge(self, other):
        return Segment(self.contig, self.left, other.right)

    def shift(self, val):
        return Segment(self.contig, self.left + val, self.right + val)

    def Seq(self):
        return self.contig.seq[self.left:self.right]

    def __str__(self):
        # type: () -> str
        return str(self.contig.id) + "[" + str(self.left) + ":" + str(self.right) + "]"

    def __len__(self):
        return self.right - self.left

    def __eq__(self, other):
        return self.contig == other.contig and self.left == other.left and self.right == other.right

    def __ne__(self, other):
        return self.contig != other.contig or self.left != other.left or self.right != other.right


class AlignedRead(NamedSequence):
    def __init__(self, rec, rc = None):
        # type: (NamedSequence, Optional[AlignedRead]) -> AlignedRead
        NamedSequence.__init__(self, rec.seq, rec.id)
        self.alignments = [] # type: list[AlignmentPiece]
        if rc is None:
            rc = AlignedRead(rec.RC(), self)
        self.rc = rc

    def __len__(self):
        # type: () -> int
        return len(self.seq)

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
        seg_to = Segment(contig, rec.pos - 1, rec.pos - 1 + rec.alen)
        if not rec.rc:
            seg_from = Segment(self, ls, len(self.seq) - rs)
        else:
            seg_from = Segment(self.rc, rs, len(self.seq) - ls).RC()
            seg_to = seg_to.RC()
        piece = AlignmentPiece(seg_from, seg_to, new_cigar)
        self.alignments.append(piece)
        self.rc.alignments.append(piece.RC())
        return piece

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.alignments.append(al)
        self.rc.alignments.append(al.RC())

    def sort(self):
        self.alignments = sorted(self.alignments, key = lambda alignment: (alignment.seg_to.contig.id, alignment.seg_from.left))

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
            if al.noncontradicting(seg):
                return True
        return False


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

class ReadCollection:
    def __init__(self, contigs = ContigCollection()):
        # type: (ContigCollection) -> ReadCollection
        self.reads = dict() # type: Dict[str, AlignedRead]
        self.contigs = contigs

    def extend(self, other_collection):
        # type: (Iterable[AlignedRead]) -> ReadCollection
        for read in other_collection:
            self.add(read)
        return self

    def addNewRead(self, rec):
        # type: (NamedSequence) -> AlignedRead
        rec.id = str(rec.id).split()[0]
        read = AlignedRead(rec)
        self.reads[read.id] = read
        return read

    def addNewAlignment(self, rec):
        # type: (sam_parser.SAMEntryInfo) -> None
        if rec.is_unmapped:
            return
        rname = rec.query_name.split()[0]
        if rname not in self.reads:
            return
        self.reads[rname].AddSamAlignment(rec, self.contigs[int(rec.tname)])

    def loadFromSam(self, sam, filter = lambda rec: True):
        # type: (sam_parser.Samfile, Callable[[AlignedRead], bool]) -> ReadCollection
        last = [] # type: list[sam_parser.SAMEntryInfo]
        for new_rec in itertools.chain(sam, [None]):
            if new_rec is not None and (len(last) == 0 or new_rec.query_name == last[-1].query_name):
                last.append(new_rec)
            else:
                if len(last) == 0:
                    return self
                seq = None
                for rec in last:
                    if "H" not in rec.cigar and rec.seq != "*":
                        seq = rec.seq
                        break
                addRC = False
                addStraight = False
                assert seq is not None, last[0].query_name + " " + str(len(last))
                new_read = AlignedRead(SeqIO.SeqRecord(seq, last[0].query_name))
                for rec in last: #type: sam_parser.SAMEntryInfo
                    if rec.is_unmapped:
                        continue
                    assert int(rec.tname) in self.contigs.contigs
                    contig = self.contigs[int(rec.tname)]
                    alignment = new_read.AddSamAlignment(rec, contig)
                    if alignment.seg_to.contig in self.contigs:
                        addStraight = True
                    if alignment.seg_to.contig.rc in self.contigs:
                        addRC = True
                if addStraight and filter(new_read):
                    self.reads[new_read.id] = new_read
                if addRC and filter(new_read.rc):
                    self.reads[new_read.rc.id] = new_read.rc
                last = [new_rec]
        return self

    def add(self, read):
        # type: (AlignedRead) -> None
        self.reads[read.id] = read

    def filter(self, condition):
        # type: (callable(AlignedRead)) -> ReadCollection
        res = ReadCollection(self.contigs)
        for read in self.reads.values():
            if condition(read):
                res.add(read)
        return res

    def copy(self):
        # type: () -> ReadCollection
        return self.filter(lambda read: True)

    def remove(self, read):
        # type: (AlignedRead) -> None
        if read.id in self.reads:
            del self.reads[read.id]

    def minus(self, other):
        # type: (ReadCollection) -> ReadCollection
        return self.filter(lambda read: read.id not in other)

    def minusAll(self, others):
        # type: (list[ReadCollection]) -> ReadCollection
        tmp = ReadCollection(self.contigs)
        for other in others:
            tmp.extend(other)
        return self.minus(tmp)

    def cap(self, other):
        # type: (ReadCollection) -> ReadCollection
        return self.filter(lambda read: read.id in other)

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

    def __contains__(self, item):
        return self.reads.__contains__(item)

    def __len__(self):
        return len(self.reads)

    def print_fasta(self, hander):
        # type: (file) -> None
        for read in self:
            SeqIO.write(read, hander, "fasta")

    def print_alignments(self, handler):
        # type: (file) -> None
        for read in self:
            handler.write(read.__str__() + "\n")

    def asSeqRecords(self):
        # type: () -> Generator[SeqIO.SeqRecord]
        for read in self.reads.values():
            yield SeqIO.SeqRecord(read.seq, read.id)

    def loadFromFasta(self, handler):
        # type: (file) -> ReadCollection
        for rec in SeqIO.parse_fasta(handler):
            self.add(AlignedRead(rec))
        return self

    def cleanCopy(self, contigs):
        res = ReadCollection(contigs)
        for read in self.reads.values():
            if read.rc.id in res.reads:
                res.add(res.reads[read.rc.id].rc)
            else:
                res.addNewRead(read)
        return res


class Consensus:
    def __init__(self, seq, cov):
        # type: (str, list[int]) -> Consensus
        self.seq = seq
        self.cov = cov

    def suffix(self, pos):
        if pos < 0:
            pos = self.__len__() + pos
        return Consensus(self.seq[pos:], self.cov[pos:])

    def printQuality(self, handler, cov_threshold = params.reliable_coverage):
        # type: (file, int) -> None
        for c, a in zip(self.seq, self.cov):
            if a < cov_threshold:
                handler.write(c.lower())
            else:
                handler.write(c.upper())
        handler.write("\n")

    def printCoverage(self, handler, step):
        # type: (file, int) -> None
        for i, val in list(enumerate(self.cov))[::step]:
            handler.write(str((i, val)) + "\n")

    def cut(self, cov_threshold = params.reliable_coverage, length = None):
        l = len(self.seq)
        while l > 0 and self.cov[l - 1] < cov_threshold:
            l -= 1
        if length is not None:
            l = min(l, length)
        return Consensus(self.seq[:l], self.cov[:l])

    def __len__(self):
        return len(self.seq)

    def RC(self):
        return Consensus(basic.RC(self.seq), self.cov[::-1])

    def merge(self, other, pos = None):
        # type: (Consensus, Optional[int]) -> None
        if pos is None:
            pos = len(self)
        self.seq = self.seq[:pos] + other.seq
        self.cov = self.cov[:pos] + other.cov


class MatchingSequence:
    def __init__(self, seq_from, seq_to, matchingPositions):
        # type: (str, str, list[Tuple[int, int]]) -> MatchingSequence
        self.seq_from = seq_from
        self.seq_to = seq_to
        self.matches = matchingPositions

    def common(self, other):
        # type: (MatchingSequence) -> Generator(Tuple[int, int])
        assert self.seq_from == other.seq_from
        cur_self = 0
        cur_other = 0
        while cur_self < len(self) and cur_other < len(other):
            if self.matches[cur_self][0] < other.matches[cur_other][0]:
                cur_self += 1
            elif self.matches[cur_self][0] < other.matches[cur_other][0]:
                cur_other += 1
            else:
                yield (cur_self, cur_other)
                cur_self += 1
                cur_other += 1

    def reduce(self, left, right):
        return MatchingSequence(self.seq_from, self.seq_to, filter(lambda match: left <= match[0] <= right, self.matches))

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

    def compose(self, other):
        # type: (MatchingSequence) -> MatchingSequence
        matchings = [(self.matches[pos_self][1], self.matches[pos_other][1]) for pos_self, pos_other in self.common(other)]
        return MatchingSequence(self.seq_to, other.seq_to, matchings)

    def combine(self, others):
        # type: (list[MatchingSequence]) -> MatchingSequence
        others = filter(lambda other: len(self.inter(other)) > 20, others)
        matchings = sorted(list(itertools.chain(others)))
        res = [self.matches[0]]
        for matching in matchings[1:]:
            if matching[0] < self.matches[-1][0] and matching[1] < self.matches[-1][1] and matching[0] > res[-1][0] and matching[1] > res[-1][1]:
                res.append(matching)
        if len(self.matches) > 1:
            res.append(self.matches[-1])
        return MatchingSequence(self.seq_from, self.seq_to, res)


    def __len__(self):
        return len(self.matches)


class AlignmentPiece:
    def __init__(self, seg_from, seg_to, cigar = None):
        # type: (Segment, Segment, str) -> AlignmentPiece
        self.seg_from = seg_from
        self.seg_to = seg_to
        assert cigar != "X"
        assert cigar.find("H") == -1 and cigar.find("S") == -1
        self.cigar = cigar

    def __str__(self):
        # type: () -> str
        return "(" + str(self.seg_from) + "->" + str(self.seg_to) + ")"

    def precedes(self, other, delta = 0):
        # type: (AlignmentPiece, int) -> bool
        return self.seg_from.precedes(other.seg_from, delta) and self.seg_to.precedes(other.seg_to, delta)

    def connects(self, other, delta = 0):
        # type: (AlignmentPiece, int) -> bool
        return self.seg_from.connects(other.seg_from, delta) and self.seg_to.connects(other.seg_to, delta)

    def contains(self, other):
        return self.seg_from.contains(other) and self.seg_to.contains(other)

    def merge(self, other):
        return AlignmentPiece(self.seg_from.merge(other.seg_from), self.seg_to.merge(other.seg_to), self.cigar + other.cigar)

    def __eq__(self, other):
        return self.seg_from == other.seg_from and self.seg_to == other.seg_to and self.cigar == other.cigar

    def __ne__(self, other):
        return self.seg_from != other.seg_from or self.seg_to != other.seg_to or self.cigar != other.cigar

    def __len__(self):
        return len(self.seg_from)

    def RC(self):
        return AlignmentPiece(self.seg_from.RC(), self.seg_to.RC(), sam_parser.RCCigar(self.cigar))

    def replaceQuery(self, seq):
        return AlignmentPiece()

    def matchingPositions(self, equalOnly = False):
        # type: (bool) -> list[Tuple[int, int]]
        if self.cigar == "=":
            for i in range(len(self.seg_from)):
                yield (self.seg_from.left + i, self.seg_to.left + i)
            return
        cur_tar = self.seg_from.left
        cur_query = self.seg_to.right
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
                cur_query += n
            elif c in "I":
                cur_tar += n

    def MatchingSequence(self):
        return MatchingSequence(self.seg_from.contig.seq, self.seg_to.contig.seq, self.matchingPositions(True))

    def percentIdentity(self):
        res = len(list(self.matchingPositions(True)))
        return float(res) / max(len(self.seg_from), len(self.seg_to))

    def noncontradicting(self, seg):
        # type: (Segment) -> bool
        if not self.seg_to.inter(seg):
            return True
        return (self.seg_from.left < 500 or self.seg_to.left < seg.left + 500) and \
                (self.seg_from.right > len(self.seg_from.contig) - 500 or self.seg_to.right > seg.right - 500)


def UniqueList(sequences):
    # type: (Iterable[NamedSequence]) -> Generator[Any]
    visited = set()
    for seq in sequences:
        if seq.id in visited or basic.Reverse(seq.id) in visited:
            pass
        else:
            yield seq
            visited.add(seq.id)
