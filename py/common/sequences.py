import itertools
import sys

import common.seq_records
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from dag_resolve import params

sys.path.append("py")
from common import sam_parser, SeqIO, basic
from typing import Generator, Iterator, Dict, Tuple, Optional, Union, Callable, Iterable, Any, BinaryIO, List


class ContigCollection():
    def __init__(self, contigs_list=None):
        # type: (Optional[Iterable[Contig]]) -> ContigCollection
        self.contigs = dict()  # type: Dict[str, Contig]
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
        # type: (str) -> Contig
        return self.contigs[contig_id]

    def loadFromFasta(self, handler, num_names=True):
        # type: (BinaryIO, bool) -> ContigCollection
        for rec in SeqIO.parse_fasta(handler):
            if num_names:
                self.add(Contig(rec.seq, str(basic.parseNegativeNumber(rec.id))))
            else:
                self.add(Contig(rec.seq, rec.id))
        return self

    def print_fasta(self, handler):
        # type: (file) -> None
        for contig in self.contigs.values():
            contig.print_fasta(handler)

    def RC(self):
        return ContigCollection(map(lambda contig: contig.rc, self.contigs.values()))

    def __len__(self):
        # type: () -> int
        return len(self.contigs)

    def __contains__(self, item):
        return item.id in self.contigs

    def containsKey(self, key):
        return key in self.contigs

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine("ContigCollection")
        handler.writeTokens(self.contigs.keys())
        handler.writeIntLine(len(list(UniqueList(self.contigs.values()))))
        for contig in UniqueList(self.contigs.values()):
            contig.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        message = handler.readToken()
        assert message == "ContigCollection", message
        keys = set(handler.readTokens())
        n = handler.readInt()
        for i in range(n):
            contig = Contig.load(handler)
            self.add(contig)
            if contig.rc.id in keys:
                self.add(contig.rc)




class TmpInfo:
    def __init__(self, l):
        self.misc = l


class Contig(NamedSequence):
    def __init__(self, seq, id, info=None, rc=None):
        # type: (str, str, Optional[TmpInfo], Optional[Contig]) -> None
        if info is None:
            info = []
        if isinstance(info, list):
            info = TmpInfo(info)
        self.info = info
        if rc is None:
            rc = Contig(basic.RC(seq), basic.Reverse(id), info, self)
        self.rc = rc
        NamedSequence.__init__(self, seq, id)

    def segment(self, left, right):
        return Segment(self, left, right)

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

    def asSegment(self):
        # type: () -> Segment
        return Segment(self, 0, len(self))

    def print_fasta(self, handler):
        # type: (file) -> None
        SeqIO.write(self, handler, "fasta")

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def __str__(self):
        return str(self.id) + "(" + str(self.__len__()) + ")"

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)

    @staticmethod
    def load(handler):
        # type: (TokenReader) -> Contig
        id = handler.readToken()
        seq = handler.readToken()
        return Contig(seq, id)


class Segment:
    def __init__(self, contig, left=None, right=None):
        # type: (Union[Contig, AlignedRead, NamedSequence], Optional[int], Optional[int]) -> Segment
        if left is None:
            left = 0
        if right is None:
            right = len(contig)
        assert 0 <= left + contig.zero_pos <= right + contig.zero_pos <= len(contig), str([0, left, right, len(contig), contig.id])
        self.contig = contig
        self.left = left
        self.right = right

    def cap(self, other):
        # type: (Segment) -> Segment
        assert self.inter(other)
        return Segment(self.contig, max(self.left, other.left), min(self.right, other.right))

    def cup(self, other):
        # type: (Segment) -> Segment
        assert self.inter(other)
        return Segment(self.contig, min(self.left, other.left), max(self.right, other.right))

    def RC(self):
        # type: () -> Segment
        l = len(self.contig)
        return Segment(self.contig.rc,
                       l - self.right - self.contig.zero_pos - self.contig.rc.zero_pos,
                       l - self.left - self.contig.zero_pos - self.contig.rc.zero_pos)

    def asNamedSequence(self):
        # type: () -> NamedSequence
        return NamedSequence(self.Seq(), self.contig.id + "[" + str(self.left) + "," + str(self.right) + "]")

    def precedes(self, other, delta=0):
        # type: (Segment, int) -> bool
        return self.contig == other.contig and self.right <= other.left + delta

    def connects(self, other, delta=0):
        # type: (Segment, int) -> bool
        return self.contig == other.contig and self.right - delta <= other.left <= self.right + delta

    def inter(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and not (self.right <= other.left or self.left >= other.right)

    def contains(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and self.left <= other.left and self.right >= other.right

    def subcontig(self):
        # type: () -> Contig
        return Contig(self.Seq(),
                      "(" + self.contig.id + ")[" + str(self.left) + "," + str(self.right) + "]", self.contig.info)

    def merge(self, other):
        return Segment(self.contig, self.left, other.right)

    def shift(self, val):
        return Segment(self.contig, self.left + val, self.right + val)

    def Seq(self):
        # type: () -> str
        return self.contig.seq[self.left + self.contig.zero_pos:self.right + self.contig.zero_pos]

    def __str__(self):
        # type: () -> str
        if self.contig.zero_pos == 0 and (self.right < len(self.contig) * 0.6 or (
                self.right < 0.9 * len(self.contig) and len(self.contig) - self.right > 500)):
            return str(self.contig.id) + "[" + str(self.left) + ":" + str(self.right) + "]"
        else:
            return str(self.contig.id) + "[" + str(self.left) + ":" + str(len(self.contig)) + "-" + str(
                len(self.contig) - self.right) + "]"

    def __len__(self):
        return self.right - self.left

    def __eq__(self, other):
        return self.contig == other.contig and self.left == other.left and self.right == other.right

    def __ne__(self, other):
        return self.contig != other.contig or self.left != other.left or self.right != other.right

    def changeContig(self, contig):
        # type: (Segment) -> Segment
        return Segment(contig, self.left, self.right)

    def expand(self, range):
        # type: (int) -> Segment
        return Segment(self.contig, max(self.left - range, 0), min(self.right + range, len(self.contig)))

    def copy(self):
        # type: () -> Segment
        return Segment(self.contig, self.left, self.right)

    def contigAsSegment(self, seg):
        # type: (Segment) -> Segment
        return Segment(seg.contig, self.left + seg.left + self.contig.zero_pos, self.right + seg.left + self.contig.zero_pos)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.contig.id)
        handler.writeToken(str(self.left))
        handler.writeToken(str(self.right))

    @staticmethod
    def load(handler, contig_or_collection):
        # type: (TokenReader, Any) -> Segment
        id = handler.readToken()
        if isinstance(contig_or_collection, NamedSequence) and contig_or_collection.id == id:
            contig = contig_or_collection
        else:
            contig = contig_or_collection[id]
        assert contig is not None and contig.id == id, id
        return Segment(contig, int(handler.readToken()), int(handler.readToken()))



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
        seg_from.contig.rc.alignments.append(piece.RC())
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
        self.rc.alignments.append(al.RC())

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
            self.rc.alignments.append(al.RC())

    def mergeAlignments(self, other):
        # type: (AlignedRead) -> AlignedRead
        for al in other.alignments:
            self.alignments.append(al.changeQuery(self))
            self.rc.alignments.append(al.RC().changeQuery(self.rc))
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


def loadFromSam(reads, sam, contigs, filter=lambda rec: True):
    # type: (ReadCollection, sam_parser.Samfile, ContigCollection, Callable[[AlignedRead], bool]) -> ReadCollection
    recs = list(sam)
    for new_rec in recs:
        if "H" not in new_rec.cigar and new_rec.seq != "*":
            if new_rec.rc:
                seq = basic.RC(new_rec.seq).upper()
            else:
                seq = new_rec.seq.upper()
            read = AlignedRead(NamedSequence(seq, new_rec.query_name))
            reads.reads[read.id] = read
    for new_rec in recs:
        read = reads.reads[new_rec.query_name]
        addRC = False
        if new_rec.is_unmapped:
            continue
        assert new_rec.tname in contigs.contigs, str(new_rec.tname) + " " + str(list(contigs.contigs.keys()))
        contig = contigs[new_rec.tname]
        alignment = read.AddSamAlignment(new_rec, contig)
        if alignment.seg_to.contig.rc in contigs:
            addRC = True
        if addRC and filter(read.rc):
            reads.reads[read.rc.id] = read.rc
    return reads

class ReadCollection:
    # type: (ContigCollection, Optional[Iterable[AlignedRead]]) -> ReadCollection
    def __init__(self, reads=None):
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
                        new_read.addAlignment(al.changeQuery(new_read))
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


class Consensus:
    def __init__(self, seq, cov):
        # type: (str, list[int]) -> Consensus
        self.seq = seq
        self.cov = cov

    def suffix(self, pos):
        if pos < 0:
            pos = self.__len__() + pos
        return Consensus(self.seq[pos:], self.cov[pos:])

    def printQuality(self, handler, cov_threshold=params.reliable_coverage):
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

    def cut(self, cov_threshold=params.reliable_coverage, length=None):
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

    def merge(self, other, pos=None):
        # type: (Consensus, Optional[int]) -> Consensus
        if pos is None:
            pos = len(self)
        return Consensus(self.seq[:pos] + other.seq, self.cov[:pos] + other.cov)


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


class AlignmentPiece:
    def __init__(self, seg_from, seg_to, cigar=None):
        # type: (Segment, Segment, str) -> None
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

    def RC(self):
        return AlignmentPiece(self.seg_from.RC(), self.seg_to.RC(), sam_parser.RCCigar(self.cigar))

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


def UniqueList(sequences):
    # type: (Iterable[NamedSequence]) -> Generator[Any]
    visited = set()
    for seq in sequences:
        if seq.id in visited or basic.Reverse(seq.id) in visited:
            pass
        else:
            yield seq
            visited.add(seq.id)
