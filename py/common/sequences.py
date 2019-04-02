import sys

import common.seq_records
from common.alignment_storage import AlignedRead
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence

sys.path.append("py")
from common import sam_parser, SeqIO, basic, params
from typing import Generator, Iterator, Dict, Optional, Union, Callable, Iterable, Any, BinaryIO

class EasyContig(NamedSequence):
    def __init__(self, seq, id, rc = None):
        # type: (str, str, Optional[EasyContig]) -> None
        NamedSequence.__init__(self, seq, id)
        if rc is None:
            rc = EasyContig(basic.RC(seq), basic.Reverse(id), self)
        self.rc = rc

class EasyContigStorage:
    def __init__(self, iter, add_rc = True):
        # type: (Iterable[EasyContig], bool) -> None
        self.items = dict() # type: Dict[str, EasyContig]
        self.add_rc = add_rc
        for item in iter:
            self.add(item)

    def add(self, item):
        # type: (EasyContig) -> None
        self.items[item.id] = item
        if self.add_rc:
            self.items[item.rc.id] = item.rc

    def __getitem__(self, item):
        # type: (str) -> Optional[EasyContig]
        if item in self.items:
            return self.items[item]
        return None


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


class Contig(EasyContig):
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
        EasyContig.__init__(self, seq, id, rc)

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
        assert 0 <= left <= right <= len(contig), str([0, left, right, len(contig), contig.id])
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
        return Segment(self.contig.rc, l - self.right, l - self.left)

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
        return self.contig.seq[self.left:self.right]

    def __str__(self):
        # type: () -> str
        if (self.right < len(self.contig) * 0.6 or (
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
        return Segment(seg.contig, self.left + seg.left, self.right + seg.left)

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


def UniqueList(sequences):
    # type: (Iterable[NamedSequence]) -> Generator[Any]
    visited = set()
    for seq in sequences:
        if seq.id in visited or basic.Reverse(seq.id) in visited:
            pass
        else:
            yield seq
            visited.add(seq.id)
