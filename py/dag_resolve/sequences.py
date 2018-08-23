import sys
sys.path.append("py")
from common import sam_parser, SeqIO, basic
from typing import Generator


class ContigCollection():
    def __init__(self, contigs_list = []):
        # type: (list[Contig]) -> ContigCollection
        self.contigs = dict() # type: dict[str, Contig]
        for contig in contigs_list:
            self.add(contig)
        self.main = None

    def add(self, contig):
        # type: (Contig) -> None
        if "main" in contig.info.misc:
            self.main = contig
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

    def main(self):
        # type: () -> Contig
        return self.main

    def print_names(self, handler):
        # type: (file) -> object
        for contig in self:
            handler.write(contig.id + " " + " ".join(contig.info.misc) + "\n")

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

class TmpInfo:
    def __init__(self, l):
        self.misc = l

class Contig:
    def __init__(self, seq, id, info = []):
        # type: (str, int, TmpInfo) -> Contig
        self.seq = seq.upper()
        self.id = id
        if isinstance(info, list):
            info = TmpInfo(info)
        self.info = info

    def suffix(self, len):
        # type: (int) -> Segment
        len = min(len, self.__len__())
        return Segment(self, self.__len__() - len, self.__len__())

    def prefix(self, len):
        # type: (int) -> Segment
        len = min(len, self.__len__())
        return Segment(self, 0, len)

    def as_segment(self):
        # type: () -> Segment
        return Segment(self, 0, len(self))

    def print_fasta(self, handler):
        # type: (file) -> None
        SeqIO.write(self, handler, "fasta")

    def __len__(self):
        # type: () -> int
        return len(self.seq)

class Segment:
    def __init__(self, contig, left, right):
        # type: (Contig, int, int) -> Segment
        self.contig = contig
        self.left = left
        self.right = right

    def inter(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and not (self.right < other.left or self.left > other.right)

    def contains(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and self.left <= other.left and self.right >= other.right

    def subcontig(self):
        # type: () -> Contig
        return Contig(self.contig.seq[self.left:self.right], str(self.contig.id) + "[" + str(self.left) + "," + str(self.right) + "]", self.contig.info)

    def __str__(self):
        # type: () -> str
        return str(self.contig.id) + "[" + str(self.left) + ":" + str(self.right) + "]"


class AlignmentPiece:
    def __init__(self, seg_from, seg_to, rc, cigar = None):
        # type: (Segment, Segment, bool, str) -> AlignmentPiece
        self.seg_from = seg_from
        self.seg_to = seg_to
        self.rc = rc
        self.cigar = cigar

    def __str__(self):
        # type: () -> str
        return "(" + str(self.seg_from) + "->" + str(self.seg_to) + ")"

class Read:
    def __init__(self, rec):
        # type: (SeqIO.SeqRecord) -> Read
        self.id = rec.id
        self.seq = rec.seq.upper()
        self.alignments = [] # type: list[AlignmentPiece]

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def AddSamAlignment(self, rec, contig):
        # type: (sam_parser.SAMEntryInfo, Contig) -> None
        cigar_list = list(sam_parser.CigarToList(rec.cigar))
        ls = 0
        rs = 0
        if cigar_list[0][0] in "HS":
            ls = cigar_list[0][1]
        if cigar_list[-1][0] in "HS":
            rs = cigar_list[-1][1]
        if not rec.rc:
            self.alignments.append(AlignmentPiece(Segment(self, ls, len(self.seq) - rs), Segment(contig, rec.pos, rec.pos + rec.alen), rec.rc, rec.cigar))
        else:
            self.alignments.append(AlignmentPiece(Segment(self, len(self.seq) - ls, rs), Segment(contig, rec.pos, rec.pos + rec.alen), rec.rc, rec.cigar))

    def __str__(self):
        return "Read:" + str(self.id) + "[" + ".".join(map(str, self.alignments)) + "]"


    def inter(self, other):
        # type: (Segment) -> bool
        for ap in self.alignments:
            if ap.seg_to.inter(other):
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
        self.reads = dict() # type: dict[str, Read]
        self.contigs = contigs

    def extend(self, other_collection):
        # type: (ReadCollection) -> None
        for read in other_collection.reads.values():
            self.add(read)
        return self

    def addNewRead(self, rec):
        # type: (SeqIO.SeqRecord) -> Read
        rec.id = rec.id.split()[0]
        read = Read(rec)
        self.reads[rec.id] = read
        return read

    def addNewAlignment(self, rec):
        # type: (sam_parser.SAMEntryInfo) -> None
        if rec.is_unmapped:
            return
        rname = rec.query_name.split()[0]
        if rname not in self.reads:
            return
        self.reads[rname].AddSamAlignment(rec, self.contigs[int(rec.tname)])

    def loadFromSam(self, sam):
        # type: (sam_parser.Samfile) -> None
        for rec in sam:
            if rec.is_unmapped:
                continue
            if rec.query_name not in self.reads:
                new_read = Read(SeqIO.SeqRecord(None, rec.query_name))
                self.add(new_read)
            else:
                new_read = self.reads[rec.query_name]
            if not rec.secondary and new_read.seq is None:
                if rec.rc:
                    new_read.seq = basic.RC(rec.seq)
                else:
                    new_read.seq = rec.seq
            self.addNewAlignment(rec)
        for read in self.reads.values():
            assert read.seq is not None, "Read " + read.id + " was aligned but had no primary alignments"

    def add(self, read):
        # type: (Read) -> None
        self.reads[read.id] = read

    def filter(self, condition):
        # type: (callable(Read)) -> ReadCollection
        res = ReadCollection(self.contigs)
        for read in self.reads.values():
            if condition(read):
                res.add(read)
        return res

    def copy(self):
        # type: () -> ReadCollection
        return self.filter(lambda read: True)

    def remove(self, read):
        # type: (Read) -> None
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

    def inter(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.inter(segment))

    def contain(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.contain(segment))

    def __iter__(self):
        # type: () -> Iterator[Read]
        return self.reads.values().__iter__()

    def __getitem__(self, read_id):
        # type: (basestring) -> Read
        return self.reads[read_id.split()[0]]

    def __len__(self):
        return len(self.reads)

    def print_fasta(self, hander):
        # type: (file) -> None
        for read in self:
            SeqIO.write(read, hander, "fasta")

    def print_alignments(self, handler):
        # type: (file) -> None
        for read in self:
            handler.write(read.id + "\n")
            for rec in read.alignments:
                handler.write(str(rec) + "\n")

    def asSeqRecords(self):
        # type: () -> Generator[SeqIO.SeqRecord]
        for read in self.reads.values():
            yield SeqIO.SeqRecord(read.seq, read.id)

    def loadFromFasta(self, handler):
        # type: (file) -> ReadCollection
        for rec in SeqIO.parse_fasta(handler):
            self.add(Read(rec))
        return self