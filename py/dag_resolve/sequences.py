import sys
sys.path.append("py")
from common import sam_parser, SeqIO
from typing import Generator


class ContigCollection():
    def __init__(self, contigs_list = []):
        self.contigs = dict()
        for contig in contigs_list:
            self.add(contig)
        self.main = None

    def add(self, contig):
        # type: (Contig) -> None
        if "main" in contig.info.misc:
            self.main = contig
        self.contigs[contig.id] = contig

    def filter(self, condition):
        # type: (function) -> ContigCollection
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

class TmpInfo:
    def __init__(self, l):
        self.misc = l

class Contig:
    def __init__(self, seq, id, info):
        # type: (basestring, int, TmpInfo) -> Contig
        self.seq = seq
        self.id = id
        if isinstance(info, list):
            info = TmpInfo(list)
        self.info = info

    def end(self, len):
        # type: (int) -> Segment
        return Segment(self, self.__len__() - len, self.__len__())

    def start(self, len):
        # type: (int) -> Segment
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
        return Contig(self.contig.seq[self.left:self.right], self.contig.id + "[" + str(self.left) + "," + str(self.right) + "]", self.contig.info)

    def __str__(self):
        # type: () -> basestring
        return self.contig.id + "[" + str(self.left) + ":" + str(self.right) + "]"


class AlignmentPiece:
    def __init__(self, seg_from, seg_to, cigar = None):
        # type: (Segment, Segment, basestring) -> AlignmentPiece
        self.seg_from = seg_from
        self.seg_to = seg_to
        self.cigar = cigar

    def __str__(self):
        # type: () -> basestring
        return "(" + str(self.seg_from) + "->" + str(self.seg_to) + ")"

class Read:
    def __init__(self, rec):
        # type: (SeqIO.SeqRecord) -> Read
        self.id = rec.id
        self.seq = rec.seq
        self.alignments = []

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def AddSamAlignment(self, rec, contig):
        # type: (sam_parser.SAM_entry, Contig) -> None
        cigar_list = list(sam_parser.CigarToList(rec.cigar))
        ls = 0
        rs = 0
        if cigar_list[0][0] in "HS":
            ls = cigar_list[0][1]
        if cigar_list[-1][0] in "HS":
            rs = cigar_list[-1][1]
        self.alignments.append(AlignmentPiece(Segment(self, ls, len(self.seq) - rs), Segment(contig, rec.pos, rec.pos + rec.alen), rec.cigar))

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
    def __init__(self, contigs):
        # type: (list[Contig]) -> ReadCollection
        self.reads = dict()
        self.contigs = contigs

    def addNewRead(self, rec):
        # type: (SeqIO.SeqRecord) -> Read
        rec.id = rec.id.split()[0]
        read = Read(rec)
        self.reads[rec.id] = read
        return read

    def addNewAlignment(self, rec):
        # type: (sam_parser.SAM_entry) -> None
        if rec.is_unmapped:
            return
        rname = rec.query_name.split()[0]
        if rname not in self.reads:
            return
        self.reads[rname].AddSamAlignment(rec, self.contigs[rec.tname])

    def add(self, read):
        # type: (Read) -> None
        self.reads[read.id] = read

    def filter(self, condition):
        # type: (function) -> ReadCollection
        res = ReadCollection(self.contigs)
        for read in self.reads.values():
            if condition(read):
                res.add(read)
        return res

    def inter(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.inter(segment))

    def contain(self, segment):
        # type: (Segment) -> ReadCollection
        return self.filter(lambda read: read.contain(segment))

    def __iter__(self):
        # type: () -> Generator[Read]
        return self.reads.values().__iter__()

    def __getitem__(self, read_id):
        # type: (int) -> Read
        return self.reads[read_id]

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
