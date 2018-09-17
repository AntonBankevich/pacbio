import sys

from dag_resolve import params

sys.path.append("py")
from common import sam_parser, SeqIO, basic
from typing import Generator, Iterator, Dict, Tuple, Optional, Union


class ContigCollection():
    def __init__(self, contigs_list = None):
        # type: (Optional[list[Contig]]) -> ContigCollection
        self.contigs = dict() # type: Dict[str, Contig]
        if contigs_list is not None:
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

class TmpInfo:
    def __init__(self, l):
        self.misc = l

class Contig:
    def __init__(self, seq, id, info = None):
        # type: (str, Union[str,int], TmpInfo) -> Contig
        if info is None:
            info = []
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

    def asSegment(self):
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
        # type: (Union[Contig, AlignedRead], int, int) -> Segment
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

    def __len__(self):
        return self.right - self.left + 1


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

class AlignedRead:
    def __init__(self, rec):
        # type: (SeqIO.SeqRecord) -> AlignedRead
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

    def sort(self):
        self.alignments = sorted(self.alignments, key = lambda alignment: alignment.seg_from.left)

    def __str__(self):
        self.sort()
        return str(self.id) + "(" + str(len(self.seq)) + ")" + "[" + ".".join(map(str, self.alignments)) + "]"

    def setSeq(self, seq):
        assert self.seq == "", "Changing assigned read sequence"
        self.seq = seq
        for al in self.alignments:
            if al.rc:
                al.seg_from.left += len(seq)
            else:
                al.seg_from.right += len(seq)

    def removeContig(self, contig):
        self.alignments = filter(lambda alignment: alignment.seg_to.contig.id != contig.id, self.alignments)

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
        # type: (ReadCollection) -> ReadCollection
        for read in other_collection.reads.values():
            self.add(read)
        return self

    def addNewRead(self, rec):
        # type: (SeqIO.SeqRecord) -> AlignedRead
        rec.id = rec.id.split()[0]
        read = AlignedRead(rec)
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

    def loadFromSam(self, sam, filter = lambda rec: True):
        # type: (sam_parser.Samfile) -> ReadCollection
        for rec in sam:
            if rec.is_unmapped or not filter(rec):
                continue
            if rec.query_name not in self.reads:
                new_read = AlignedRead(SeqIO.SeqRecord("", rec.query_name))
                self.add(new_read)
            else:
                new_read = self.reads[rec.query_name]
            if "H" not in rec.cigar and rec.seq != "*" and new_read.seq == "":
                if rec.rc:
                    new_read.setSeq(basic.RC(rec.seq))
                else:
                    new_read.setSeq(rec.seq)
            self.addNewAlignment(rec)
        for read in self.reads.values():
            assert read.seq != "", "Read " + read.id + " was aligned but had no primary alignments"
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


class Consensus:
    def __init__(self, seq, cov, full_seq = None):
        # type: (str, list[int], Optional[Consensus]) -> Consensus
        self.seq = seq
        self.cov = cov
        if full_seq is None:
            self.full = self # type: Consensus
        else:
            self.full = full_seq # type: Consensus

    def suffix(self, pos):
        if self == self.full:
            return Consensus(self.seq[pos:], self.cov[pos:])
        else:
            return Consensus(self.seq[pos:], self.cov[pos:], self.full.suffix(pos))

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
        l = min(10, len(self.seq))
        while l < len(self.seq) and self.cov[l] >= cov_threshold:
            l += 1
        if length is not None:
            l = min(l, length)
        return Consensus(self.seq[:l], self.cov[:l], self.full)

    def __len__(self):
        return len(self.seq)