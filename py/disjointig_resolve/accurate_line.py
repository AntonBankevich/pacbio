from typing import Optional, Iterable, List, Iterator, BinaryIO
from common import basic
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, AlignmentPiece, UniqueList, ReadCollection, ContigCollection, AlignedRead
from disjointig_resolve.disjointigs import DisjointigCollection, UniqueMarker


class NewLine(NamedSequence):
    def __init__(self, seq, id, rc = None, zero_pos = 0):
        # type: (str, str, Optional[NewLine], int) -> None
        NamedSequence.__init__(self, seq, id, zero_pos)
        if rc is None:
            self.rc = NewLine(basic.RC(seq), basic.Reverse(self.id), len(seq) - zero_pos) #type: NewLine
        else:
            self.rc = rc #type: NewLine
        self.initial = [] #type:List[AlignmentPiece]
        self.correct_segments = []
        self.read_alignments = [] #type: List[AlignmentPiece]

    def addReads(self, alignments):
        # type: (Iterable[AlignmentPiece]) -> None
        self.read_alignments.extend(alignments)
        self.rc.read_alignments.extend([al.RC() for al in alignments])
        self.read_alignments = sorted(self.read_alignments, key = lambda al: al.seg_to.left)
        self.rc.read_alignments = sorted(self.rc.read_alignments, key = lambda al: al.seg_to.left)

    def left(self):
        return -self.zero_pos

    def right(self):
        return len(self) - self.zero_pos

    def segment(self, start, end):
        return Segment(self, start, end)

    def asSegment(self):
        return self.segment(-self.zero_pos, len(self) - self.zero_pos)

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def suffix(self, pos = None, length = None):
        # type: (Optional[int], Optional[int]) -> Segment
        assert (pos is None) != (length is None), str([pos, length])
        if length is not None:
            pos = len(self) - length - self.zero_pos
        return self.segment(pos, len(self) - self.zero_pos)

    def prefix(self, pos = None, length = None):
        # type: (Optional[int], Optional[int]) -> Segment
        assert (pos is None) != (length is None), str([pos, length])
        if length is not None:
            pos = length - self.zero_pos
        return self.segment(pos, len(self) - self.zero_pos)

    def extendRight(self, seq):
        # type: (str) -> None
        self.seq = self.seq + seq
        self.rc.seq = basic.RC(seq) + self.rc.seq
        self.rc.zero_pos += len(seq)

    # def extendLeft(self, seq):
    #     # type: (str) -> None
    #     self.seq = seq + self.seq
    #     self.rc.seq = self.rc.seq + basic.RC(seq)
    #     self.zero_pos += len(seq)

    def cutRight(self, pos):
        assert pos > 0 and pos + self.zero_pos <= len(self)
        cut_length = len(self) - pos - self.zero_pos
        if cut_length == 0:
            return
        self.seq = self.seq[:-cut_length]
        self.rc.seq = self.rc.seq[cut_length:]
        self.rc.zero_pos -= cut_length

    # def cutLeft(self, pos):
    #     assert pos < 0 and 0 <= pos + self.zero_pos <= len(self)
    #     cut_length = pos + self.zero_pos
    #     if cut_length == 0:
    #         return
    #     self.seq = self.seq[cut_length:]
    #     self.rc.seq = self.rc.seq[:-cut_length]
    #     self.zero_pos -= cut_length

    def correctSequence(self, alignment):
        # type: (AlignmentPiece) -> None
        #IMPLEMENT fix read alignments, correct segments, initial alignments, sequence and zero_pos
        pass

    def extendRightWithAlignment(self, alignment):
        # type: (AlignmentPiece) -> None
        assert alignment.seg_from.contig == self
        self.cutRight(alignment.seg_from.right)
        self.correctSequence(alignment)
        self.extendRight(alignment.seg_to.contig.suffix(alignment.seg_to.right))

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)
        handler.writeIntLine(len(self.initial))
        for al in self.initial:
            al.save(handler)
        handler.writeIntLine(len(self.correct_segments))
        for seg in self.correct_segments:
            seg.save(handler)
        handler.writeIntLine(len(self.read_alignments))
        for al in self.read_alignments :
            al.save(handler)

    def load(self, handler, disjointigs, reads, contigs):
        # type: (TokenReader, DisjointigCollection, ReadCollection, ContigCollection) -> None
        self.id = handler.readToken()
        self.rc.id = basic.RC(self.id)
        n = handler.readInt()
        for al in range(n):
            self.initial.append(AlignmentPiece.load(handler, disjointigs, self))
        n = handler.readInt()
        for i in range(n):
            self.correct_segments.append(Segment.load(handler, self))
        n = handler.readInt()
        for i in range(n):
            al = AlignmentPiece.load(handler, reads, self)
            self.read_alignments.append(al)
            read = al.seg_from.contig #type: AlignedRead
            read.addAlignment(al)

    def __str__(self):
        points = [self.left()]
        points.extend(self.initial)
        points.append(self.right())
        points = map(str, points)
        return "Line:" + str(self.id) + ":" + "[" + ":".join(points) +"]"


class NewLineStorage:
    def __init__(self, disjointigs):
        # type: (DisjointigCollection) -> None
        self.disjointigs = disjointigs
        self.lines = dict()
        self.cnt = 1

    def __iter__(self):
        # type: () -> Iterator[NewLine]
        return self.lines.values().__iter__()

    def __contains__(self, item):
        # type: (NamedSequence) -> bool
        return item.id in self.lines

    def __getitem__(self, item):
        # type: (str) -> NewLine
        return self.lines[item]

    def containsKey(self, key):
        # type: (str) -> bool
        return key in self.lines

    def add(self, seq, name = None):
        # type: (str, Optional[str]) -> NewLine
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        new_line = NewLine(seq, name)
        self.lines[name] = new_line
        self.lines[new_line.rc.id] = new_line.rc
        return new_line

    def fillFromContigs(self, contigs):
        for contig in UniqueList(contigs):
            marker = UniqueMarker(self.disjointigs)
            segments = list(marker.findUnique(contig))
            if len(segments) == 0:
                print "No unique segments found in contig", contig
                continue
            line = self.add(contig.seq)
            line.initial.append(AlignmentPiece(line.asSegment(), contig.asSegment(), str(len(line)) + "M"))
            for seg in segments:
                line.correct_segments.append(seg.contigAsSegment(line.asSegment()))

    def fillFromDisjointigs(self):
        for seg in UniqueMarker(self.disjointigs).findAllUnique(self.disjointigs):
            line = self.add(seg.Seq())
            line.initial.append(AlignmentPiece(line.asSegment(), seg, str(len(line)) + "M"))
        #IMPLEMENT Filter all lines already present in the collection

    def mergeLines(self, alignment):
        # type: (AlignmentPiece) -> NewLine
        line1 = alignment.seg_from.contig #type: NewLine
        line2 = alignment.seg_to.contig #type: NewLine
        line = self.add(line1.seq, name = "(" + line1.id + "," + line2.id + ")")
        line.extendRightWithAlignment(alignment.changeQuery(line))
        # IMPLEMENT merge lines into a new line based on the given alignment. Polysh the result. Transfer all reads, segments and alignments
        del self.lines[line1.id]
        del self.lines[line2.id]
        return line

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(str(self.cnt))
        handler.writeTokens(map(lambda line: line.id, UniqueList(self.lines.values())))
        for line in UniqueList(self.lines.values()):
            line.save(handler)

    def load(self, handler, reads, contigs):
        # type: (TokenReader, ReadCollection, ContigCollection) -> None
        self.cnt = int(handler.readToken())
        keys = handler.readTokens()
        for key in keys:
            line = self.add("", key)
            line.load(handler, self.disjointigs, reads, contigs)
            assert key == line.id

    def printToFile(self, handler):
        # type: (BinaryIO) -> None
        for line in self.lines:
            handler.write(line.__str__() + "\n")

    def printToFasta(self, handler):
        # type: (BinaryIO) -> None
        # IMPLEMENT
        pass
