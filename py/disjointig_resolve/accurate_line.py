import itertools

from typing import Optional, Iterable, List, Iterator, BinaryIO, Dict, Any, Generator

from alignment.align_tools import Aligner
from common import basic, SeqIO
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, UniqueList, ReadCollection, ContigCollection, Contig, Contig, ContigStorage
from common.alignment_storage import AlignmentPiece, AlignedRead, Correction
from disjointig_resolve.disjointigs import DisjointigCollection, UniqueMarker, Disjointig
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage, LineListener

class ReadAlignmentListener(LineListener):

    def __init__(self, line, rc = None):
        # type: (NewLine, Optional[ReadAlignmentListener]) -> None
        self.line = line # type: NewLine
        if rc is None:
            rc = ReadAlignmentListener(line.rc, self)
        LineListener.__init__(self, rc)
        self.rc = rc # type: ReadAlignmentListener

    def refreshReadAlignments(self):
        for al in self.line.read_alignments:
            read = al.seg_from.contig  # type: AlignedRead
            read.removeContig(self.line)
        for al in self.line.read_alignments:
            read = al.seg_from.contig  # type: AlignedRead
            read.addAlignment(al)

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        assert self.line == line, str((line.id, self.line.id))
        self.refreshReadAlignments()

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        assert self.line == line, str((line.id, self.line.id))
        self.refreshReadAlignments()

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.refreshReadAlignments()

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.refreshReadAlignments()

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.refreshReadAlignments()

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        self.refreshReadAlignments()

# This class supports alignment of reads and disjointigs to extended sequences
class ExtensionHandler(LineListener):
    def __init__(self, disjointigs, aligner):
        # type: (DisjointigCollection, Aligner) -> None
        LineListener.__init__(self, self)
        self.disjointigs = disjointigs
        self.aligner = aligner

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (NewLine, str, Optional[List[AlignmentPiece]]) -> None
        line = line # type: NewLine
        if relevant_als is not None:
            line.removeInter(line.suffix(1000 + len(seq))) # alignments of all reads to the end of line are manually removed
            for al in relevant_als:
                line.addReadAlignment(al.changeTargetContig(line)) # and then substituted with precalculated alignments to extendedd sequence
        new_seg = line.asSegment().suffix(length = len(seq) + 1000)
        for al in self.aligner.alignClean(new_seg.asContig(), self.disjointigs):
            al = al.reverse().targetAsSegment(new_seg)
            line.disjointig_alignments.addAndMergeRight(al)




# Reads only store alignments to lines they were selected.
class NewLine(Contig):
    def __init__(self, seq, id, extension_handler, rc = None):
        # type: (str, str, ExtensionHandler, Optional[NewLine]) -> None
        self.extensionHandler = extension_handler
        self.seq = seq
        self.id = id
        if rc is None:
            self.initial = AlignmentStorage()
            self.correct_segments = SegmentStorage()
            self.completely_resolved = SegmentStorage()
            self.disjointig_alignments = AlignmentStorage()
            self.read_alignments = AlignmentStorage()
            self.read_listener = ReadAlignmentListener(self)
            self.listeners = [self.initial, self.correct_segments, self.disjointig_alignments, self.read_alignments, self.read_listener] # type: List[LineListener]
            rc = NewLine(basic.RC(seq), basic.Reverse(self.id), extension_handler.rc, self) #type: NewLine
        else:
            self.initial = rc.initial.rc # type: AlignmentStorage
            self.correct_segments = rc.correct_segments.rc # type: SegmentStorage
            self.completely_resolved = rc.completely_resolved.rc # type: SegmentStorage
            self.disjointig_alignments = rc.disjointig_alignments.rc # type: AlignmentStorage
            self.read_alignments = rc.read_alignments.rc # type: AlignmentStorage
            self.read_listener = rc.read_listener.rc # type: ReadAlignmentListener
            self.listeners = [listener.rc for listener in rc.listeners] # type: List[LineListener]
        Contig.__init__(self, seq, id, rc)
        self.rc = rc #type: NewLine

    def addReads(self, alignments):
        # type: (Iterable[AlignmentPiece]) -> None
        self.read_alignments.addAll(alignments)

    def getReadAlignmentsTo(self, seg):
        # type: (Segment) -> Iterable[AlignmentPiece]
        return self.read_alignments.getAlignmentsTo(seg)

    def getPotentialAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        result = []
        for alDL in self.disjointig_alignments.getAlignmentsTo(seg):
            reduced = alDL.reduce(target=seg)
            dt = alDL.seg_from.contig # type: Disjointig
            for alRD in dt.getAlignmentsTo(reduced.seg_from):
                result.append(alRD.compose(alDL))
        result = sorted(result, key = lambda al: (al.seg_from.contig.id, -len(al.seg_from)))
        for read, iter in itertools.groupby(result, key = lambda al: al.seg_from.contig):
            readRes = []
            for al in iter:
                found = False
                for al1 in readRes:
                    inter = al.matchingSequence(True).inter(al1.matchingSequence(True))
                    if len(inter.matches) != 0:
                        found = True
                if not found:
                    yield al
                    readRes.append(al)

    def position(self, pos):
        # type: (int) -> LinePosition
        return LinePosition(self, pos)

    def extendRight(self, seq, relevant_als = None):
        # type: (str, List[AlignmentPiece]) -> None
        if relevant_als is None:
            relevant_als = []
        new_seq = Contig(self.seq + seq, "TMP_" + self.id)
        self.notifyBeforeExtendRight(new_seq, seq)
        self.seq = self.seq + seq
        self.rc.seq = basic.RC(seq) + self.rc.seq
        self.notifyAfterExtendRight(seq, relevant_als)
        # TODO: put the handling of read and disjointig alignments to extended sequence into listeners one way or another

    def notifyBeforeExtendRight(self, new_seq, seq):
        # type: (Contig, str) -> None
        for listener in self.listeners:
            listener.fireBeforeExtendRight(self, new_seq, seq)

    def notifyAfterExtendRight(self, seq, relevant_als):
        # type: (str, Optional[List[AlignmentPiece]]) -> None
        for listener in self.listeners:
            listener.fireAfterExtendRight(self, seq, relevant_als)


    def cutRight(self, pos):
        assert pos > 0 and pos <= len(self)
        cut_length = len(self) - pos
        if cut_length == 0:
            return
        new_seq = Contig(self.seq[:pos], "TMP_" + self.id)
        self.notifyBeforeCutRight(new_seq, pos)
        self.seq = self.seq[:-cut_length]
        self.rc.seq = self.rc.seq[cut_length:]
        self.notifyAfterCutRight(pos)

    def notifyBeforeCutRight(self, new_seq, pos):
        # type: (Contig, int) -> None
        for listener in self.listeners:
            listener.fireBeforeCutRight(self, new_seq, pos)

    def notifyAfterCutRight(self, pos):
        # type: (int) -> None
        for listener in self.listeners:
            listener.fireAfterCutRight(self, pos)

    def correctSequence(self, alignments):
        # type: (List[AlignmentPiece]) -> None
        assert len(alignments) > 0
        correction = Correction.constructCorrection(alignments)
        self.notifyBeforeCorrect(correction)
        self.seq = correction.seq_from.seq
        self.rc.seq = basic.RC(self.seq)
        self.notifyAfterCorrect()

    def notifyBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        for listener in self.listeners:
            listener.fireBeforeCorrect(alignments)

    def notifyAfterCorrect(self):
        # type: () -> None
        for listener in self.listeners:
            listener.fireAfterCorrect(self)


    # def extendRightWithAlignment(self, alignment):
    #     # type: (AlignmentPiece) -> None
    #     assert alignment.seg_from.contig == self
    #     self.cutRight(alignment.seg_from.right)
    #     self.correctSequence([alignment])
    #     self.extendRight(alignment.seg_to.contig.suffix(alignment.seg_to.right))

    def addReadAlignment(self, al):
        # type: (AlignmentPiece) -> AlignmentPiece
        self.read_alignments.add(al)
        read = al.seg_from.contig # type: AlignedRead
        read.addAlignment(al)
        return al

    def addListener(self, listener):
        self.listeners.append(listener)
        self.rc.addListener(listener.rc)

    def removeListener(self, listener):
        self.listeners.remove(listener)
        self.rc.listeners.remove(listener.rc)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)
        self.initial.save(handler)
        self.correct_segments.save(handler)
        self.completely_resolved.save(handler)
        self.disjointig_alignments.save(handler)
        self.read_alignments.save(handler)

    def loadLine(self, handler, disjointigs, reads, contigs):
        # type: (TokenReader, DisjointigCollection, ReadCollection, ContigCollection) -> None
        self.id = handler.readToken()
        self.rc.id = basic.RC(self.id)
        self.initial.load(handler, disjointigs, self)
        self.correct_segments.load(handler, self)
        self.completely_resolved.load(handler, self)
        self.disjointig_alignments.load(handler, disjointigs, self)
        self.read_alignments.load(handler, reads, self)
        for al in self.read_alignments:
            read = al.seg_from.contig #type: AlignedRead
            read.addAlignment(al)

    def __str__(self):
        points = [self.left()]
        points.extend(self.initial)
        points.append(self.right())
        points = map(str, points)
        return "Line:" + str(self.id) + ":" + "[" + ":".join(points) +"]"

    def removeInter(self, seg):
        als = self.read_alignments.allInter(seg)
        for al in als:
            read = al.seg_from.contig # type: AlignedRead
            read.alignments.remove(al)
        self.read_alignments.removeInter(seg)


class NewLineStorage(ContigStorage):
    def __init__(self, disjointigs, aligner):
        # type: (DisjointigCollection, Aligner) -> None
        ContigStorage.__init__(self, [], False)
        self.disjointigs = disjointigs
        self.aligner = aligner
        self.items = dict() # type: Dict[str, NewLine]
        self.cnt = 1

    def __iter__(self):
        # type: () -> Iterator[NewLine]
        return self.items.values().__iter__()

    def __getitem__(self, item):
        # type: (str) -> NewLine
        return self.items[item]

    def addNew(self, seq, name = None):
        # type: (str, Optional[str]) -> NewLine
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        new_line = NewLine(seq, name, ExtensionHandler(self.disjointigs, self.aligner))
        self.items[name] = new_line
        self.items[new_line.rc.id] = new_line.rc
        return new_line

    def fillFromContigs(self, contigs):
        # type: (Iterable[Contig]) -> None
        for contig in UniqueList(contigs):
            marker = UniqueMarker(self.disjointigs)
            segments = list(marker.findUnique(contig))
            if len(segments) == 0:
                print "No unique segments found in contig", contig
                continue
            line = self.addNew(contig.seq)
            line.initial.add(AlignmentPiece.Identical(line.asSegment(), contig.asSegment()))
            for seg in segments:
                line.correct_segments.add(seg.contigAsSegment(line.asSegment()))
            for seg in line.correct_segments:
                line.completely_resolved.add(seg)

    def fillFromDisjointigs(self):
        # type: () -> None
        for seg in UniqueMarker(self.disjointigs).findAllUnique(self.disjointigs):
            line = self.addNew(seg.Seq())
            line.initial.add(AlignmentPiece.Identical(line.asSegment(), seg))
        #TODO Filter all lines already present in the collection

    def mergeLines(self, alignment):
        # type: (AlignmentPiece) -> NewLine
        line1 = alignment.seg_from.contig #type: NewLine
        line2 = alignment.seg_to.contig #type: NewLine
        line = self.addNew(line1.seq, name ="(" + line1.id + "," + line2.id + ")")
        line.extendRightWithAlignment(alignment.changeQueryContig(line))
        # IMPLEMENT merge lines into a new line based on the given alignment. Polysh the result. Transfer all reads, segments and alignments
        del self.items[line1.id]
        del self.items[line2.id]
        return line

    def printToFile(self, handler):
        # type: (BinaryIO) -> None
        for line in self.items:
            handler.write(line.__str__() + "\n")

    def printToFasta(self, handler):
        # type: (BinaryIO) -> None
        for line in UniqueList(self.items.values()):
            SeqIO.write(line, handler, "fasta")

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(str(self.cnt))
        line_ids = map(lambda line: line.id, UniqueList(self.items.values()))
        handler.writeTokens(line_ids)
        for line_id in line_ids:
            line = self.items[line_id]
            handler.writeTokenLine(line.seq)
        for line_id in line_ids:
            line = self.items[line_id]
            line.save(handler)

    def load(self, handler, reads, contigs):
        # type: (TokenReader, ReadCollection, ContigCollection) -> None
        self.cnt = int(handler.readToken())
        keys = handler.readTokens()
        for key in keys:
            self.addNew(handler.readToken(), key)
        for key in keys:
            line = self.items[key]
            line.loadLine(handler, self.disjointigs, reads, contigs)

class LinePosition(LineListener):
    def __init__(self, line, pos, rc=None):
        # type: (NewLine, int, Optional[LinePosition]) -> None
        self.line = line
        self.pos = pos
        line.addListener(self)
        if rc is None:
            rc = LinePosition(line.rc, len(line) - 1 - pos)
        LineListener.__init__(self, rc)

    def suffix(self):
        return self.line.suffix(self.pos)

    def prefix(self):
        return self.line.prefix(self.pos)

    def fixRC(self):
        self.rc.pos = len(self.line) - 1 - self.pos

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        pass

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (NewLine, str, Optional[List[AlignmentPiece]]) -> None
        self.fixRC()

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        if pos >= line.right():
            self.pos = line.right() - 1
            self.fixRC()

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.fixRC()

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.pos = alignments.mapPositionsUp([self.pos])[0] # type: int

