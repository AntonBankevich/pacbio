import itertools

from typing import Optional, Iterable, List, Any, Generator

from alignment.align_tools import Aligner
from common import basic, params
from common.save_load import TokenWriter, TokenReader
from common.sequences import Segment, ContigCollection, Contig
from common.alignment_storage import AlignmentPiece, AlignedRead, ReadCollection
from disjointig_resolve.correction import Correction
from disjointig_resolve.disjointigs import DisjointigCollection, Disjointig
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
            tmp = line.read_alignments.merge(AlignmentStorage().addAll(relevant_als).targetAsSegment(line.asSegment()))
            line.read_alignments.clean()
            line.read_alignments.addAll(tmp)
        new_seg = line.asSegment().suffix(length = len(seq) + 1000)
        for al in self.aligner.alignClean([new_seg.asContig()], self.disjointigs):
            al = al.reverse().targetAsSegment(new_seg)
            line.disjointig_alignments.addAndMergeRight(al)


# Reads only store alignments to lines they were selected.
class NewLine(Contig):
    def __init__(self, seq, id, extension_handler, rc = None):
        # type: (str, str, ExtensionHandler, Optional[NewLine]) -> None
        self.extensionHandler = extension_handler
        self.seq = seq
        self.id = id # type: str
        self.circular = False
        if rc is None:
            # TODO: move all these to separate classes
            self.initial = AlignmentStorage()
            self.correct_segments = SegmentStorage()
            self.completely_resolved = SegmentStorage()
            self.disjointig_alignments = AlignmentStorage()
            self.read_alignments = AlignmentStorage()
            self.listeners = [self.initial, self.correct_segments, self.completely_resolved, self.disjointig_alignments, self.read_alignments, extension_handler] # type: List[LineListener]
            rc = NewLine(basic.RC(seq), basic.Reverse(self.id), extension_handler.rc, self) #type: NewLine
            self.rc = rc
            self.addListener(ReadAlignmentListener(self))
        else:
            self.initial = rc.initial.rc # type: AlignmentStorage
            self.correct_segments = rc.correct_segments.rc # type: SegmentStorage
            self.completely_resolved = rc.completely_resolved.rc # type: SegmentStorage
            self.disjointig_alignments = rc.disjointig_alignments.rc # type: AlignmentStorage
            self.read_alignments = rc.read_alignments.rc # type: AlignmentStorage
            self.listeners = [listener.rc for listener in rc.listeners] # type: List[LineListener]
        Contig.__init__(self, seq, id, rc)
        self.rc = rc #type: NewLine

    def updateCorrectSegments(self, seg, threshold = params.reliable_coverage):
        # type: (Segment, int) -> None
        segs = AlignmentStorage().addAll(self.read_alignments.allInter(seg)).filterByCoverage(mi=threshold)
        self.correct_segments.addAll(segs)
        self.correct_segments.mergeSegments()

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

    def getRelevantAlignmentsFor(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        result = []
        for alDL in self.disjointig_alignments.allInter(seg):
            reduced = alDL.reduce(target=seg)
            dt = alDL.seg_from.contig # type: Disjointig
            for alRD in dt.allInter(reduced.seg_from):
                al = alRD.compose(alDL)
                if len(al.seg_to) >= params.k:
                    result.append(al)
        result = sorted(result, key = lambda al: (al.seg_from.contig.id, -len(al.seg_from)))
        for read, iter in itertools.groupby(result, key = lambda al: al.seg_from.contig): # type: AlignedRead, Generator[AlignmentPiece]
            readRes = []
            for al in iter:
                found = False
                for al1 in readRes:
                    inter = al.matchingSequence(True).inter(al1.matchingSequence(True))
                    if len(inter.matches) != 0:
                        found = True
                        break
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
        new_seq = Contig(self.seq + seq, "TMP2_" + self.id)
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
        new_seq = Contig(self.seq[:pos], "TMP3_" + self.id)
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
        # type: (Iterable[AlignmentPiece]) -> None
        alignments = list(alignments)
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
        self.rc.listeners.append(listener.rc)

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
        self.initial.load(handler, contigs, self)
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

    def setCircular(self):
        self.circular = True
        self.rc.circular = True

    def cleanReadAlignments(self):
        for read in self.read_alignments:
            read.seg_from.contig.removeContig(self)


class LinePosition(LineListener):
    def __init__(self, line, pos, rc=None):
        # type: (NewLine, int, Optional[LinePosition]) -> None
        self.line = line
        self.pos = pos
        if rc is None:
            rc = LinePosition(line.rc, len(line) - 1 - pos, self)
            self.rc = rc
            line.addListener(self)
        LineListener.__init__(self, rc)
        self.rc = rc # type: LinePosition

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

