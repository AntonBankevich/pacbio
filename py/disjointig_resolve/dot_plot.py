from typing import Dict, List, Any, Optional, Generator

from alignment.align_tools import Aligner
from common import basic, params
from common.alignment_storage import AlignmentPiece, Correction
from disjointig_resolve.accurate_line import NewLineStorage, NewLine, LineStorageListener
from disjointig_resolve.smart_storage import LineListener, AlignmentStorage
from common.save_load import TokenWriter, TokenReader
from common.sequences import Segment, UniqueList, Contig, ContigStorage


# Unlike all other storages AutoAlignmentStorage stores only half of the alignments. The other half are the reversals of the stored half.
# Also identical alignment can not be stored here but is returned as a part of iteration (see __iter__)
class AutoAlignmentStorage(LineListener):

    def __init__(self, line, rc = None):
        # type: (Contig, Optional[AutoAlignmentStorage]) -> None
        self.line = line
        if rc is None:
            self.content = AlignmentStorage()
            rc = AutoAlignmentStorage(line.rc, self)
            rc.content = self.content.rc
        LineListener.__init__(self, rc)
        self.rc = rc # type: AutoAlignmentStorage
        self.state = 1 # from precedes to

    def makeCanonical(self, al):
        if (self.state == 1) == (al.seg_from.left < al.seg_from.right):
            return al
        else:
            return al.reverse()

    def isCanonical(self, al):
        return (self.state == 1) == (al.seg_from.left < al.seg_from.right)

    def add(self, al):
        # type: (AlignmentPiece) -> None
        if al.isIdentical():
            return
        self.content.add(self.makeCanonical(al))

    def addAndMergeRight(self, al):
        if al.isIdentical():
            return
        if self.isCanonical(al):
            self.content.addAndMergeRight(al)
        else:
            self.content.addAndMergeLeft(al.reverse())

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        for al in self.content:
            yield al
        for al in self.content:
            yield al.reverse()
        yield AlignmentPiece.Identical(self.line.asSegment())

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        for al in self:
            if al.seg_to.contains(seg):
                yield al

    def allInter(self, seg):
        for al in self:
            if al.seg_to.inter(seg):
                yield al

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        self.content.fireBeforeExtendRight(line, new_seq, seq)
        self.reverse()
        self.content.fireBeforeExtendRight(line, new_seq, seq)

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        self.content.fireBeforeCutRight(line, new_seq, pos)
        self.reverse()
        self.content.fireBeforeCutRight(line, new_seq, pos)
    # alignments from new sequence to new sequence

    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.content.fireBeforeCorrect(alignments)
        self.reverse()
        self.content.fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.content.fireAfterExtendRight(line, seq)
        self.reverse()
        self.content.fireAfterExtendRight(line, seq)

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.content.fireAfterCutRight(line, pos)
        self.reverse()
        self.content.fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        self.content.fireAfterCorrect(line)
        self.reverse()
        self.content.fireAfterCorrect(line)

    def reverse(self):
        self.state = -self.state
        self.content = self.content.reverse()
        self.rc.content = self.content.rc

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        self.content.load(handler, self.line, self.line)


class RCAlignmentStorage(LineListener):

    def __init__(self, line, rc = None):
        # type: (Contig, Optional[RCAlignmentStorage]) -> None
        self.line = line
        if rc is None:
            self.content = AlignmentStorage()
            rc = RCAlignmentStorage(line.rc, self)
            rc.content = self.content.rc
        LineListener.__init__(self, rc)
        self.rc = rc # type: AutoAlignmentStorage

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        return self.content.__iter__()

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        return self.content.getAlignmentsTo(seg)

    def allInter(self, seg):
        return self.content.allInter(seg)

    def add(self, alignment):
        self.content.add(alignment)
        self.content.add(alignment.reverse().rc)

    def addAndMergeRight(self, al):
        # type: (AlignmentPiece) -> None
        self.content.addAndMergeRight(al)
        self.content.addAndMergeLeft(al.reverse().rc)

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        self.content.fireBeforeExtendRight(line, new_seq, seq)
        self.reverse()
        self.content.fireBeforeExtendRight(line, new_seq, seq)

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        self.content.fireBeforeCutRight(line, new_seq, pos)
        self.reverse()
        self.content.fireBeforeCutRight(line, new_seq, pos)
    # alignments from new sequence to new sequence

    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.content.fireBeforeCorrect(alignments)
        self.reverse()
        self.content.fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.content.fireAfterExtendRight(line, seq)
        self.reverse()
        self.content.fireAfterExtendRight(line, seq)

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.content.fireAfterCutRight(line, pos)
        self.reverse()
        self.content.fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        self.content.fireAfterCorrect(line)
        self.reverse()
        self.content.fireAfterCorrect(line)
    # This is CRAAAZY!!! But correct.

    def reverse(self):
        self.rc.content = self.content.reverse()
        self.content = self.rc.content.rc

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        self.content.load(handler, self.line, self.line)


class TwoLineAlignmentStorage(LineListener):

    def __init__(self, line_from, line_to, rc = None, reverse = None):
        # type: (Contig, Contig, Optional[TwoLineAlignmentStorage], TwoLineAlignmentStorage) -> None
        assert line_from.id != line_to.id and line_from.rc.id != line_to.id
        self.line_from = line_from
        self.line_to = line_to
        if rc is None:
            self.content = AlignmentStorage()
            if reverse is None:
                rc = TwoLineAlignmentStorage(line_from.rc, line_to.rc, self, reverse)
            else:
                rc = TwoLineAlignmentStorage(line_from.rc, line_to.rc, self, reverse.rc)
        else:
            self.content = rc.content.rc # type: AlignmentStorage
        LineListener.__init__(self, rc)
        self.rc = rc # type: TwoLineAlignmentStorage
        if reverse is None:
            reverse = TwoLineAlignmentStorage(line_to, line_from, None, self)
        self.reverse = reverse
        self.rc.reverse = self.reverse.rc

    def add(self, al):
        # type: (AlignmentPiece) -> None
        self.content.add(al)
        reverse = al.reverse()
        self.reverse.content.add(reverse)

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        return self.content.__iter__()

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        return self.content.getAlignmentsTo(seg)

    def allInter(self, seg):
        return self.content.allInter(seg)

    def normalizeReverse(self):
        self.reverse.content = self.content.reverse()

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        self.content.fireBeforeExtendRight(line, new_seq, seq)
        self.normalizeReverse()

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        self.content.fireBeforeCutRight(line, new_seq, pos)
        self.normalizeReverse()
    # alignments from new sequence to new sequence

    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.content.fireBeforeCorrect(alignments)
        self.normalizeReverse()

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.content.fireAfterExtendRight(line, seq)
        self.normalizeReverse()

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.content.fireAfterCutRight(line, pos)
        self.normalizeReverse()

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        self.content.fireAfterCorrect(line)
        self.normalizeReverse()

    def addAndMergeRight(self, al):
        self.content.addAndMergeRight(al)
        self.normalizeReverse()

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler, lines):
        # type: (TokenReader, Any) -> None
        self.content.load(handler, lines, lines)
        self.normalizeReverse()


class DotPlot:
    def __init__(self, lines):
        # type: (ContigStorage) -> None
        self.lines = lines
        self.auto_alignments = dict() # type: Dict[str, AutoAlignmentStorage] # here we store alignments of line to itself
        self.alignmentsToFrom = dict() # type: Dict[str, Dict[str, TwoLineAlignmentStorage]] # here we store all the remaining alignments
        self.rc_alignments = dict() # type: Dict[str, RCAlignmentStorage] # here we stora alignments of line to its RC
        for line in lines.unique():
            self.addLine(line)

    def addLine(self, line):
        self.alignmentsToFrom[line.id] = dict()
        self.alignmentsToFrom[line.rc.id] = dict()
        self.addRCAlignmentStorage(line)
        self.addSelfAlignmentStorage(line)

    def addTwoLineStorage(self, line1, line2):
        # type: (Contig, Contig) -> TwoLineAlignmentStorage
        storage = TwoLineAlignmentStorage(line1, line2)
        if line2.id not in self.alignmentsToFrom:
            self.alignmentsToFrom[line2.id] = dict()
        self.alignmentsToFrom[line2.id][line1.id] = storage
        self.alignmentsToFrom[line2.rc.id][line1.rc.id] = storage.rc
        self.alignmentsToFrom[line1.id][line2.id] = storage.reverse
        self.alignmentsToFrom[line1.rc.id][line2.rc.id] = storage.rc.reverse
        return storage

    def simpleAddAlignment(self, al):
        # type: (AlignmentPiece) -> None
        to_line = al.seg_to.contig # type: NewLine
        from_line = al.seg_from.contig # type: NewLine
        to_id = al.seg_to.contig.id
        from_id = al.seg_from.contig.id
        if from_id == to_id:
            self.auto_alignments[from_id].add(al)
        elif to_line == from_line.rc:
            self.rc_alignments[to_id].add(al)
        else:
            if from_id not in self.alignmentsToFrom[from_id]:
                self.addTwoLineStorage(from_line, to_line)
            self.alignmentsToFrom[to_id][from_id].add(al)

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.simpleAddAlignment(al)
        self.simpleAddAlignment(al.reverse())

    def addRCAlignmentStorage(self, line):
        # type: (Contig) -> RCAlignmentStorage
        storage = RCAlignmentStorage(line)
        self.rc_alignments[line.id] = storage
        self.rc_alignments[line.rc.id] = storage.rc
        return storage

    def addSelfAlignmentStorage(self, line):
        # type: (Contig) -> AutoAlignmentStorage
        storage = AutoAlignmentStorage(line)
        self.auto_alignments[line.id] = storage
        self.auto_alignments[line.rc.id] = storage.rc
        return storage

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        for al in self.auto_alignments[seg.contig.id].getAlignmentsTo(seg):
            yield al
        for al in self.rc_alignments[seg.contig.id].getAlignmentsTo(seg):
            yield al
        for storage in self.alignmentsToFrom[seg.contig.id].values():
            for al in storage.getAlignmentsTo(seg):
                yield al

    def allInter(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        for al in self.auto_alignments[seg.contig.id].allInter(seg):
            yield al
        for al in self.rc_alignments[seg.contig.id].allInter(seg):
            yield al
        for storage in self.alignmentsToFrom[seg.contig.id].values():
            for al in storage.allInter(seg):
                yield al

    def construct(self, aligner):
        # type: (Aligner) -> None
        sequences = UniqueList(self.lines)
        for al in aligner.alignClean(sequences, ContigStorage(sequences)):
            if len(al) > params.k and al.percentIdentity() > 0.8:
                self.addAlignment(al)

    def save(self, handler):
        # type: (TokenWriter) -> None
        keys = [key for key in self.lines.items.keys() if basic.isCanonocal(key)]
        handler.writeTokens(keys)
        for l1, d1 in self.alignmentsToFrom.items():
            if not basic.isCanonocal(l1):
                continue
            for l2, als in d1.items():
                if l1 < basic.Normalize(l2):
                    handler.writeTokens([l1, l2])
                    als.save(handler)
        handler.writeTokens(["0", "0"])
        for lid in keys:
            storage = self.rc_alignments[lid]
            storage.save(handler)
        for lid in keys:
            storage = self.auto_alignments[lid]
            storage.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        keys = handler.readTokens()
        while True:
            l1 = handler.readToken()
            l2 = handler.readToken()
            if l1 == "0" and l2 == "0":
                break
            storage = self.addTwoLineStorage(self.lines[l1], self.lines[l2])
            storage.load(handler, self.lines)
        for lid in keys:
            storage = self.rc_alignments[lid]
            storage.load(handler)
        for lid in keys:
            storage = self.auto_alignments[lid]
            storage.load(handler)


class LineDotPlot(LineListener, LineStorageListener, DotPlot):

    def __init__(self, lines, aligner):
        # type: (NewLineStorage, Aligner) -> None
        DotPlot.__init__(self, lines)
        LineListener.__init__(self, self)
        self.lines = lines # type: NewLineStorage
        self.lines.addListener(self)
        for line in lines.unique():
            line.addListener(self)
        self.aligner = aligner

    def FireMergedLines(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        new_line = al1.seg_to.contig
        line1 = al1.seg_from.contig.id
        line2 = al2.seg_from.contig.id
        self.addLine(new_line)
        auto_alignments = self.auto_alignments[new_line.id]
        for al in self.auto_alignments[line1]:
            auto_alignments.add(al1.composeTargetDifference(al.compose(al1)))
        for al in self.auto_alignments[line2]:
            auto_alignments.add(al2.composeTargetDifference(al.compose(al2)))
        rc_alignments = self.rc_alignments[new_line.id]
        for al in self.rc_alignments[line1]:
            rc_alignments.add(al1.rc.composeTargetDifference(al.compose(al1)))
        for al in self.rc_alignments[line2]:
            rc_alignments.add(al2.rc.composeTargetDifference(al.compose(al2)))
        for al in self.getAlignmentsToFrom(line1, line2):
            auto_alignments.add(al2.composeTargetDifference(al.compose(al1)))
        for al in self.getAlignmentsToFrom(line1, line2.rc):
            rc_alignments.add(al2.rc.composeTargetDifference(al.compose(al1)))
        for storage in self.alignmentsToFrom[line1.id].values():
            for al in storage:
                self.addAlignment(al.compose(al1))
        for storage in self.alignmentsToFrom[line2.id].values():
            for al in storage:
                self.addAlignment(al.compose(al2))
        self.removeLine(al1.seg_from.contig)
        self.removeLine(al2.seg_from.contig)

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeExtendRight(line, new_seq, seq)
        self.auto_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)
        self.rc_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)

    def getAlignmentsToFrom(self, line_to, line_from):
        # type: (NewLine, NewLine) -> Generator[AlignmentPiece]
        if line_to == line_from:
            return self.auto_alignments[line_to.id].__iter__()
        elif line_to == line_from.rc:
            return self.rc_alignments[line_to.id].__iter__()
        elif line_from.id in self.alignmentsToFrom[line_to.id]:
            return self.alignmentsToFrom[line_to.id][line_from.id].__iter__()

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeCutRight(line, new_seq, pos)
        self.auto_alignments[line.id].fireBeforeCutRight(line, new_seq, pos)
        self.rc_alignments[line.id].fireBeforeCutRight(line, new_seq, pos)

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        line = alignments.seq_to # type: NewLine
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeCorrect(alignments)
        self.auto_alignments[line.id].fireBeforeCorrect(alignments)
        self.rc_alignments[line.id].fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (NewLine, str, Optional[List[AlignmentPiece]]) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterExtendRight(line, seq)
        self.auto_alignments[line.id].fireAfterExtendRight(line, seq)
        self.rc_alignments[line.id].fireAfterExtendRight(line, seq)
        new_seg = line.asSegment().suffix(length=len(seq) + 1000)
        for al in self.aligner.alignClean([new_seg.asContig()], self.lines):
            al = al.queryAsSegment(new_seg)
            self.addAndMergeRight(al)


    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterCutRight(line, pos)
        self.auto_alignments[line.id].fireAfterCutRight(line, pos)
        self.rc_alignments[line.id].fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterCorrect(line)
        self.auto_alignments[line.id].fireAfterCorrect(line)
        self.rc_alignments[line.id].fireAfterCorrect(line)

    def addAndMergeRight(self, al):
        # type: (AlignmentPiece) -> None
        if al.seg_to.contig == al.seg_from.contig:
            self.auto_alignments[al.seg_to.contig.id].addAndMergeRight(al)
        elif al.seg_to.contig == al.seg_from.contig.rc:
            self.rc_alignments[al.seg_to.contig.id].addAndMergeRight(al)
        else:
            self.alignmentsToFrom[al.seg_to.contig.id][al.seg_from.contig.id].addAndMergeRight(al)

    def removeLine(self, line):
        # type: (NewLine) -> None
        for storage in self.alignmentsToFrom[line.id].values():
            line_from = storage.line_from # type: NewLine
            del self.alignmentsToFrom[line_from.id][line.id]
        del self.alignmentsToFrom[line.id]
        del self.alignmentsToFrom[line.rc.id]
        self.deleteRCAlignmentStorage(line)
        self.deleteSelfAlignmentStorage(line)

    def deleteRCAlignmentStorage(self, line):
        # type: (NewLine) -> None
        del self.rc_alignments[line.id]
        del self.rc_alignments[line.rc.id]

    def deleteSelfAlignmentStorage(self, line):
        # type: (NewLine) -> None
        del self.auto_alignments[line.id]
        del self.auto_alignments[line.rc.id]

