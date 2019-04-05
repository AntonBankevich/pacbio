import itertools

from typing import Dict, List, Any, Optional, Iterable, Generator

from alignment.align_tools import Aligner
from common import basic, params
from common.alignment_storage import AlignmentPiece, Correction
from common.seq_records import NamedSequence
from disjointig_resolve.accurate_line import NewLineStorage, NewLine
from disjointig_resolve.smart_storage import LineListener, AlignmentStorage
from common.save_load import TokenWriter, TokenReader
from common.sequences import Segment, Contig, UniqueList, Contig


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

    def add(self, alignment):
        # type: (AlignmentPiece) -> None
        self.content.add(alignment)

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

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
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

    def add(self, alignment):
        self.content.add(alignment)
        self.content.add(alignment.reverse().rc)

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

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
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

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
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

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler, lines):
        # type: (TokenReader, Any) -> None
        self.content.load(handler, lines, lines)
        self.normalizeReverse()


class DotPlot:
    def __init__(self, lines):
        # type: (Iterable[Contig]) -> None
        self.lines = dict() # type: Dict[str, Contig]
        for line in lines:
            self.lines[line.id] = line
        for line in lines:
            self.alignmentsToFrom[line.id] = dict()
        for line in UniqueList(lines):
            self.addRCAlignmentStorage(line)
            self.addSelfAlignmentStorage(line)
        self.auto_alignments = dict() # type: Dict[str, AutoAlignmentStorage] # here we store alignments of line to itself
        self.alignmentsToFrom = dict() # type: Dict[str, Dict[str, TwoLineAlignmentStorage]] # here we store all the remaining alignments
        self.rc_alignments = dict() # type: Dict[str, RCAlignmentStorage] # here we stora alignments of line to its RC

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

    def construct(self, aligner):
        # type: (Aligner) -> None
        sequences = UniqueList(self.lines.values())
        for al in aligner.alignClean(sequences, sequences):
            if len(al) > params.k and al.percentIdentity() > 0.8:
                self.addAlignment(al)

    def save(self, handler):
        # type: (TokenWriter) -> None
        keys = [key for key in self.lines.keys() if basic.isCanonocal(key)]
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


class LineDotPlot(LineListener, DotPlot):

    def __init__(self, lines):
        # type: (Iterable[NewLine]) -> None
        DotPlot.__init__(self, lines)
        LineListener.__init__(self, self)
        for line in lines:
            line.addListener(self)

    # IMPLEMENT construct additional alignments when expanding maybe when correcting too. Maybe additional operation find new alignments of a segment
    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeExtendRight(line, new_seq, seq)
        self.auto_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)
        self.rc_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)


    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeExtendRight(line, new_seq, seq)
        self.auto_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)
        self.rc_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        line = alignments.seq_to # type: NewLine
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeCorrect(alignments)
        self.auto_alignments[line.id].fireBeforeCorrect(alignments)
        self.rc_alignments[line.id].fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
        for d in self.alignmentsToFrom[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterExtendRight(line, seq)
        self.auto_alignments[line.id].fireAfterExtendRight(line, seq)
        self.rc_alignments[line.id].fireAfterExtendRight(line, seq)

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

