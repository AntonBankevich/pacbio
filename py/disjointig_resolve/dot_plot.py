import itertools

from typing import Dict, List, Any

from alignment.align_tools import Aligner
from common import basic
from common.alignment_storage import AlignmentPiece, Correction
from disjointig_resolve.accurate_line import NewLineStorage, NewLine
from disjointig_resolve.smart_storage import LineListener, AlignmentStorage
from common.save_load import TokenWriter, TokenReader
from common.sequences import Segment, Contig, UniqueList


class LineDotPlot(LineListener):
    def __init__(self, lines):
        # type: (NewLineStorage) -> None
        LineListener.__init__(self, self)
        self.lines = lines
        for line in lines:
            line.addListener(self)
        self.auto_alignments = dict() # type: Dict[str, AutoAlignmentStorage] # here we store alignments of line to itself
        self.rc_alignments = dict() # type: Dict[str, RCAlignmentStorage] # here we stora alignments of line to its RC
        self.alignmentsTo = dict() # type: Dict[str, Dict[str, TwoLineAlignmentStorage]] # here we store all the remaining alignments
        for line in lines:
            self.alignmentsTo[line.id] = dict()
        for line in UniqueList(lines):
            self.auto_alignments[line.id] = AutoAlignmentStorage(line)
            self.auto_alignments[line.rc.id] = self.auto_alignments[line.id].rc
            self.rc_alignments[line.id] = RCAlignmentStorage(line)
            self.rc_alignments[line.rc.id] = self.rc_alignments[line.id].rc

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.simpleAddAlignment(al)
        self.simpleAddAlignment(al.reverse())

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
            if from_id not in self.alignmentsTo[from_id]:
                self.addTwoLineStorage(from_line, to_line)
            self.alignmentsTo[to_id][from_id].add(al)

    def addTwoLineStorage(self, line1, line2):
        # type: (NewLine, NewLine) -> TwoLineAlignmentStorage
        storage = TwoLineAlignmentStorage(line1, line2)
        if line2.id not in self.alignmentsTo:
            self.alignmentsTo[line2.id] = dict()
        self.alignmentsTo[line2.id][line1.id] = storage
        self.alignmentsTo[line2.rc.id][line1.rc.id] = storage.rc
        self.alignmentsTo[line1.id][line2.id] = storage.reverse
        self.alignmentsTo[line1.rc.id][line2.rc.id] = storage.rc.reverse
        return storage

    def addRCLineStorage(self, line):
        # type: (NewLine) -> RCAlignmentStorage
        storage = RCAlignmentStorage(line)
        self.rc_alignments[line.id] = storage
        self.rc_alignments[line.rc.id] = storage.rc
        return storage

    def addSelfAlignmentStorage(self, line):
        # type: (NewLine) -> AutoAlignmentStorage
        storage = AutoAlignmentStorage(line)
        self.auto_alignments[line.id] = storage
        self.auto_alignments[line.rc.id] = storage.rc
        return storage

    def getAllAlignments(self, seg):
        # type: (Segment) -> List[AlignmentPiece]
        # IMPLEMENT
        pass

    def construct(self, aligner):
        # type: (Aligner) -> None
        # IMPLEMENT
        pass

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        for d in self.alignmentsTo[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeExtendRight(line, new_seq, seq)
        self.auto_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)
        self.rc_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)


    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        for d in self.alignmentsTo[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeExtendRight(line, new_seq, seq)
        self.auto_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)
        self.rc_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        line = alignments.seq_to # type: NewLine
        for d in self.alignmentsTo[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireBeforeCorrect(alignments)
        self.auto_alignments[line.id].fireBeforeCorrect(alignments)
        self.rc_alignments[line.id].fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
        for d in self.alignmentsTo[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterExtendRight(line, seq)
        self.auto_alignments[line.id].fireAfterExtendRight(line, seq)
        self.rc_alignments[line.id].fireAfterExtendRight(line, seq)

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        for d in self.alignmentsTo[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterCutRight(line, pos)
        self.auto_alignments[line.id].fireAfterCutRight(line, pos)
        self.rc_alignments[line.id].fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        for d in self.alignmentsTo[line.id]: # type: Dict[str, TwoLineAlignmentStorage]
            for storage in d.values():
                storage.fireAfterCorrect(line)
        self.auto_alignments[line.id].fireAfterCorrect(line)
        self.rc_alignments[line.id].fireAfterCorrect(line)

    def save(self, handler):
        # type: (TokenWriter) -> None
        for l1, d1 in self.alignmentsTo.items():
            if not basic.isCanonocal(l1):
                continue
            for l2, als in d1.items():
                if l1 < basic.Normalize(l2):
                    handler.writeTokens([l1, l2])
                    als.save(handler)
        handler.writeTokens(["0", "0"])
        for lid in sorted(self.lines.lines.keys()):
            storage = self.rc_alignments[lid]
            storage.save(handler)
        for lid in sorted(self.lines.lines.keys()):
            storage = self.auto_alignments[lid]
            storage.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        while True:
            l1 = handler.readToken()
            l2 = handler.readToken()
            if l1 == "0" and l2 == "0":
                break
            storage = self.addTwoLineStorage(self.lines[l1], self.lines[l2])
            storage.load(handler, self.lines)
        for lid in sorted(self.lines.lines.keys()):
            storage = self.rc_alignments[lid]
            storage.load(handler, self.lines)
        for lid in sorted(self.lines.lines.keys()):
            storage = self.auto_alignments[lid]
            storage.load(handler, self.lines)


class TwoLineAlignmentStorage(LineListener):
    def __init__(self, line_from, line_to, rc = None, reverse = None):
        # type: (NewLine, NewLine, TwoLineAlignmentStorage, TwoLineAlignmentStorage) -> None
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

    def normalizeReverse(self):
        self.reverse.content.clean()
        for al in self.content:
            self.reverse.content.add(al.reverse())

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
        # type: (TokenReader, NewLineStorage) -> None
        self.content.load(handler, lines, lines)
        self.normalizeReverse()
