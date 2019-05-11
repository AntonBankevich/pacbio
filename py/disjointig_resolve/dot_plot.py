from typing import Dict, List, Any, Optional, Generator, BinaryIO

from alignment.align_tools import Aligner
from common import basic, params
from common.alignment_storage import AlignmentPiece
from disjointig_resolve.correction import Correction
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_alignments import AutoAlignmentStorage, RCAlignmentStorage, TwoLineAlignmentStorage
from disjointig_resolve.line_storage import NewLineStorage, LineStorageListener
from disjointig_resolve.smart_storage import LineListener, AlignmentStorage
from common.save_load import TokenWriter, TokenReader
from common.sequences import Segment, Contig, ContigStorage


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

    def addAlignment(self, al):
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
            if from_id not in self.alignmentsToFrom[to_id]:
                self.addTwoLineStorage(from_line, to_line)
            self.alignmentsToFrom[to_id][from_id].add(al)

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
        for al in aligner.alignClean(self.lines.unique(), self.lines):
            # print al, len(al) > params.k, al.percentIdentity() > 0.8, al.seg_from.contig.id < al.seg_to.contig.id, al.seg_from <= al.seg_to
            if len(al) > params.k and al.percentIdentity() > 0.8 and (
                    al.seg_from.contig.id < al.seg_to.contig.id or al.seg_from <= al.seg_to):
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
                    handler.writeToken(l1)
                    handler.writeToken(l2)
                    handler.newLine()
                    als.save(handler)
        handler.writeToken("0")
        handler.writeToken("0")
        handler.newLine()
        for lid in keys:
            storage = self.rc_alignments[lid]
            storage.save(handler)
        for lid in keys:
            storage = self.auto_alignments[lid]
            storage.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        keys = list(handler.readTokens())
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
        LineListener.__init__(self, self)
        DotPlot.__init__(self, lines)
        self.lines = lines # type: NewLineStorage
        self.lines.addListener(self)
        # for line in lines.unique():
        #     line.addListener(self)
        self.aligner = aligner

    def addLine(self, line):
        # type: (NewLine) -> None
        DotPlot.addLine(self, line)
        line.addListener(self)

    def getAlignmentsToFrom(self, line_to, line_from):
        # type: (NewLine, NewLine) -> Generator[AlignmentPiece]
        content = []
        if line_to == line_from:
            content = self.auto_alignments[line_to.id]
        elif line_to == line_from.rc:
            content = self.rc_alignments[line_to.id]
        elif line_from.id in self.alignmentsToFrom[line_to.id]:
            content = self.alignmentsToFrom[line_to.id][line_from.id]
        for al in content:
            yield al

    def FireMergedLines(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        print "Fire merged lines", al1, al2
        new_line = al1.seg_to.contig
        line1 = al1.seg_from.contig
        line2 = al2.seg_from.contig
        print list(self.allInter(line1.asSegment()))
        print list(self.allInter(line2.asSegment()))
        self.addLine(new_line)
        self.auto_alignments[line1.id].setState(-1)
        self.auto_alignments[line2.id].setState(-1)
        auto1 = AutoAlignmentStorage(new_line).addAll([al1.composeTargetDifference(al.compose(al1)) for al in self.auto_alignments[line1.id].content])
        auto2 = AutoAlignmentStorage(new_line).addAll([al2.composeTargetDifference(al.compose(al2)) for al in self.auto_alignments[line2.id].content])
        auto3 = AutoAlignmentStorage(new_line).addAll([al2.composeTargetDifference(al.compose(al1)) for al in self.getAlignmentsToFrom(line1, line2)])
        self.auto_alignments[new_line.id].addAll(auto1.merge(auto3).merge(auto2).content)
        rc1 = RCAlignmentStorage(new_line).addAll([al1.rc.composeTargetDifference(al.compose(al1)) for al in self.rc_alignments[line1.id]])
        rc2 = RCAlignmentStorage(new_line).addAll([al2.rc.composeTargetDifference(al.compose(al2)) for al in self.rc_alignments[line2.id]])
        rc3 = RCAlignmentStorage(new_line).addAll([al2.rc.composeTargetDifference(al.compose(al1)) for al in self.getAlignmentsToFrom(line1, line2.rc)])
        self.rc_alignments[new_line.id].addAll(rc1.merge(rc3).merge(rc2))
        common = set(self.alignmentsToFrom[line1.id].keys()).intersection(set(self.alignmentsToFrom[line2.id].keys()))
        for storage in self.alignmentsToFrom[line1.id].values():
            if storage.line_from.id not in common and storage.line_from != line2:
                self.addTwoLineStorage(storage.line_from, new_line).addAll([al.compose(al1) for al in storage])
        for storage in self.alignmentsToFrom[line2.id].values():
            if storage.line_from.id not in common and storage.line_from != line1:
                self.addTwoLineStorage(storage.line_from, new_line).addAll([al.compose(al2) for al in storage])
        for c in common:
            storage1 = self.alignmentsToFrom[line1.id][c]
            storage2 = self.alignmentsToFrom[line2.id][c]
            als1 = AlignmentStorage().addAll([al.compose(al1) for al in storage1])
            als2 = AlignmentStorage().addAll([al.compose(al2) for al in storage2])
            self.addTwoLineStorage(storage1.line_from, new_line).addAll(als1.merge(als2))
        self.removeLine(al1.seg_from.contig)
        self.removeLine(al2.seg_from.contig)
        print list(self.allInter(new_line.asSegment()))

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        for storage in self.alignmentsToFrom[line.id].values():
            storage.fireBeforeExtendRight(line, new_seq, seq)
        self.auto_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)
        self.rc_alignments[line.id].fireBeforeExtendRight(line, new_seq, seq)

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        for storage in self.alignmentsToFrom[line.id].values():
            storage.fireBeforeCutRight(line, new_seq, pos)
        self.auto_alignments[line.id].fireBeforeCutRight(line, new_seq, pos)
        self.rc_alignments[line.id].fireBeforeCutRight(line, new_seq, pos)

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        line = alignments.seq_to # type: NewLine
        for storage in self.alignmentsToFrom[line.id].values():
            storage.fireBeforeCorrect(alignments)
        self.auto_alignments[line.id].fireBeforeCorrect(alignments)
        self.rc_alignments[line.id].fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (NewLine, str, Optional[List[AlignmentPiece]]) -> None
        for storage in self.alignmentsToFrom[line.id].values():
            storage.fireAfterExtendRight(line, seq)
        self.auto_alignments[line.id].fireAfterExtendRight(line, seq)
        self.rc_alignments[line.id].fireAfterExtendRight(line, seq)
        new_seg = line.asSegment().suffix(length=len(seq) + 1000)
        # print "Aligning new extension"
        for al in self.aligner.alignClean([new_seg.asContig()], self.lines):
            al = al.queryAsSegment(new_seg)
            self.addAndMergeRight(al)



    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        for storage in self.alignmentsToFrom[line.id].values():
            storage.fireAfterCutRight(line, pos)
        self.auto_alignments[line.id].fireAfterCutRight(line, pos)
        self.rc_alignments[line.id].fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        for storage in self.alignmentsToFrom[line.id].values():
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
        line.removeListener(self)

    def deleteRCAlignmentStorage(self, line):
        # type: (NewLine) -> None
        del self.rc_alignments[line.id]
        del self.rc_alignments[line.rc.id]

    def deleteSelfAlignmentStorage(self, line):
        # type: (NewLine) -> None
        del self.auto_alignments[line.id]
        del self.auto_alignments[line.rc.id]

    def printAll(self, out):
        # type: (BinaryIO) -> None
        for line in self.lines:
            out.write(str(list(self.allInter(line.asSegment()))) + "\n")
