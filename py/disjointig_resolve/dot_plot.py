from typing import Dict, List

from alignment.align_tools import Aligner
from common.alignment_storage import AlignmentPiece
from disjointig_resolve.accurate_line import NewLineStorage
from disjointig_resolve.smart_storage import LineListener
from common.save_load import TokenWriter, TokenReader
from common.sequences import Segment


class LineDotPlot(LineListener):
    def __init__(self, lines):
        # type: (NewLineStorage) -> None
        LineListener.__init__(self)
        self.lines = lines
        for line in lines:
            line.addListener(self)
        self.alignmentsToFrom = dict() # type: Dict[str, Dict[str, List[AlignmentPiece]]]
        self.alignmentsFromTo = dict() # type: Dict[str, Dict[str, List[AlignmentPiece]]]

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        to_id = al.seg_to.contig.id
        from_id = al.seg_from.contig.id
        if to_id not in self.alignmentsToFrom:
            self.alignmentsToFrom[to_id] = dict()
        if from_id not in self.alignmentsFromTo:
            self.alignmentsFromTo[from_id] = dict()
        if from_id not in self.alignmentsToFrom[to_id]:
            self.alignmentsToFrom[to_id][from_id] = []
        if to_id not in self.alignmentsFromTo[from_id]:
            self.alignmentsFromTo[from_id][to_id] = []
        self.alignmentsToFrom[to_id][from_id].append(al)
        if to_id != from_id:
            self.alignmentsFromTo[from_id][to_id].append(al)

    def getAllAlignments(self, seg):
        # type: (Segment) -> List[AlignmentPiece]
        # IMPLEMENT
        pass

    def construct(self, aligner):
        # type: (Aligner) -> None
        # IMPLEMENT
        pass

    def fireExtendRight(self, line, seq):
        # IMPLEMENT
        pass

    def fireCutRight(self, line, pos):
        # IMPLEMENT
        pass

    def fireCorrect(self, alignments):
        # IMPLEMENT
        pass

    def save(self, handler):
        # type: (TokenWriter) -> None
        n = sum(map(lambda als: sum(map(len, als.values())), self.alignmentsToFrom.values()))
        handler.writeIntLine(n)
        for d1 in self.alignmentsToFrom.values():
            for als in d1.values():
                for al in als:
                    al.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        n = handler.readInt()
        for i in range(n):
            al = AlignmentPiece.load(handler, self.lines, self.lines)
            self.addAlignment(al)

