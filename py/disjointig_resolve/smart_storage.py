from typing import Optional, Iterator, List, Any

from common.alignment_storage import AlignmentPiece, Correction
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment


class LineListener:
    def __init__(self, rc):
        # type: (LineListener) -> None
        self.rc = rc

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, NamedSequence, str) -> None
        pass

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, NamedSequence, int) -> None
        pass

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        pass

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
        pass

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        pass

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        pass


class SmartStorage(LineListener):
    def __init__(self, rc=None):
        # type: (Optional[SmartStorage]) -> None
        self.items = None  # List
        if rc is None:
            rc = SmartStorage(self)
            self.items = []
        LineListener.__init__(self, rc)
        self.sorted = False
        self.key = lambda item: item

    def isCanonical(self):
        return self.items is not None

    def isSorted(self):
        if self.isCanonical():
            return self.sorted
        else:
            return self.rc.isSorted()

    def add(self, seg):
        assert False
        pass

    def addAll(self, items):
        # type: (List) -> None
        for item in items:
            self.add(item)

    def sort(self):
        # type: () -> None
        if self.isCanonical():
            if not self.sorted:
                self.items = sorted(self.items, key=self.key)
                self.sorted = True
        else:
            self.rc.sort()

    def save(self, handler):
        # type: (TokenWriter) -> None
        if self.isCanonical():
            handler.writeIntLine(len(self.items))
            for item in self.items:
                item.save(handler)
        else:
            handler.writeIntLine(0)

    def __len__(self):
        if self.isCanonical():
            return len(self.items)
        else:
            return len(self.rc)


# Collection of segments where no segment contains another.
# The reason is that sortings are the same for strait and rc
class SegmentStorage(SmartStorage):
    def __init__(self, rc=None):
        # type: (Optional[SegmentStorage]) -> None
        self.items = None  # List[Segment]
        if rc is None:
            rc = SegmentStorage(self)
            self.items = []  # List[Segment]
        SmartStorage.__init__(self, rc)
        self.key = lambda seg: seg.left

    def __getitem__(self, item):
        # type: (int) -> Segment
        if self.isCanonical():
            return self.items[item]
        else:
            return self.rc.__getitem__(-1 - item).RC()

    def __iter__(self):
        # type: () -> Iterator[Segment]
        self.sort()
        if self.isCanonical():
            for item in self.items:
                yield item
        else:
            for item in self.rc.items[::-1]:
                yield item.RC()

    def add(self, seg):
        if self.isCanonical():
            self.items.append(seg)
            self.sorted = False
        else:
            self.rc.add(seg.RC())

    def isIn(self, seg):
        # type: (Segment) -> bool
        if self.isCanonical():
            for seg1 in self.items:
                if seg1.contains(seg):
                    return True
            return False
        else:
            return self.rc.isIn(seg.RC())

    # Merge all segments that have intersection of at least inter_size and remove all subsegments
    # After this operation segment ordering is the same when we sort by left bound as if we sort by right bounds
    def mergeSegments(self, inter_size = 0):
        # type: (int) -> None
        if self.isCanonical():
            self.sort()
            if len(self.items) == 0:
                return
            res = [self.items[0]]  # type: List[Segment]
            for seg in self.items[1:]:  # type: Segment
                if seg.left <= res[-1].right - inter_size or seg.right <= res[-1].right:
                    res[-1].right = max(res[-1].right, seg.right)
                else:
                    res.append(seg)
            self.items = res
        else:
            self.rc.mergeSegments(inter_size)

    def cap(self, seg):
        # type: (Segment) -> List[Segment]
        if self.isCanonical():
            res = []
            for seg1 in self.items:
                if seg.inter(seg1):
                    res.append(seg.cap(seg1))
            return res
        else:
            return [seg1.RC() for seg1 in self.rc.cap(seg.RC())]

    def makeCanonical(self):
        if self.isCanonical():
            return
        self.items = [seg.RC() for seg in self.items[::-1]]
        self.rc.items = None
        self.sorted = self.rc.sorted

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, NamedSequence, str) -> None
        self.makeCanonical()


    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, NamedSequence, int) -> None
        self.makeCanonical()
        self.items = [seg for seg in self.items if seg.left < pos] # type: List[Segment]
        if len(self.items) > 0 and self.items[-1].right > pos:
            assert self.items[-1].contig == line
            self.items[-1] = self.items[-1].cap(line.segment(line.left(), pos))

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, correction):
        # type: (Correction) -> None
        self.makeCanonical()
        self.items = correction.mapSegmentsUp(self.items)

    def fireAfterExtendRight(self, line, seq):
        # type: (Any, str) -> None
        pass

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        pass

    def fireAfterCorrect(self, line):
        # type: (Any) -> None
        self.makeCanonical()
        self.items = [seg.contigAsSegment(line.segment(line.left(), line.right())) for seg in self.items]

    def load(self, handler, contig):
        # type: (TokenReader, NamedSequence) -> None
        n = handler.readInt()
        for i in range(n):
            self.add(Segment.load(handler, contig))


class AlignmentStorage(SmartStorage):
    def __init__(self, rc=None):
        # type: (Optional[AlignmentStorage]) -> None
        self.items = None  # type: List[AlignmentPiece]
        if rc is None:
            rc = AlignmentStorage(self)
            self.items = []
        SmartStorage.__init__(self, rc)
        self.key = lambda al: al.seg_to.left

    def __getitem__(self, item):
        # type: (int) -> AlignmentPiece
        if self.isCanonical():
            return self.items[item]
        else:
            return self.rc.__getitem__(-1 - item).rc

    def __iter__(self):
        # type: () -> Iterator[AlignmentPiece]
        self.sort()
        if self.isCanonical():
            for item in self.items:
                yield item
        else:
            for item in self.rc.items[::-1]:
                yield item.rc

    def add(self, al):
        if self.isCanonical():
            self.items.append(al)
            self.sorted = False
        else:
            self.rc.add(al.rc)

    def fireBeforeExtendRight(self, line, seq):
        # IMPLEMENT
        pass

    def fireBeforeCutRight(self, line, pos):
        # IMPLEMENT
        pass

    def fireBeforeCorrect(self, alignments):
        # IMPLEMENT
        pass

    def allInter(self, seg):
        # type: (Segment) -> List[AlignmentPiece]
        result = []
        for al in self: # type: AlignmentPiece
            if al.seg_to.inter(seg):
                result.append(al)
        return result

    def allContaining(self, seg):
        # type: (Segment) -> List[AlignmentPiece]
        result = []
        for al in self: # type: AlignmentPiece
            if al.seg_to.contains(seg):
                result.append(al)
        return result

    def load(self, handler, collection_from, collection_to):
        # type: (TokenReader, Any, Any) -> None
        n = handler.readInt()
        for i in range(n):
            self.add(AlignmentPiece.load(handler, collection_from, collection_to))