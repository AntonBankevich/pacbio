from typing import Generator, Dict, Set, Optional, Iterable

from alignment.align_tools import Aligner
from dag_resolve.repeat_graph import Edge, Graph, Vertex
from common.sequences import Consensus, ReadCollection, Contig, ContigCollection, AlignmentPiece, AlignedRead


# class PseudoLineSegment:
#     def __init__(self, line, eseg, seq, reads):
#         # type: (Line, EdgeSegment, str, AlignmentCollection) -> PseudoLineSegment
#         self.line = line
#         self.edgeSegment = eseg
#         self.seq = seq
#         self.reads = reads
#         self.edgeSegment.psegments.append(self)

# class EdgeSegment(NamedSequence):
#     def __init__(self, seq, id, edge, alignment = None, rc = None):
#         # type: (str, int, Edge, AlignmentPiece, Optional[EdgeSegment]) -> EdgeSegment
#         NamedSequence.__init__(self, seq, id)
#         self.edge = edge
#         if alignment is not None:
#             assert seq == alignment.seg_from.contig.seq
#             alignment.seg_from.contig = self
#         self.alignment = alignment
#         self.psegments = []
#         if rc is None:
#             if alignment is None:
#                 rc_alignment = None
#             else:
#                 rc_alignment = alignment.RC()
#             self.rc = EdgeSegment(basic.RC(seq), -id, edge.rc, rc_alignment, self)
#         else:
#             self.rc = rc # type: EdgeSegment
#
#     def setAlignment(self, alignment):
#         # type: (AlignmentPiece) -> None
#         self.alignment = alignment
#         self.rc.alignment = alignment.RC()


class Line(Contig):
    def __init__(self, edge, aligner, rc = None):
        # type: (Edge, Aligner, Optional[Line]) -> None
        assert edge.info.unique
        self.reads = None # type: ReadCollection
        self.baseEdge = edge
        self.rc = None # type: Line
        self.knot = None # type: Knot
        self.id = edge.id
        if rc is None:
            rc = Line(edge.rc, aligner, self)
        self.rc = rc # type: Line
        self.consensus = Consensus(edge.seq, [1000] * len(edge.seq))
        Contig.__init__(self, edge.seq, edge.id, None, self.rc)
        self.chain = [AlignmentPiece(self.asSegment(), edge.asSegment(), "=")] # type: list[AlignmentPiece]
        self.reads = ReadCollection(ContigCollection([self]))
        self.invalidated_reads = [] # type: list[AlignedRead]
        self.listeners = []
        self.centerPos = LinePosition(self, len(edge) / 2)
        self.discarded = []
        self.aligner = aligner

    def setConsensus(self, consensus):
        self.consensus = consensus
        self.seq = consensus.seq
        self.rc.consensus = consensus.RC()
        self.rc.seq = self.rc.consensus.seq

    def extendRight(self, consensus, pos = None):
        # type: (Consensus, Optional[int]) -> None
        if pos is None:
            pos = len(self.seq)
        elif pos < 0:
            pos = len(self.seq) + pos
        assert pos >= self.centerPos.pos
        if pos != len(self.seq):
            self.notifyCutRightBefore(pos)
            cut_len = len(self) - pos
            self.notifyCutLeftBefore(cut_len)
            self.setConsensus(self.consensus.cut(length=pos))
            assert len(self) == pos
            self.notifyCutRightAfter(pos)
            self.rc.notifyCutLeftAfter(cut_len)
        if len(consensus) != 0:
            self.rc.notifyExtendLeftBefore(len(consensus))
            self.notifyExtendRightBefore(len(consensus))
            self.setConsensus(self.consensus.merge(consensus))
            self.rc.notifyExtendLeftAfter(len(consensus))
            self.notifyExtendRightAfter(len(consensus))

    def addRead(self, read):
        # type: (AlignedRead) -> None
        self.reads.add(read)
        self.rc.reads.add(read.rc)

    def addReads(self, reads):
        # type: (Iterable[AlignedRead]) -> None
        for read in reads:
            self.addRead(read)

    def invalidateRead(self, read):
        # type: (AlignedRead) -> None
        read.invalidate(self.centerPos.suffix())
        self.invalidated_reads.append(read)

    def removeRead(self, read):
        # type: (AlignedRead) -> None
        read.removeContig(self)
        self.reads.remove(read)
        self.rc.reads.remove(read.rc)
        self.discarded.append(read)
        self.rc.discarded.append(read.rc)

    def cutAlignments(self, pos):
        while len(self.chain) > 0 and self.chain[-1].seg_from.right > pos:
            if self.chain[-1].seg_from.left >= pos:
                self.chain.pop()
            else:
                matchingSequence = self.chain[-1].matchingSequence(False).reduceQuery(0, pos)
                self.chain[-1] = AlignmentPiece(matchingSequence.SegFrom(self),
                                                matchingSequence.SegTo(self.chain[-1].seg_to.contig),
                                                matchingSequence.cigar())
        self.rc.chain = [al.RC() for al in self.chain[::-1]]

    def addAlignment(self, piece):
        self.cutAlignments(piece.seg_from.left)
        if self.chain[-1].connects(piece, 5):
            self.chain[-1] = self.chain[-1].merge(piece)
        else:
            self.chain.append(piece)
        self.rc.chain = [al.RC() for al in self.chain[::-1]]

    def isSimpleLoop(self):
        return self.knot is not None and len(self.chain) == 1 and len(self.rc.chain) == 1 and self.knot.line2.rc.id == self.id

    def findEdge(self, edge):
        # type: (Edge) -> Generator[AlignmentPiece]
        for al in self.chain:
            if al.seg_to.contig == edge:
                yield al

    def __str__(self):
        return "Line:" + str(self.id) + "(" + str(len(self)) + "," + str(self.chain[-1].seg_from.right) + \
               "):[" + ",".join(map(str, self.chain)) + "]"

    def fixLineAlignments(self):
        # type: () -> None
        print "Fixing alignments."
        to_fix = ReadCollection(ContigCollection([Contig(self.centerPos.suffix().Seq(), "half")]))
        to_fix.extendClean(self.invalidated_reads)
        still_invalidated_reads = []
        self.aligner.alignReadCollection(to_fix)
        seg_dict = {"half": self.centerPos.suffix()}
        to_fix.contigsAsSegments(seg_dict)
        to_fix_center = []
        for read in self.invalidated_reads: # type: AlignedRead
            if read.inter(self.centerPos.suffix()):
                continue
            aligned_read = to_fix[read.id]
            has_contradicting = False
            noncontradicting_al = None
            has_center = False
            for al in aligned_read.alignmentsTo(self.asSegment()):
                if not al.contradicting(self.asSegment()) and al.seg_to.left - al.seg_from.left > self.centerPos.pos:
                    noncontradicting_al = al
                elif not al.contradicting(self.centerPos.suffix()):
                    has_center = True
                else:
                    has_contradicting = True
            if noncontradicting_al is not None:
                print "Adding read", read, "to line", self, "with alignment", noncontradicting_al
                assert read in self.reads
                print "Adding read-to-line alignment", noncontradicting_al.changeQuery(read)
                read.addAlignment(noncontradicting_al.changeQuery(read))
            else:
                if has_contradicting:
                    print "REMOVING READ!!!", aligned_read, "since it only has contradicting alignments to the line"
                    self.removeRead(read)
                elif has_center:
                    to_fix_center.append(read)
                else:
                    still_invalidated_reads.append(read)
        self.invalidated_reads = still_invalidated_reads
        if len(to_fix_center) == 0:
            return
        to_fix = ReadCollection(ContigCollection([self]))
        to_fix.extendClean(to_fix_center)
        self.aligner.alignReadCollection(to_fix)
        for read in to_fix_center:
            aligned_read = to_fix[read.id]
            has_contradicting = False
            noncontradicting_al = None
            for al in aligned_read.alignmentsTo(self.asSegment()): # type: AlignmentPiece
                if al.contradicting(self.asSegment()):
                    has_contradicting = True
                elif al.seg_to.inter(self.centerPos.suffix()):
                    noncontradicting_al = al
            if noncontradicting_al is not None:
                print "Adding read", read, "to line", self, "with alignment", noncontradicting_al
                assert read in self.reads
                print "Adding read-to-line alignment", noncontradicting_al.changeQuery(read)
                read.addAlignment(noncontradicting_al.changeQuery(read))
            else:
                if has_contradicting:
                    print "REMOVING READ!!!", aligned_read, "since it only has contradicting alignments to the line"
                    self.removeRead(read)


    def notifyExtendLeftBefore(self, l):
        for listener in self.listeners:
            listener.fireExtendLeftBefore(l)

    def notifyExtendLeftAfter(self, l):
        for listener in self.listeners:
            listener.fireExtendLeftAfter(l)
        for read in self.reads:
            for al in read.alignments:
                if al.seg_to.contig == self:
                    al.seg_to = al.seg_to.shift(l)
        self.chain = [al.RC() for al in self.rc.chain[::-1]]

    def notifyCutRightBefore(self, pos):
        for listener in self.listeners:
            listener.fireCutRightBefore(pos)
        for read in self.reads.inter(self.suffix(pos)):
            self.invalidateRead(read)

    def notifyCutRightAfter(self, pos):
        for listener in self.listeners:
            listener.fireCutRightAfter(pos)
        self.cutAlignments(pos)

    def notifyExtendRightBefore(self, l):
        for listener in self.listeners:
            listener.fireExtendRightBefore(l)

    def notifyExtendRightAfter(self, l):
        for listener in self.listeners:
            listener.fireExtendRightAfter(l)
        self.fixLineAlignments()

    def notifyCutLeftBefore(self, pos):
        for listener in self.listeners:
            listener.fireCutLeftBefore(pos)

    def notifyCutLeftAfter(self, pos):
        for listener in self.listeners:
            listener.fireCutLeftAfter(pos)
        # print "CutLeftAfter. Position:", pos
        for read in self.reads:
            # print read
            for al in read.alignments:
                if al.seg_to.contig == self:
                    # print al
                    al.seg_to = al.seg_to.shift(-pos)


class LinePosition:
    def __init__(self, line, pos, invalidateOnCut = False):
        # type: (Line, int, bool) -> LinePosition
        self.line = line
        self.pos = pos
        self.line.listeners.append(self)
        self.invalidateOnCut = invalidateOnCut

    def fireCutRightBefore(self, pos):
        pass

    def fireCutRightAfter(self, pos):
        if self.pos is None:
            return
        if self.pos > pos:
            if self.invalidateOnCut:
                self.pos = None
            else:
                self.pos = pos

    def fireCutLeftBefore(self, pos):
        pass

    def fireCutLeftAfter(self, pos):
        if self.pos is None:
            return
        if self.pos < pos:
            if self.invalidateOnCut:
                self.pos = None
            else:
                self.pos = 0
        else:
            self.pos = self.pos - pos

    def valid(self):
        return self.pos is not None

    def fireExtendLeftBefore(self, l):
        pass

    def fireExtendLeftAfter(self, l):
        if self.valid():
            self.pos += l

    def fireExtendRightBefore(self, l):
        pass

    def fireExtendRightAfter(self, l):
        pass

    def suffix(self):
        assert self.valid()
        return self.line.suffix(self.pos)

    def prefix(self):
        assert self.valid()
        return self.line.prefix(self.pos)

    def close(self):
        self.line.listeners.remove(self)

class LineStorage:
    def __init__(self, g, aligner):
        # type: (Graph, Aligner) -> LineStorage
        self.g = g
        self.lines = [] #type: list[Line]
        self.edgeLines = dict() # type: Dict[str, list[Line]]
        for edge in g.E.values():
            self.edgeLines[edge.id] = []
        self.resolved_edges = set() #type: Set[int]
        for edge1, edge2 in g.unorientedEdges():
            if edge1.info.unique:
                assert edge2.info.unique
                self.resolved_edges.add(edge1.id)
                self.resolved_edges.add(edge2.id)
                line = Line(edge1, aligner)
                self.edgeLines[edge1.id].append(line)
                self.edgeLines[edge2.id].append(line.rc)
                self.lines.append(line)
                self.lines.append(line.rc)
        self.reads = ReadCollection(ContigCollection(self.lines))
        for read in self.g.reads:
            self.reads.addNewRead(read)
        # for line in self.lines:
        #     for read in self.g.E[line.id].reads:
        #         line.new_reads.append(self.reads[read.id])

    def getLine(self, eid):
        for line in self.lines:
            if line.id == eid:
                return line
        return None

    def getEdgeLines(self, edge):
        # type: (Edge) -> list[Line]
        if edge.id not in self.edgeLines:
            return []
        return self.edgeLines[edge.id]

    def extendLine(self, line, edge):
        if edge.id not in self.edgeLines:
            self.edgeLines[edge.id] = []
        if not edge.unique():
            self.edgeLines[edge.id].append(line)


    def isResolvableLeft(self, v):
        # type: (Vertex) -> bool
        if len(v.inc) == 0 or len(v.out) == 0:
            return False
        for e in v.inc:
            if e.id not in self.resolved_edges:
                return False
        return True

    def printToFile(self, handler):
        # type: (file) -> None
        for line in self.lines:
            handler.write(line.__str__() + "\n")


# class LineTail:
#     def __init__(self, line, edge, tail_consensus, read_collection):
#         # type: (Line, Edge, Consensus, ReadCollection) -> LineTail
#         self.line = line
#         self.edge = edge
#         self.tail_consensus = tail_consensus
#         self.reads = read_collection
#         self.alignment = None # type: AlignedSequences
#         self.reverse_alignment = None # type: AlignedSequences
#         self.phasing = None
#
#     def edgeSegment(self):
#         assert self.alignment is not None
#         return Segment(self.edge, self.alignment.first, self.alignment.last + 1)
#
#     def __len__(self):
#         # type: () -> int
#         return self.tail_consensus.__len__()


class Knot:
    def __init__(self, line1, line2, seq, reads):
        # type: (Line, Line, str, ReadCollection) -> Knot
        self.line1 = line1
        self.line2 = line2
        self.seq = seq
        self.reads = reads