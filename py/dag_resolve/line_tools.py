from typing import Generator, Dict, Set, Optional

from common import basic
from common.SeqIO import NamedSequence
from dag_resolve.phasing import Phasing
from dag_resolve.repeat_graph import Edge, Graph, Vertex
from dag_resolve.sequences import Consensus, ReadCollection, Contig, ContigCollection, AlignmentPiece, AlignedRead, \
    Segment


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
    def __init__(self, edge, rc = None):
        # type: (Edge, Optional[Line]) -> None
        assert edge.info.unique
        self.reads = None # type: ReadCollection
        self.rc = None # type: Line
        self.knot = None # type: Knot
        self.id = edge.id
        if rc is None:
            rc = Line(edge.rc, self)
        self.rc = rc # type: Line
        self.consensus = Consensus(edge.seq, [1000] * len(edge.seq))
        Contig.__init__(self, edge.seq, edge.id, None, self.rc)
        self.chain = [AlignmentPiece(self.asSegment(), edge.asSegment(), "=")] # type: list[AlignmentPiece]
        self.reads = ReadCollection(ContigCollection([self]))
        self.new_reads = [] # type: list[AlignedRead]
        self.listeners = []
        self.centerPos = LinePosition(self, len(edge) / 2)
        self.discarded = []

    def setConsensus(self, consensus):
        self.consensus = consensus
        self.seq = consensus.seq
        self.rc.consensus = consensus.RC()
        self.rc.seq = self.rc.consensus.seq

    def cut(self, pos):
        self.extendRight(Consensus("", []), pos)

    def extendRight(self, consensus, pos = None):
        # type: (Consensus, Optional[int]) -> None
        if pos is None:
            pos = len(self.seq)
        elif pos < 0:
            pos = len(self.seq) + pos
        # assert pos >= self.centerPos.pos
        if pos != len(self.seq):
            self.notifyCutRight(pos)
            self.rc.notifyCutLeft(len(self) - pos)
        if len(consensus) != 0:
            self.notifyExtendRight(len(consensus))
            self.rc.notifyExtendLeft(len(consensus))
        old_len = len(self)
        self.setConsensus(self.consensus.merge(consensus, pos))
        for read in self.rc.reads:
            for al in read.alignments:
                if al.seg_to.contig == self.rc:
                    al.seg_to = al.seg_to.shift(len(self) - old_len)
        self.invalidateReadsAfter(max(self.centerPos.pos, pos - 500))
        self.cutAlignments(pos)
        self.fixRC()

    def fixRC(self):
        self.rc.seq = basic.RC(self.seq)
        self.rc.chain = [al.RC() for al in self.chain[::-1]]

    def addRead(self, read):
        # type: (AlignedRead) -> None
        self.reads.add(read)
        self.rc.reads.add(read.rc)

    def invalidateReadsAfter(self, pos):
        assert pos >= self.centerPos.pos
        for read in self.reads.inter(self.suffix(pos)):
            self.invalidateRead(read, self.suffix(pos))

    def invalidateRead(self, read, seg):
        # type: (AlignedRead, Segment) -> None
        read.invalidate(seg)
        self.new_reads.append(read)
        # self.rc.new_reads.append(read.rc)

    def removeRead(self, read):
        # type: (AlignedRead) -> None
        read.removeContig(self)
        self.reads.remove(read)
        self.rc.reads.remove(read.rc)
        self.discarded.append(read)
        self.rc.discarded.append(read)

    def cutAlignments(self, pos):
        while len(self.chain) > 0 and self.chain[-1].seg_from.right > pos:
            if self.chain[-1].seg_from.left >= pos:
                self.chain.pop()
            else:
                matchingSequence = self.chain[-1].matchingSequence(False).reduceQuery(0, pos)
                self.chain[-1] = AlignmentPiece(matchingSequence.SegFrom(self),
                                                matchingSequence.SegTo(self.chain[-1].seg_to.contig),
                                                matchingSequence.cigar())
        self.fixRC()

    def addAlignment(self, piece):
        self.cutAlignments(piece.seg_from.left)
        if self.chain[-1].connects(piece, 5):
            self.chain[-1] = self.chain[-1].merge(piece)
        else:
            self.chain.append(piece)
        self.fixRC()

    def isSimpleLoop(self):
        return self.knot is not None and len(self.chain) == 1 and len(self.rc.chain) == 1 and self.knot.line2.rc.id == self.id

    # def freezeTail(self):
    #     tail_seq = self.tail.alignment.alignedSequence()
    #     new_segment = LineSegment(self, self.chain[-1].pos + 1, self.tail.edge, tail_seq,
    #                               self.tail.phasing, self.tail.reads)
    #     self.chain.append(new_segment)
    #     self.tail = None

    def findEdge(self, edge):
        # type: (Edge) -> Generator[AlignmentPiece]
        for al in self.chain:
            if al.seg_to.contig == edge:
                yield al

    def __str__(self):
        return "Line:" + str(self.id) + "(" + str(len(self)) + "," + str(self.chain[-1].seg_from.right) + \
               "):[" + ",".join(map(str, self.chain)) + "]"

    # def rcStr(self):
    #     return "[" + ",".join(map(lambda seg: str(seg.edge.rc.id), self.chain[::-1])) + "]"

    def notifyExtendLeft(self, l):
        for listener in self.listeners:
            listener.fireExtendLeft(l)

    def notifyCutRight(self, pos):
        for listener in self.listeners:
            listener.fireCutRight(pos)

    def notifyExtendRight(self, l):
        for listener in self.listeners:
            listener.fireExtendRight(l)


    def notifyCutLeft(self, pos):
        for listener in self.listeners:
            listener.fireCutLeft(pos)


class LinePosition:
    def __init__(self, line, pos, invalidateOnCut = False):
        # type: (Line, int, bool) -> LinePosition
        self.line = line
        self.pos = pos
        self.line.listeners.append(self)
        self.invalidateOnCut = invalidateOnCut
        print "Line position:", line.id, pos

    def fireCutRight(self, pos):
        print "Changing position. Old:", self.line.id, self.pos
        if self.pos is None:
            return
        if self.pos > pos:
            if self.invalidateOnCut:
                self.pos = None
            else:
                self.pos = pos
        print "Changed position. New:", self.line.id, self.pos

    def fireCutLeft(self, pos):
        print "Changing position. Old:", self.line.id, self.pos
        if self.pos is None:
            return
        if self.pos < pos:
            if self.invalidateOnCut:
                self.pos = None
            else:
                self.pos = max(0, self.pos - pos)
        print "Changed position. New:", self.line.id, self.pos

    def valid(self):
        return self.pos is not None

    def fireExtendLeft(self, l):
        if self.valid():
            self.pos += l

    def fireExtendRight(self, l):
        pass

    def suffix(self):
        assert self.valid()
        print "Suffix:", self.line.id, self.pos
        return self.line.suffix(self.pos)

    def prefix(self):
        assert self.valid()
        print "Prefix:", self.line.id, self.pos
        return self.line.prefix(self.pos)

    def close(self):
        self.line.listeners.remove(self)

class LineStorage:
    def __init__(self, g):
        # type: (Graph) -> LineStorage
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
                line = Line(edge1)
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