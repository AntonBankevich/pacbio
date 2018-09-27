from typing import Generator, Dict, Set, Optional

from alignment.align_tools import AlignedSequences
from common import basic
from common.SeqIO import NamedSequence
from dag_resolve.phasing import Phasing
from dag_resolve.repeat_graph import Edge, Graph, Vertex
from dag_resolve.sequences import Segment, Consensus, ReadCollection, AlignmentCollection, AlignmentPiece, Contig, \
    AlignedRead, ContigCollection


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
        # type: (Edge) -> Line
        assert edge.info.unique
        self.chain = [AlignmentPiece(self.asSegment(), edge.asSegment(), "=")] # type: list[AlignmentPiece]
        self.reads = edge.reads
        self.rc = None # type: Line
        self.knot = None # type: Knot
        if rc is None:
            rc = Line(edge.rc, self)
        self.rc = rc
        self.consensus = Consensus(edge.seq, [1000] * len(edge.seq))
        Contig.__init__(self, edge.seq, edge.id, self.rc)
        if self.rc.reads is None:
            self.reads = edge.reads.cleanCopy(ContigCollection([self, self.rc]))
            self.rc.reads = self.reads

    def extendRight(self, consensus, pos = None):
        # type: (Consensus, Optional(int)) -> None
        if pos is None:
            pos = len(self.seq)
        self.consensus.merge(consensus, pos)
        self.seq = self.consensus.seq
        self.removeAlignments(pos)
        self.fixRC()

    def fixRC(self):
        self.rc.seq = basic.RC(self.seq)
        self.rc.chain = [al.RC() for al in self.chain[::-1]]

    def removeAlignments(self, pos):
        while len(self.chain) > 0 and self.chain[-1].seg_from.right > pos:
            self.chain.pop()
        self.fixRC()

    def addAlignment(self, piece):
        self.removeAlignments(piece.seg_from.left)
        if self.chain[-1].connects(piece, 0):
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
        return "[" + ",".join(map(lambda al: str(al.seg_to), self.chain)) + "]"

    def rcStr(self):
        return "[" + ",".join(map(lambda seg: str(seg.edge.rc.id), self.chain[::-1])) + "]"

    def __getitem__(self, item):
        # type: (int) -> AlignmentPiece
        return self.chain[item]

class LineStorage:
    def __init__(self, g):
        # type: (Graph) -> LineStorage
        self.g = g
        self.lines = [] #type: list[Line]
        self.edgeLines = dict() # type: Dict[int, list[Line]]
        self.resolved_edges = set() #type: Set[int]
        for edge1, edge2 in g.unorientedEdges():
            if edge1.info.unique:
                assert edge2.info.unique
                self.resolved_edges.add(edge1.id)
                self.resolved_edges.add(edge2.id)
                line = Line(edge1)
                self.lines.append(line)
                self.lines.append(line.rc)

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