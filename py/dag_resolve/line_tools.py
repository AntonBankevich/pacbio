from typing import Generator, Dict, Set

from alignment.align_tools import AlignedSequences
from common import basic
from dag_resolve.phasing import Phasing
from dag_resolve.repeat_graph import Edge, Graph, Vertex
from dag_resolve.sequences import Segment, Consensus, ReadCollection


class LineSegment:
    def __init__(self, line, pos, edge, seq, phasing, reads):
        # type: (Line, int, Edge, str, Phasing, ReadCollection) -> LineSegment
        self.line = line
        self.pos = pos
        self.edge = edge
        self.seq = seq
        self.phasing = phasing
        self.reads = reads

class Line:
    def __init__(self, edge):
        # type: (Edge) -> Line
        assert edge.info.unique
        self.chain = [LineSegment(self, 0, edge, edge.seq, Phasing(), edge.reads)] # type: list[LineSegment]
        # self.nextLine = None
        self.id = edge.id
        self.tail = None # type: LineTail
        self.rc = None # type: Line
        self.knot = None # type: Knot

    def extendRight(self, edge, seq, phasing, reads):
        # type: (Edge, str, Phasing, ReadCollection) -> None
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, phasing, reads))

    def extendLeft(self, edge, seq, phasing, reads):
        # type: (Edge, str, Phasing, ReadCollection) -> None
        self.chain.insert(0, LineSegment(self, self.chain[0].pos - 1, edge, seq, phasing, reads))

    def isSimpleLoop(self):
        return self.knot is not None and len(self.chain) == 1 and len(self.rc.chain) == 1 and self.knot.line2.rc.id == self.id

    def freezeTail(self):
        tail_seq = self.tail.alignment.alignedSequence()
        new_segment = LineSegment(self, self.chain[-1].pos + 1, self.tail.edge, tail_seq,
                                  self.tail.phasing, self.tail.reads)
        self.chain.append(new_segment)
        self.tail = None

    def rightSegment(self):
        # type: () -> LineSegment
        return self.chain[-1]

    def leftSegment(self):
        # type: () -> LineSegment
        return self.chain[0]

    def find(self, edge):
        # type: (Edge) -> Generator[LineSegment]
        for seg in self.chain:
            if seg.edge.id == edge.id:
                yield seg

    def __str__(self):
        return "[" + ",".join(map(lambda seg: str(seg.edge.id), self.chain)) + "]"

    def rcStr(self):
        return "[" + ",".join(map(lambda seg: str(seg.edge.rc.id), self.chain[::-1])) + "]"

    def __getitem__(self, item):
        # type: (int) -> LineSegment
        return self.chain[item - self.chain[0].pos]
    
    def merge(self, other):
        # type: (Line) -> None
        self.nextLine = other

class LineStorage:
    def __init__(self, g):
        # type: (Graph) -> LineStorage
        self.g = g
        self.lines = [] #type: list[Line]
        self.resolved_edges = set() #type: Set[int]
        for edge1, edge2 in g.unorientedEdges():
            if edge1.info.unique:
                self.resolved_edges.add(edge1.id)
                self.resolved_edges.add(edge2.id)
                line1 = Line(edge1)
                line2 = Line(edge2)
                line1.rc = line2
                line2.rc = line1
                self.lines.append(line1)
                self.lines.append(line2)

    def isResolvableLeft(self, v):
        # type: (Vertex) -> bool
        if len(v.inc) == 0:
            return False
        for e in v.inc:
            if e.id not in self.resolved_edges:
                return False
        return True

    def edgeLines(self, edge):
        # type: (Edge) -> Generator[LineSegment]
        assert edge.id in self.resolved_edges
        for line in self.lines:
            for segment in line.chain:
                if segment.edge == edge:
                    yield segment

    def edgeTails(self, edge):
        # type: (Edge) -> Generator[LineTail]
        for line in self.lines:
            if line.tail is not None and line.tail.edge == edge:
                yield line.tail

    def printToFile(self, handler):
        # type: (file) -> None
        for line in self.lines:
            handler.write(line.__str__() + "\n")


class LineTail:
    def __init__(self, line, edge, tail_consensus, read_collection):
        # type: (Line, Edge, Consensus, ReadCollection) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus
        self.reads = read_collection
        self.alignment = None # type: AlignedSequences
        self.reverse_alignment = None # type: AlignedSequences
        self.phasing = None

    def edgeSegment(self):
        assert self.alignment is not None
        return Segment(self.edge, self.alignment.first, self.alignment.last + 1)

    def __len__(self):
        # type: () -> int
        return self.tail_consensus.__len__()


class Knot:
    def __init__(self, line1, line2, seq, reads):
        # type: (Line, Line, str, ReadCollection) -> Knot
        self.line1 = line1
        self.line2 = line2
        self.seq = seq
        self.reads = reads