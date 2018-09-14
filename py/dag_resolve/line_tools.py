from typing import Generator, Dict

from alignment import align_tools
from common import basic
from dag_resolve import repeat_graph
from dag_resolve import sequences
from dag_resolve.phasing import Phasing


class LineSegment:
    def __init__(self, line, pos, edge, seq, phasing, reads):
        # type: (Line, int, repeat_graph.Edge, str, Phasing, sequences.ReadCollection) -> LineSegment
        self.line = line
        self.pos = pos
        self.edge = edge
        self.seq = seq
        self.phasing = phasing
        self.reads = reads

class Line:
    def __init__(self, edge):
        # type: (repeat_graph.Edge) -> Line
        assert edge.info.unique
        self.chain = [LineSegment(self, 0, edge, edge.seq, Phasing(), edge.reads)] # type: list[LineSegment]
        self.nextLine = None
        self.id = edge.id
        self.tail = None # type: LineTail
        self.rc = None

    def extendRight(self, edge, seq, phasing, reads):
        # type: (repeat_graph.Edge, str, Phasing, sequences.ReadCollection) -> None
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, phasing, reads))

    def extendLeft(self, edge, seq, phasing, reads):
        # type: (repeat_graph.Edge, str, Phasing, sequences.ReadCollection) -> None
        self.chain.insert(0, LineSegment(self, self.chain[0].pos - 1, edge, seq, phasing, reads))

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
        # type: (repeat_graph.Edge) -> Generator[LineSegment]
        for seg in self.chain:
            if seg.edge.id == edge.id:
                yield seg

    def __str__(self):
        return "[" + ",".join(map(lambda seg: str(seg.edge.id), self.chain)) + "]"

    def __getitem__(self, item):
        # type: (int) -> LineSegment
        return self.chain[item - self.chain[0].pos]
    
    def merge(self, other):
        # type: (Line) -> None
        self.nextLine = other

class LineStorage:
    def __init__(self, g):
        # type: (repeat_graph.Graph) -> LineStorage
        self.g = g
        self.lines = [] #type: list[Line]
        self.resolved_edges = set() #type: set[int]
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
        # type: (repeat_graph.Vertex) -> bool
        if len(v.inc) == 0:
            return False
        for e in v.inc:
            if e.id not in self.resolved_edges:
                return False
        return True

    def edgeLines(self, edge):
        # type: (repeat_graph.Edge) -> Generator[LineSegment]
        assert edge.id in self.resolved_edges
        for line in self.lines:
            for segment in line.chain:
                if segment.edge == edge:
                    yield segment

    def edgeTails(self, edge):
        # type: (repeat_graph.Edge) -> Generator[LineTail]
        for line in self.lines:
            if line.tail is not None and line.tail.edge == edge:
                yield line.tail

    def printToFile(self, handler):
        # type: (file) -> None
        for line in self.lines:
            handler.write(line.__str__() + "\n")


class LineTail:
    def __init__(self, line, edge, tail_consensus, read_collection):
        # type: (Line, repeat_graph.Edge, sequences.Consensus, sequences.ReadCollection) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus
        self.reads = read_collection
        self.alignment = None # type: align_tools.AlignedSequences
        self.phasing = None

    def __len__(self):
        # type: () -> int
        return self.tail_consensus.__len__()