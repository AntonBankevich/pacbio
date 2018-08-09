from dag_resolve import graph
from dag_resolve import sequences
from typing import Generator

class LineSegment:
    def __init__(self, line, pos, edge, seq, reads):
        # type: (Line, int, graph.Edge, basestring, sequences.ReadCollection) -> LineSegment
        self.line = line
        self.pos = pos
        self.edge = edge
        self.seq = seq
        self.reads = reads

class Line:
    def __init__(self, edge):
        # type: (graph.Edge) -> Line
        assert edge.info.unique
        self.chain = [LineSegment(self, 0, edge, edge.seq, edge.reads)]

    def extendRight(self, edge, seq, reads):
        # type: (graph.Edge, basestring, sequences.ReadCollection) -> None
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, reads))

    def extendLeft(self, edge, seq, reads):
        # type: (graph.Edge, basestring, sequences.ReadCollection) -> None
        self.chain.insert(0, LineSegment(self, self.chain[0].pos - 1, edge, seq, reads))

    def rightSegment(self):
        # type: () -> LineSegment
        return self.chain[-1]

    def leftSegment(self):
        # type: () -> LineSegment
        return self.chain[0]

    def find(self, edge):
        # type: (graph.Edge) -> Generator[LineSegment]
        for seg in self.chain:
            if seg.edge.id == edge.id:
                yield seg

    def __getitem__(self, item):
        # type: (int) -> LineSegment
        return self.chain[item - self.chain[0].pos]

class LineStorage:
    def __init__(self, g):
        # type: (graph.Graph) -> LineStorage
        self.g = g
        self.lines = []
        self.resolved_edges = dict()
        for edge in g.E.values():
            if edge.info.unique:
                self.addLine(Line(edge))

    def addLine(self, line):
        # type: (Line) -> None
        for seg in line.chain:
            self.addSegment(seg)
        self.lines.append(line)

    def addSegment(self, seg):
        # type: (LineSegment) -> None
        if seg.edge.id not in self.segments:
            self.resolved_edges = []
        self.resolved_edges[seg.edge.id].append(seg)

    def resolvableVerticesLeft(self):
        # type: () -> Generator[graph.Vertex]
        for v in self.g.V:
            resolved = True
            for e in v.inc:
                if e not in self.resolved_edges:
                    resolved = False
                    break
            if resolved:
                yield v

    def extendRight(self, line, edge, seq, reads):
        # type: (Line, graph.Edge, basestring, sequences.ReadCollection) -> None
        line.extendRight(edge, seq, reads)
        self.addSegment(line.rightSegment())





