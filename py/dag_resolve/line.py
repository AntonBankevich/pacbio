from dag_resolve import repeat_graph
from dag_resolve import sequences
from typing import Generator

class DivergenceState:
    def __init__(self, div, neighborhood):
        # type: (Divergence, str) -> DivergenceState
        self.divergence = div
        self.neighborhood = neighborhood

    def isAmbiguous(self):
        # type: () -> bool
        return self.neighborhood is None

class Divergence:
    def __init__(self, edge, pos, states):
        # type: (repeat_graph.Edge, int, list[basestring]) -> Divergence
        self.edge = edge
        self.pos = pos
        self.states = []
        self.ambiguous = DivergenceState(self, None)
        for state in states:
            self.states.append(DivergenceState(self, state))


class LineSegment:
    def __init__(self, line, pos, edge, seq, states, reads):
        # type: (Line, int, repeat_graph.Edge, str, list[DivergenceState], sequences.ReadCollection) -> LineSegment
        self.line = line
        self.pos = pos
        self.edge = edge
        self.seq = seq
        self.states = states
        self.reads = reads

class Line:
    def __init__(self, edge):
        # type: (repeat_graph.Edge) -> Line
        assert edge.info.unique
        self.chain = [LineSegment(self, 0, edge, edge.seq, [], edge.reads)] # type: list[LineSegment]

    def extendRight(self, edge, seq, reads):
        # type: (repeat_graph.Edge, basestring, sequences.ReadCollection) -> None
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, reads))

    def extendLeft(self, edge, seq, reads):
        # type: (repeat_graph.Edge, basestring, sequences.ReadCollection) -> None
        self.chain.insert(0, LineSegment(self, self.chain[0].pos - 1, edge, seq, reads))

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

    def __getitem__(self, item):
        # type: (int) -> LineSegment
        return self.chain[item - self.chain[0].pos]

class LineStorage:
    def __init__(self, g):
        # type: (repeat_graph.Graph) -> LineStorage
        self.g = g
        self.lines = [] #type: list[Line]
        self.resolved_edges = dict() #type: dict[int, LineSegment]
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
        if seg.edge.id not in self.resolved_edges:
            self.resolved_edges[seg.edge.id] = []
        self.resolved_edges[seg.edge.id].append(seg)

    def isResolvableLeft(self, v):
        # type: (repeat_graph.Vertex) -> bool
        if len(v.inc) == 0:
            return False
        for e in v.inc:
            if e.id not in self.resolved_edges:
                return False
        return True

    def extendRight(self, line, edge, seq, reads):
        # type: (Line, repeat_graph.Edge, basestring, sequences.ReadCollection) -> None
        line.extendRight(edge, seq, reads)
        self.addSegment(line.rightSegment())

    def printToFile(self, handler):
        # type: (file) -> None
        for line in self.lines:
            for seg in line.chain:
                handler.write(str(seg.edge.id) + " ")
            handler.write("\n")





