from dag_resolve import repeat_graph
from dag_resolve import sequences
from typing import Generator

class DivergenceState:
    def __init__(self, div, seq, id):
        # type: (Divergence, str) -> DivergenceState
        self.divergence = div
        self.seq = seq
        self.id = id

    def isAmbiguous(self):
        # type: () -> bool
        return self.seq is None

    def char(self):
        if self.isAmbiguous():
            return "-"
        return str(self.id)

    def __eq__(self, other):
        # type: (DivergenceState) -> bool
        return self.seq == other.seq and self.divergence == other.divergence

    def __ne__(self, other):
        return not self.__eq__(other)

class Phasing:
    def __init__(self, states = None):
        # type: (list[DivergenceState]) -> Phasing
        if states is None:
            states = []
        self.states = states

    def add(self, state):
        # type: (DivergenceState) -> DivergenceState
        self.states.append(state)
        return self.states[-1]

    def printToFile(self, handler):
        # type: (file) -> None
        for state in self.states:
            handler.write(state.char())
        handler.write("\n")

    def ambibuousRate(self):
        res = 0
        for state in self.states:
            if state.isAmbiguous():
                res += 1
        return float(state) / len(self.states)

    def called(self):
        res = 0
        for state in self.states:
            if not state.isAmbiguous():
                res += 1
        return res

    def __getitem__(self, item):
        # type: (int) -> DivergenceState
        return self.states[item]

    def __iter__(self):
        # type: () -> Iterator[DivergenceState]
        return self.states.__iter__()

    def __len__(self):
        # type: () -> int
        return self.states.__len__()

class Statistics:
    def __init__(self):
        self.correct = 0 # type: int
        self.wrong = 0 # type: int
        self.ambig = 0 # type: int
        self.called = 0 # type: int

    def __str__(self):
        return str((self.correct, self.wrong, self.ambig, self.called))

class Divergence:
    def __init__(self, edge, pos, states = []):
        # type: (repeat_graph.Edge, tuple[int, int], list[str]) -> Divergence
        self.edge = edge
        self.pos = pos
        self.states = [] # type: list[DivergenceState]
        self.ambiguous = DivergenceState(self, None, -1)
        for state in states:
            self.addState(state)
        self.statistics = Statistics()

    def addState(self, state):
        # type: (str) -> DivergenceState
        self.states.append(DivergenceState(self, state, len(self.states)))
        return self.states[-1]

    def __eq__(self, other):
        # type: (Divergence) -> bool
        return self.edge.id == other.edge.id and self.pos == other.pos

    def __ne__(self, other):
        return not self.__eq__(other)

    def printToFile(self, handler, delim = " "):
        # type: (file) -> None
        handler.write(str(self.edge.id) + ": " + str(self.pos) + "\n")
        for state in self.states:
            handler.write(state.seq + " ")
        handler.write("\n")

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
        self.chain = [LineSegment(self, 0, edge, edge.seq, [], edge.reads)] # type: list[LineSegment]

    def extendRight(self, edge, seq, phasing, reads):
        # type: (repeat_graph.Edge, str, sequences.ReadCollection) -> None
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, phasing, reads))

    def extendLeft(self, edge, seq, reads):
        # type: (repeat_graph.Edge, str, sequences.ReadCollection) -> None
        self.chain.insert(0, LineSegment(self, self.chain[0].pos - 1, edge, seq, phasing, reads))

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

    def shortStr(self):
        return "[" + ",".join(map(lambda seg: str(seg.edge.id), self.chain)) + "]"

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

    def extendRight(self, line, edge, seq, phasing, reads):
        # type: (Line, repeat_graph.Edge, str, Phasing, sequences.ReadCollection) -> None
        line.extendRight(edge, seq, phasing, reads)
        self.addSegment(line.rightSegment())

    def printToFile(self, handler):
        # type: (file) -> None
        for line in self.lines:
            handler.write(line.shortStr() + "\n")





