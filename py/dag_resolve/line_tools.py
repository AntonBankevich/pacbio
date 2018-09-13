from typing import Generator, Dict

from alignment import align_tools
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

    def extendRight(self, edge, seq, phasing, reads):
        # type: (repeat_graph.Edge, str, Phasing, sequences.ReadCollection) -> None
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, phasing, reads))

    def extendLeft(self, edge, seq, phasing, reads):
        # type: (repeat_graph.Edge, str, Phasing, sequences.ReadCollection) -> None
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
    
    def merge(self, other):
        # type: (Line) -> None
        self.nextLine = other

class LineStorage:
    def __init__(self, g):
        # type: (repeat_graph.Graph) -> LineStorage
        self.g = g
        self.lines = [] #type: list[Line]
        self.resolved_edges = dict() #type: Dict[int, list[LineSegment]]
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


class LineTail:
    def __init__(self, line, edge, tail_consensus, read_collection):
        # type: (Line, repeat_graph.Edge, sequences.Consensus, sequences.ReadCollection) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus
        self.reads = read_collection
        self.alignment = None # type: align_tools.AlignedSequences

    def __len__(self):
        # type: () -> int
        return self.tail_consensus.__len__()