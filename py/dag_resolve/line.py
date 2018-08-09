from dag_resolve import graph


class LineSegment:
    def __init__(self, line, pos, edge, seq, reads):
        self.line = line
        self.pos = pos
        self.edge = edge
        self.seq = seq
        self.reads = reads

class Line:
    def __init__(self, edge):
        assert edge.info.unique
        self.chain = [LineSegment(self, 0, edge, edge.seq, edge.reads)]

    def extendRight(self, edge, seq, reads):
        self.chain.append(LineSegment(self, self.chain[-1].pos + 1, edge, seq, reads))

    def extendLeft(self, edge, seq, reads):
        self.chain.insert(0, LineSegment(self, self.chain[0].pos - 1, edge, seq, reads))

    def rightSegment(self):
        return self.chain[-1]

    def leftSegment(self):
        return self.chain[0]

    def find(self, edge):
        for seg in self.chain:
            if seg.edge.id == edge.id:
                yield seg

    def __getitem__(self, item):
        return self.chain[item - self.chain[0].pos]

class LineStorage:
    def __init__(self, g):
        self.g = g
        self.lines = []
        self.resolved_edges = dict()
        for edge in g.E.values():
            if edge.info.unique:
                self.addLine(Line(edge))

    def addLine(self, line):
        for seg in line.chain:
            self.addSegment(seg)
        self.lines.append(line)

    def addSegment(self, seg):
        if seg.edge.id not in self.segments:
            self.resolved_edges = []
        self.resolved_edges[seg.edge.id].append(seg)

    def resolvableVerticesLeft(self):
        for v in self.g.V:
            resolved = True
            for e in v.inc:
                if e not in self.resolved_edges:
                    resolved = False
                    break
            if resolved:
                yield v

    def extendRight(self, line, edge, seq, reads):
        line.extendRight(edge, seq, reads)
        self.addSegment(line.rightSegment())





