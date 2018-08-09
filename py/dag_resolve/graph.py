from common import basic
from dag_resolve import sequences

class EdgeInfo:
    def __init__(self, label, unique):
        self.label = label
        self.unique = unique
        self.misc = []


class Vertex:
    def __init__(self, id, label = None):
        self.id = id
        self.inc = []
        self.out = []
        self.label = label

class Edge(sequences.Contig):
    def __init__(self, id, start, end, consensus, info = None):
        sequences.Contig.__init__(self, consensus, id, info)
        self.start = start
        self.end = end
        self.reads = sequences.ReadCollection(sequences.ContigCollection([self]))

class Graph:
    def __init__(self):
        self.V = dict()
        self.E = dict()
        self.source = self.addVertex(10000, "source")
        self.sink = self.addVertex(10001, "sink")

    def addVertex(self, v_id, label = None):
        vertex = Vertex(v_id, label)
        self.V[v_id] = vertex
        return vertex

    def addEdge(self, edge_id, start_id, end_id, consensus, info = None):
        if start_id not in self.V:
            self.addVertex(start_id)
        if end_id not in self.V:
            self.addVertex(end_id)
        start = self.V[start_id]
        end = self.V[end_id]
        edge = Edge(edge_id, start, end, consensus, info)
        start.out.append(edge)
        end.inc.append(edge)

    def loadFromDot(self, contigs, dot):
        dot = DotParser("")
        contigs = sequences.ContigCollection()
        for eid, start, end, l, info in dot.parse():
            if start == "source":
                start = self.source.id
            if start == "sink":
                start = self.sink.id
            seq = contigs[abs(eid)].seq
            if eid < 0:
                seq = basic.RC(seq)
            self.addEdge(eid, start, end, seq, info)




class DotParser:
    def __init__(self, dot):
        self.dot = dot

    def parse(self):
        for s in open(self.dot, "r").readlines():
            tmp =s.strip().split()
            if len(s) < 2 or s[1] != "->":
                continue
            v_from = basic.parseNumber(s)
            v_to = basic.parseNumber(s, s.find("->"))
            eid = basic.parseNumber(s, s.find("id"))
            l = basic.parseNumber(s, s.find("\\l"))
            unique = (s.find("black") != -1)
            yield eid, v_from, v_to, l, EdgeInfo(s, unique)


