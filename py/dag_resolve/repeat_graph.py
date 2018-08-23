from common import basic, sam_parser, SeqIO
from dag_resolve import sequences
from typing import Generator

class EdgeInfo:
    def __init__(self, label, unique):
        # type: (basestring, bool) -> EdgeInfo
        self.label = label
        self.unique = unique
        self.misc = []


class Vertex:
    def __init__(self, id, label = None):
        # type: (int, basestring) -> Vertex
        self.id = id
        self.inc = [] # type: list[Edge]
        self.out = [] # type: list[Edge]
        self.label = label
    def __eq__(self, other):
        # type: (Vertex) -> bool
        return self.id == other.id

    def __ne__(self, other):
        return not self.__eq__(other)

class Edge(sequences.Contig):
    def __init__(self, id, start, end, consensus, info = None):
        # type: (int, Vertex, Vertex, basestring, EdgeInfo) -> Edge
        sequences.Contig.__init__(self, consensus, id, info)
        self.start = start
        self.end = end
        self.reads = sequences.ReadCollection(sequences.ContigCollection([self]))

    def __eq__(self, other):
        # type: (Edge) -> bool
        return self.id == other.id

    def __ne__(self, other):
        return not self.__eq__(other)

# class VertexPort:
#     def __init__(self, v):
#         # type: (Vertex) -> object
#         self.v = v
#         self.adjacent = []
#
# class Orientation:
#     def __init__(self, reverse):
#         # type: (object) -> object
#         self.reverse = reverse
# class OrientedEdge:
#     def __init__(self, edge, orientation):

class Graph:
    def __init__(self):
        # type: () -> Graph
        self.V = dict() # type: dict[int, Vertex]
        self.E = dict() # type: dict[int, Edge]
        self.source = self.addVertex(10000, "source") # type: Vertex
        self.sink = self.addVertex(10001, "sink") # type: Vertex

    def addVertex(self, v_id, label = None):
        # type: (int, str) -> Vertex
        vertex = Vertex(v_id, label)
        self.V[v_id] = vertex
        return vertex

    def addEdge(self, edge_id, start_id, end_id, consensus, info = None):
        # type: (int, int, int, basestring, EdgeInfo) -> Edge
        if start_id not in self.V:
            self.addVertex(start_id)
        if end_id not in self.V:
            self.addVertex(end_id)
        start = self.V[start_id]
        end = self.V[end_id]
        edge = Edge(edge_id, start, end, consensus, info)
        self.E[edge.id] = edge
        start.out.append(edge)
        end.inc.append(edge)
        return edge

    def loadFromDot(self, contigs, dot):
        # type: (sequences.ContigCollection, Generator[tuple]) -> Graph
        for eid, start, end, l, info in dot:
            if start == "source":
                start = self.source.id
            if start == "sink":
                start = self.sink.id
            seq = contigs[abs(eid)].seq
            if eid < 0:
                seq = basic.RC(seq)
            self.addEdge(eid, start, end, seq, info)
        return self

    def fillAlignments(self, read_recs, sam, fill_unique = True):
        # type: (Generator[SeqIO.SeqRecord], sam_parser.Samfile) -> None
        reads = sequences.ReadCollection(sequences.ContigCollection(self.E.values()))
        for rec in read_recs:
            reads.addNewRead(rec)
        for rec in sam:
            if rec.is_unmapped:
                continue
            edge_id = int(rec.tname)
            if not fill_unique and rec.pos > 5000 and rec.pos + rec.alen + 5000 <= self.E[edge_id].__len__():
                continue
            self.E[edge_id].reads.add(reads[rec.query_name])
            self.E[edge_id].reads.addNewAlignment(rec)



class DotParser:
    def __init__(self, dot):
        # type: (file) -> DotParser
        self.dot = dot

    def parse(self, edge_ids = None):
        # type: (dict[int, list[str]]) -> Generator[tuple]
        for s in self.dot.readlines():
            if s.find("->") == -1:
                continue
            v_from = basic.parseNumber(s)
            v_to = basic.parseNumber(s, s.find("->"))
            eid = basic.parseNegativeNumber(s, s.find("id"))
            l = basic.parseNumber(s, s.find("\\l"))
            unique = (s.find("black") != -1)
            if edge_ids is None or eid in edge_ids:
                if "sink" in edge_ids[eid]:
                    yield eid, v_from, -1, l, EdgeInfo(s, unique)
                if "source" in edge_ids[eid]:
                    yield eid, v_from, -2, l, EdgeInfo(s, unique)
                yield eid, v_from, v_to, l, EdgeInfo(s, unique)


