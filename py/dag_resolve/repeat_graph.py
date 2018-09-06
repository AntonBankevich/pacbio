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

    def removeEdge(self, edge):
        self.inc = filter(lambda e: e.id != edge.id, self.inc)
        self.out = filter(lambda e: e.id != edge.id, self.out)

    def __eq__(self, other):
        # type: (Vertex) -> bool
        return self.id == other.id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return str(self.id)

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
        self.min_new_vid = 5000
        self.min_new_eid = 6000
        self.newEdges = []
        self.reads = sequences.ReadCollection(self.edgeCollection())

    def addVertex(self, v_id = None, label = None):
        # type: (int, str) -> Vertex
        if v_id is None:
            v_id = self.min_new_vid
            self.min_new_vid += 1
        vertex = Vertex(v_id, label)
        self.V[v_id] = vertex
        return vertex

    def printToFile(self, handler):
        # type: (file) -> None
        handler.write("Graph:\n")
        for edge in self.E.values():
            handler.write(str(edge.id) + ": " + str(edge.start.id) + " -> " + str(edge.end.id) + "\n")

    def edgeCollection(self):
        return sequences.ContigCollection(self.E.values())

    def addEdge(self, edge_id, start_id, end_id, consensus, info = None):
        # type: (int, int, int, str, EdgeInfo) -> Edge
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
        self.newEdges.append(edge)
        return edge

    def addNewEdge(self, start_id, end_id, consensus, info = EdgeInfo("new", False)):
        # type: (int, int, str, EdgeInfo) -> Edge
        eid = self.min_new_eid
        self.min_new_eid += 1
        return self.addEdge(eid, start_id, end_id, consensus, info)

    def splitEdge(self, edge, pos_list):
        # type: (Edge, list[int]) -> list
        print "Splitting edge", edge.id, "at positions", pos_list
        self.removeEdge(edge)

        pos_list.append(0)
        pos_list.append(len(edge))
        pos_ind = sorted([(pos, i) for i, pos in enumerate(pos_list)])
        groups = [[pos_ind[0]]]
        prev = pos_ind[0][0]
        for pos, ind in pos_ind[1:]:
            if pos - prev > 500:
                groups.append([])
            groups[-1].append((pos, ind))
            prev = pos
        res = [None] * (len(pos_list))
        new_vertices = []
        vertex_positions = []
        for i, group in enumerate(groups):
            left = min([pos for pos, ind in group])
            right = max([pos for pos, ind in group])
            if left == 0:
                new_vertices.append(edge.start)
                vertex_positions.append(0)
            elif right == len(edge):
                new_vertices.append(edge.end)
                vertex_positions.append(len(edge))
            else:
                vertex_positions.append((left + right) // 2)
                new_vertices.append(self.addVertex())
            for pos, ind in group:
                res[ind] = new_vertices[-1]
            if i > 0:
                self.addNewEdge(new_vertices[i - 1].id, new_vertices[i].id, edge.seq[vertex_positions[i-1]:vertex_positions[i]])
        return res[:-2]

    def addCuttingEdge(self, edge1, pos1, edge2, pos2, seq):
        # type: (Edge, int, Edge, int, str) -> Edge
        print "Edding cutting edge from position", pos1,"on edge", edge1.id,"to position", pos2, "on edge", edge2.id, "of length", len(seq)
        if edge1.id == edge2.id:
            vertices = self.splitEdge(edge1, [pos1, pos2])
        else:
            vertices = self.splitEdge(edge1, [pos1])
            vertices.extend(self.splitEdge(edge2, [pos2]))
        return self.addNewEdge(vertices[0].id, vertices[1].id, seq)


    def removeEdge(self, edge):
        # type: (Edge) -> None
        for read in edge.reads:
            read.removeContig(edge)
        edge.start.removeEdge(edge)
        edge.end.removeEdge(edge)
        del self.E[edge.id]

    def loadFromDot(self, contigs, dot):
        # type: (sequences.ContigCollection, Generator[tuple]) -> Graph
        for eid, start, end, l, info in dot:
            if start == "source":
                start = self.source.id
            if end == "sink":
                end = self.sink.id
            seq = contigs[abs(eid)].seq
            if eid < 0:
                seq = basic.RC(seq)
            self.addEdge(eid, start, end, seq, info)
        return self

    def fillAlignments(self, read_recs, sam, fill_unique = True):
        # type: (Generator[SeqIO.SeqRecord], sam_parser.Samfile) -> None
        for rec in read_recs:
            self.reads.addNewRead(rec)
        for rec in sam:
            if rec.is_unmapped:
                continue
            edge_id = int(rec.tname)
            if not fill_unique and self.E[edge_id].info.unique and rec.pos > 5000 and rec.pos + rec.alen + 5000 <= self.E[edge_id].__len__():
                continue
            self.E[edge_id].reads.add(self.reads[rec.query_name])
            self.E[edge_id].reads.addNewAlignment(rec)
        self.newEdges = []



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
                if edge_ids is not None:
                    if "sink" in edge_ids[eid]:
                        v_to = "sink"
                    if "source" in edge_ids[eid]:
                        v_from = "source"
                yield eid, v_from, v_to, l, EdgeInfo(s, unique)


