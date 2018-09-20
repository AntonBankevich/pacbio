from common import basic, sam_parser, SeqIO
from typing import Generator, Dict, Optional, Tuple
from dag_resolve.sequences import Contig, ReadCollection, ContigCollection


class EdgeInfo:
    def __init__(self, label, unique, cov, selfrc = False):
        # type: (str, bool, int, bool) -> EdgeInfo
        self.label = label
        self.unique = unique
        self.misc = []
        self.cov = cov
        self.selfrc = selfrc


class Vertex:
    def __init__(self, id, label = None):
        # type: (int, basestring) -> Vertex
        self.id = id
        self.inc = [] # type: list[Edge]
        self.out = [] # type: list[Edge]
        self.label = label
        self.rc = None # type: Vertex

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

class Edge(Contig):
    def __init__(self, id, start, end, consensus, info):
        # type: (int, Vertex, Vertex, str, EdgeInfo) -> Edge
        Contig.__init__(self, consensus, id, info)
        self.start = start
        self.end = end
        self.reads = ReadCollection(ContigCollection([self]))
        self.rc = None # type: Edge

    def __eq__(self, other):
        # type: (Edge) -> bool
        return self.id == other.id

    def __ne__(self, other):
        return not self.__eq__(other)

class Graph:
    def __init__(self):
        # type: () -> Graph
        self.V = dict() # type: Dict[int, Vertex]
        self.E = dict() # type: Dict[int, Edge]
        self.source = self.addVertex(10000, label = "source") # type: Vertex
        self.sink = self.source.rc # type: Vertex
        self.sink.label = "sink"
        self.min_new_vid = 1
        self.min_new_eid = 6000
        self.newEdges = []
        self.reads = ReadCollection(self.edgeCollection())

    def addVertex(self, v_id = None, selfrc = False, label = None):
        # type: (Optional[int], bool, str) -> Vertex
        if v_id is None:
            v_id = self.min_new_vid
            self.min_new_vid += 1
        if v_id in self.V:
            return self.V[v_id]
        vertex = Vertex(v_id, label)
        self.V[v_id] = vertex
        if selfrc:
            vertex.rc = vertex
            return vertex
        vertex_rc = Vertex(-v_id, label)
        self.V[-v_id] = vertex_rc
        vertex.rc = vertex_rc
        vertex_rc.rc = vertex
        return vertex

    def printToFile(self, handler):
        # type: (file) -> None
        handler.write("Graph:\n")
        for edge in self.E.values():
            handler.write(str(edge.id) + "(" + str(len(edge)) + ")" + ": " + str(edge.start.id) + " -> " + str(edge.end.id) + "\n")

    def edgeCollection(self):
        return ContigCollection(self.E.values())

    def addEdge(self, edge_id, start_id, end_id, consensus, info):
        # type: (Optional[int], int, int, str, EdgeInfo) -> Edge
        if edge_id is None:
            edge_id = self.min_new_eid
            self.min_new_eid += 1
        if edge_id in self.E:
            return self.E[edge_id]
        start = self.addVertex(start_id)
        end = self.addVertex(end_id)
        edge = Edge(edge_id, start, end, consensus, info)
        edge_rc = Edge(-edge_id, end.rc, start.rc, basic.RC(consensus), info)
        edge.rc = edge_rc
        edge_rc.rc = edge
        self.E[edge.id] = edge
        self.E[edge_rc.id] = edge_rc
        start.out.append(edge)
        start.rc.inc.append(edge_rc)
        end.inc.append(edge)
        end.rc.out.append(edge_rc)
        self.newEdges.append(edge)
        self.newEdges.append(edge.rc)
        return edge

    def splitEdge(self, edge, pos_list):
        # type: (Edge, list[int]) -> list
        print "Splitting edge", edge.id, "at positions", pos_list
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
        assert len(groups) >= 2
        if len(groups) == 2:
            for pos, ind in groups[0]:
                res[ind] = edge.start
            for pos, ind in groups[1]:
                res[ind] = edge.end
        else:
            self.removeEdge(edge)
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
                    self.addEdge(None, new_vertices[i - 1].id, new_vertices[i].id, edge.seq[vertex_positions[i-1]:vertex_positions[i]], edge.info)
        return res[:-2]

    def addCuttingEdge(self, edge1, pos1, edge2, pos2, seq):
        # type: (Edge, int, Edge, int, str) -> Edge
        print "Adding cutting edge from position", pos1,"on edge", edge1.id,"to position", pos2, "on edge", edge2.id, "of length", len(seq)
        if edge1.id == edge2.id:
            vertices = self.splitEdge(edge1, [pos1, pos2])
        else:
            vertices = self.splitEdge(edge1, [pos1])
            vertices.extend(self.splitEdge(edge2, [pos2]))
        return self.addEdge(None, vertices[0].id, vertices[1].id, seq, EdgeInfo("", False))


    def removeEdge(self, edge):
        # type: (Edge) -> None
        for read in edge.reads:
            read.removeContig(edge)
            read.removeContig(edge.rc)
        edge.start.removeEdge(edge)
        edge.end.removeEdge(edge)
        edge.rc.start.removeEdge(edge.rc)
        edge.rc.end.removeEdge(edge.rc)
        del self.E[edge.id]
        del self.E[edge.rc.id]

    def loadFromDot(self, contigs, dot):
        # type: (ContigCollection, Generator[Tuple[int, int, int, int, EdgeInfo]]) -> Graph
        recs = dict() # type: Dict[int, Tuple[int, int, int, int, EdgeInfo]]
        for rec in dot:
            recs[rec[0]] = rec
        v_rc = dict()
        for eid in recs:
            if -eid in recs:
                rec1 = recs[eid]
                rec2 = recs[-eid]
                v_rc[rec1[1]] = rec2[2]
                v_rc[rec1[2]] = rec2[1]
                v_rc[rec2[1]] = rec1[2]
                v_rc[rec2[2]] = rec1[1]
        v_map = dict()
        for v in v_rc:
            if v < v_rc[v] and not v in ["source", "sink"]:
                v_map[v] = self.addVertex()
                v_map[v_rc[v]] = v_map[v].rc
            elif v == v_rc[v]:
                v_map[v] = self.addVertex(None, True)
        v_map["source"] = self.source
        v_map["sink"] = self.sink
        for eid, start, end, l, info in recs.values():
            if start not in v_map:
                v_map[start] = self.addVertex()
            if end not in v_map:
                v_map[end] = self.addVertex()
            seq = contigs[abs(eid)].seq
            if eid < 0:
                seq = basic.RC(seq)
            if info.selfrc:
                end = 40000 + abs(eid)
                v_map[end] = self.addVertex()
                seq = seq[:len(seq) / 2]
            self.addEdge(eid, v_map[start].id, v_map[end].id, seq, info)
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

    def unorientedEdges(self):
        for edge in self.E.values():
            if edge.id <= edge.rc.id:
                yield (edge, edge.rc)

class DotParser:
    def __init__(self, dot):
        # type: (file) -> DotParser
        self.dot = dot

    def parse(self, edge_ids = None):
        # type: (Optional[Dict[int, list[str]]]) -> Generator[Tuple[int, int, int, int, EdgeInfo]]
        for s in self.dot.readlines():
            if s.find("->") == -1:
                continue
            v_from = basic.parseNumber(s)
            v_to = basic.parseNumber(s, s.find("->"))
            eid = basic.parseNegativeNumber(s, s.find("id"))
            cov = basic.parseNumber(s, s.find("k "))
            # l = basic.parseNumber(s, s.find("\\l"))
            unique = (s.find("black") != -1)
            src = (s.find("dir = both") != -1)
            if edge_ids is None or eid in edge_ids:
                if edge_ids is not None:
                    if "sink" in edge_ids[eid]:
                        v_to = "sink"
                    if "source" in edge_ids[eid]:
                        v_from = "source"
                yield eid, v_from, v_to, cov, EdgeInfo(s, unique, cov, src)


