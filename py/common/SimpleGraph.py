from typing import Dict, List, Generator

from common import basic, SeqIO
from common.disjoint_set import DisjointSet


def toint(s):
    if s.startswith("\""):
        s = s[1:]
    if s.endswith("\""):
        s = s[:-1]
    return int(s)

class Edge:
    def __init__(self, id, start, fin, len, label, cov = None, seq = None):
        # type: (str, str, str, int, str, int, str) -> None
        self.id = id # type: str
        self.seq = seq # type: str
        self.label = label # type: str
        self.start = start # type: str
        self.len = len # type: int
        self.end = fin # type: str
        self.unique = (self.label.find("black") != -1)
        self.cov = cov

    def __len__(self):
        return self.len


class SimpleGraph:
    def __init__(self, graph = None):
        if graph is None:
            self.g = self
        else:
            self.g = graph
        self.v = dict() #type: Dict[str, Vertex]
        self.visited = set()
        self.e = dict() #type: Dict[str, Edge]

    def __len__(self):
        return len(self.v)

    def isBorder(self, vid):
        thisv = self.v[vid]
        otherv = self.g.v[vid]
        return thisv.inc.__len__() + thisv.out.__len__() != otherv.inc.__len__() + otherv.out.__len__()

    def isHub(self, vid):
        otherv = self.g.v[vid]
        return otherv.inc.__len__() + otherv.out.__len__() > 50

    def AddVertex(self, vid, label = ""):
        if vid not in self.v:
            self.v[vid] = Vertex(vid, label)
        return self.v[vid]

    def AddEdge(self, id, start, fin, len, label, cov = None, seq = None):
        v1 = self.AddVertex(start)
        v2 = self.AddVertex(fin)
        edge = Edge(id, start, fin, len, label, cov, seq)
        self.e[id] = edge
        v1.out.append(edge)
        v2.inc.append(edge)
        return edge

    def rcEdge(self, edge):
        # type: (Edge) -> Edge
        if basic.isCanonocal(edge.id) and basic.Reverse(edge.id) not in self.e:
            start = self.v[edge.start]
            end = self.v[edge.end]
            assert len(start.inc) == len(end.out) and len(start.out) == len(end.inc)
            return edge
        return self.e[basic.Reverse(edge.id)]

    def isSelfRc(self, edge):
        return edge.id == self.rcEdge(edge).id

    def rcVertex(self, vertex):
        # type: (Vertex) -> Vertex
        if vertex.out.__len__() > 0:
            return self.v[self.rcEdge(vertex.out[0]).end]
        elif vertex.inc.__len__() > 0:
            return self.v[self.rcEdge(vertex.inc[0]).start]
        assert False


    def ReadDot(self, f):
        for s in open(f, "r").readlines():
            tmp = s
            if s[0] != "\"":
                continue
            s =s.strip().split()
            if len(s) < 2 or s[1] != "->":
                vid = s[0][1:-1]
                self.AddVertex(vid, tmp)
                continue
            v_from = s[0][1:-1]
            v_to = s[2][1:-1]
            label = tmp
            cov = basic.parseNumber(tmp, tmp.find("k "))
            l = int(float(tmp[tmp.find("\\l") + 2: tmp.find("k ")]) * 1000)
            id = tmp[tmp.find("id ") + 3: tmp.find("\\l")]
            edge = self.AddEdge(id, v_from, v_to, l, label, cov)
        return self

    def Merge(self):
        res = SimpleGraph()
        for v in self.v.values():
            if len(v.inc) != 1 or len(v.out) != 1:
                if v.id not in res.v:
                    res.AddVertex(v.id)
                for e in v.out:
                    eids = [e.id]
                    cov = e.cov * e.len
                    l = e.len
                    seqs = [e.seq]
                    labels = [e.label]
                    while not self.v[e.end].isJunction():
                        e = self.v[e.end].out[0]
                        eids.append(e.id)
                        cov += e.cov * e.len
                        l += e.len
                        seqs.append(e.seq)
                    endid = self.v[e.end].id
                    if endid not in res.v:
                        res.AddVertex(endid)
                    if len(eids) == 1:
                        new_id = eids[0]
                    else:
                        new_id = "(" + ",".join(eids) + ")"
                        id_cand = "(" + ",".join(map(basic.Reverse, eids[::-1])) + ")"
                        if id_cand > new_id:
                            new_id = basic.Reverse(id_cand)
                    res.AddEdge(new_id, v.id, endid, l, "".join(labels), int(float(cov) / l),
                                   None if seqs[0] is None else "".join(seqs))
        return res


    def FillSeq(self, f, numeric = True):
        for s in SeqIO.parse_fasta(open(f, "r")):
            if numeric:
                s.id = str(basic.parseNumber(s.id))
            if s.id in self.e:
                self.e[s.id].seq = s.seq
                self.e[s.id].len = len(s.seq)
            if basic.Reverse(s.id) in self.e:
                self.e[basic.Reverse(s.id)].seq = basic.RC(s.seq)
                self.e[basic.Reverse(s.id)].len = len(s.seq)
        for edge in self.e.values():
            assert(edge.seq is not None)
        return self

    def ReadGFA(self, f):
        seqs = dict()
        v = DisjointSet()
        for s in open(f, "r").readlines():
            s = s.split()
            if s[0] == "S":
                seqs[s[1]] = s[2]
                v.add((True, s[1], True))
                v.add((True, s[1],False))
                v.add((False, s[1], True))
                v.add((False, s[1],False))
            elif s[0] == "L":
                v1 = (s[2] == "+", s[1], True)
                v2 = (s[4] == "+", s[3], False)
                v.union(v1, v2)
                v1, v2 = ((not v2[0], v2[1], not v2[2]), (not v1[0], v1[1], not v1[2]))
                v.union(v1, v2)
        ids = dict()
        cnt = 1
        for vid, vl in v.listComponenets():
            self.AddVertex(str(cnt), str(vid))
            ids[vid] = str(cnt)
            cnt += 1
        for eid, seq in seqs.items():
            self.AddEdge(eid, ids[v.get((True, eid, False))], ids[v.get((True, eid, True))], len(seq), eid, seq=seq)
            self.AddEdge("-" + eid, ids[v.get((False, eid, False))], ids[v.get((False, eid, True))], len(seq), "-" + eid, seq=basic.RC(seq))
        return self

    def Find(self, minlen, v):
        self.visited.add(v)
        for e in self.v[v].out:
            if e.len <= minlen and e.end not in self.visited:
                self.Find(minlen, e.end)
        for e in self.v[v].inc:
            if e.len <= minlen and e.start not in self.visited:
                self.Find(minlen, e.start)

    def Split(self, minlen):
        # type: (int) -> Generator[List[Vertex]]
        visited = set()
        for v in self.v.values():
            if v.id in visited:
                continue
            self.visited = set()
            self.Find(minlen, v.id)
            visited = visited.union(self.visited)
            yield list(self.visited)

    def Component(self, comp):
        res = SimpleGraph(self.g)
        for vid in comp:
            vert = self.v[vid]
            res.AddVertex(vid)
            for e in vert.out: #type: Edge
                if e.end not in res.v:
                    res.AddVertex(e.end)
                res.AddEdge(e.id, e.start, e.end, e.len, e.label, e.cov, e.seq)
            for e in vert.inc:  # type: Edge
                if e.start not in comp:
                    res.AddVertex(e.start)
                    res.AddEdge(e.id, e.start, e.end, e.len, e.label, e.cov, e.seq)
        for vert in res.v.values():
            vert.core = vert.id in comp
        return res

    def Draw(self, comp, out):
        out.write("digraph {\n")
        out.write("nodesep = 0.5;\n")
        for vid in comp:
            v = self.v[vid]
            if v.label != "":
                out.write(v.label)
        for v in self.v.values():
            for e in v.out:
                if e.start in comp or e.end in comp:
                    out.write(e.label + "\n")
        out.write("}\n")

    def covPerc(self, perc, max = 1000000):
        edges = []
        size = 0
        for vert in self.v:
            for edge in self.v[vert].out:
                if edge.cov < max:
                    edges.append(edge)
                    size += edge.len
        if len(edges) == 0:
            return self.covPerc(perc)
        edges.sort(key = lambda edge: -edge.cov)
        tmp = 0
        for edge in edges:
            if tmp + edge.len >= size * perc:
                return edge.cov
            tmp += edge.len

    def Print(self, out, vertex_coloring = lambda v: "white", edge_coloring = lambda e: "black"):
        out.write("digraph {\n")
        out.write("nodesep = 0.5;\n")
        for vid in self.v:
            v = self.v[vid]
            out.write(str(vid) + " [style = filled, fillcolor = " + vertex_coloring(v) + "];\n")
        for v in self.v.values():
            for e in v.out:
                col = edge_coloring(e)
                out.write("\"" + str(e.start) + "\" -> \"" + str(e.end) + "\" [label = \"id " + str(e.id) + "\\n" +
                          str(e.len /100 * 0.1) + "k " + str(e.cov) + "x\", color = \"" + col + "\"];\n")
        out.write("}\n")

    def PrintFasta(self, out):
        for e in self.e.values():
            SeqIO.write(e, out, "fasta")


class Vertex:
    def __init__(self, id, label):
        self.id = id
        self.inc = [] # type: List[Edge]
        self.out = [] # type: List[Edge]
        self.label = label
        self.core = True

    def isIsolated(self):
        return len(self.inc) == 0 and len(self.out) == 0
    def isJunction(self):
        return len(self.inc) != 1 or len(self.out) != 1
