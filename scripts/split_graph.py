import sys

import os


class Edge:
    def __init__(self, start, fin, len, label):
        self.label = label
        self.start = start
        self.len = len
        self.fin = fin

class Vertex:
    def __init__(self, id):
        self.id = id
        self.inc = []
        self.out = []

def toint(s):
    if s.startswith("\""):
        s = s[1:]
    if s.endswith("\""):
        s = s[:-1]
    return int(s)

class Graph:
    def __init__(self):
        self.v = dict()
        self.visited = set()

    def AddVertex(self, vid):
        if vid not in self.v:
            self.v[vid] = Vertex(vid)
        return self.v[vid]

    def AddEdge(self, start, fin, len, label):
        v1 = self.AddVertex(start)
        v2 = self.AddVertex(fin)
        edge = Edge(start, fin, len, label)
        v1.out.append(edge)
        v2.inc.append(edge)


    def Read(self, f):
        for s in open(f, "r").readlines():
            tmp = s
            s =s.strip().split()
            if len(s) < 2 or s[1] != "->":
                continue
            v_from = toint(s[0])
            v_to = toint(s[2])
            label = tmp
            l = float(tmp[tmp.find("\\l") + 2: tmp.find("k ")])
            self.AddEdge(v_from, v_to, l, label)

    def Find(self, minlen, v):
        if v == 18 or v == -18:
            return
        self.visited.add(v)
        for e in self.v[v].out:
            if e.len <= minlen and e.fin not in self.visited:
                self.Find(minlen, e.fin)
        for e in self.v[v].inc:
            if e.len <= minlen and e.start not in self.visited:
                self.Find(minlen, e.start)

    def Split(self, minlen):
        visited = set()
        for v in self.v.values():
            if v.id in visited:
                continue
            self.visited = set()
            self.Find(minlen, v.id)
            visited = visited.union(self.visited)
            yield self.visited

    def Draw(self, comp, out):
        out.write("digraph {\n")
        out.write("nodesep = 0.5;\n")
        for v in self.v.values():
            for e in v.out:
                if e.start in comp or e.fin in comp:
                    out.write(e.label + "\n")
        out.write("}\n")


g = Graph()
g.Read(sys.argv[1])
os.makedirs(sys.argv[2])
cnt = 1
for comp in g.Split(50):
    print len(comp)
    f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w   ")
    g.Draw(comp, f)
    f.close()
    cnt += 1





