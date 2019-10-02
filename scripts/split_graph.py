import sys

import os
sys.path.append("py")
from common import basic


class Edge:
    def __init__(self, start, fin, len, label):
        self.label = label
        self.start = start
        self.len = len
        self.fin = fin

class Vertex:
    def __init__(self, id, label):
        self.id = id
        self.inc = []
        self.out = []
        self.label = label

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

    def AddVertex(self, vid, label = ""):
        if vid not in self.v:
            self.v[vid] = Vertex(vid, label)
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
            if s[0] != "\"":
                continue
            s =s.strip().split()
            if len(s) < 2 or s[1] != "->":
                vid = toint(s[0])
                self.AddVertex(vid, tmp)
                continue
            v_from = toint(s[0])
            v_to = toint(s[2])
            label = tmp
            l = float(tmp[tmp.find("\\l") + 2: tmp.find("k ")])
            self.AddEdge(v_from, v_to, l, label)

    def Find(self, minlen, v):
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
        for vid in comp:
            v = self.v[vid]
            if v.label != "":
                out.write(v.label)
        for v in self.v.values():
            for e in v.out:
                if e.start in comp or e.fin in comp:
                    out.write(e.label + "\n")
        out.write("}\n")


g = Graph()
g.Read(sys.argv[1])
basic.ensure_dir_existance(sys.argv[2])
for cnt, comp in enumerate(g.Split(50)):
    if len(comp) < 3:
        continue
    print len(comp)
    f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w   ")
    g.Draw(comp, f)
    f.close()
    cnt += 1





