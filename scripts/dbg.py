import json
import sys
import os

import pickle

sys.path.append("py")
from common import sam_parser, basic
from common.basic import RC
from flye_tools.alignment import make_alignment
from flye_tools import polysh_job

import common.SeqIO as SeqIO

def ConstructEdges(reads, k, cov):
    cnt = 0
    edges = dict()
    for read in reads:
        read.seq = read.seq.upper()
        cnt += 1
        if cnt % 400 == 0:
            sys.stderr.write(str(cnt) + " " + str(len(edges)) + "\n")
        for i in range(k + 1, len(read) + 1):
            kmer = read.seq[i - k - 1:i]
            if kmer not in edges:
                edges[kmer] = 0
            edges[kmer] += 1
            kmer = RC(kmer)
            if kmer not in edges:
                edges[kmer] = 0
            edges[kmer] += 1
    res = dict()
    for kmer in edges:
        if edges[kmer] >= cov:
            res[kmer] = edges[kmer]
    sys.stderr.write(str(len(res)) + "\n")
    return res

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.next = []
        self.prev = []

    def simple(self):
        return len(self.next) == 1 and len(self.prev) == 1

class DBG:
    def __init__(self, edges, k):
        self.v = set()
        self.e = edges
        self.k = k
        for edge in edges:
            self.v.add(edge[1:])
            self.v.add(edge[:-1])

    def out(self, v):
        for c in "ACGT":
            edge = v + c
            if edge in self.e:
                yield edge[1:]

    def inc(self, v):
        for c in "ACGT":
            edge = c + v
            if edge in self.e:
                yield edge[:-1]

    def incdeg(self, v):
        return len(list(self.inc(v)))
    def outdeg(self, v):
        return len(list(self.out(v)))
    def isolated(self, v):
        return self.incdeg(v) == 0 and self.outdeg(v) == 0
    def simple(self, v):
        return self.incdeg(v) == 1 and self.outdeg(v) == 1

    def clean(self, l):
        for v in self.v:
            if self.incdeg(v) == 0 and self.outdeg(v) == 1:
                for next in self.out(v):
                    e = self.unique_forward(v[0], next)
                    vend = e[-self.k:]
                    if len(e) < l + self.k or (self.incdeg(vend) == 1 and self.outdeg(vend) == 0):
                        # sys.stderr.write(str(len(e)) + "\n")
                        for i in range(self.k + 1, len(e) + 1):
                            # sys.stderr.write(e[i - self.k - 1: i] + "\n")
                            edge = e[i - self.k - 1: i]
                            if edge in self.e:
                                del self.e[edge]
                                rc = RC(edge)
                                if rc in self.e:
                                    del self.e[rc]

    def unique_forward(self, c, v):
        res = [c]
        curv = v
        while self.simple(curv):
            res.append(curv[0])
            curv = self.out(curv).next()
        return "".join(res) + curv

    def cov(self, e):
        s = 0.
        for i in range(self.k + 1, len(e) + 1):
            s += self.e[e[i - self.k - 1: i]]
        return int(s / (len(e) - self.k))

    def ntostr(self, n):
        return "\"" + str(n) + "\""
    def print_dot(self, handler):
        handler.write("digraph {\n")
        handler.write("nodesep = 0.5;\n")
        handler.write("node [shape = circle, label = \"\", height = 0.3];\n")
        vert = dict()
        cnt = 0
        for v in self.v:
            inc = list(self.inc(v))
            out = list(self.out(v))
            if (self.incdeg(v) != 1 or self.outdeg(v) != 1) and self.incdeg(v) + self.outdeg(v) != 0:
                vert[v] = cnt
                cnt += 1
        for v in vert:
            for next in self.out(v):
                e = self.unique_forward(v[0], next)
                handler.write(self.ntostr(vert[v]) + " -> " + self.ntostr(vert[e[-self.k:]]) +  "[label = " +  "\"" + str(len(e) - self.k) + ":" + str(self.cov(e)) + "\"" + ", color = \"black\"] ;\n")
        handler.write("}\n")

sys.stderr.write("starting\n")
k = int(sys.argv[2])
edges = ConstructEdges(SeqIO.parse_fasta(open(sys.argv[1], "r")), k, int(sys.argv[3]))
sys.stderr.write("constructing dbg\n")
g = DBG(edges, k)
sys.stderr.write("cleaning")
g.clean(20)
sys.stderr.write("printing dbg\n")
g.print_dot(sys.stdout)


