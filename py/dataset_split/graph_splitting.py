import sys
import os
from typing import List, Callable, Generator

sys.path.append("py")

from common.SimpleGraph import SimpleGraph, Edge
from common import basic

class DipolidCalculator:
    def __init__(self, edge_length):
        self.edge_length = edge_length

    def calculateComponentCoverage(self, comp, max_cov):
        low = []
        high = []
        for vid in comp.v:
            for e in comp.v[vid].out:
                if e.len >= self.edge_length:
                    if e.cov < max_cov * 3 / 4 and e.cov >= max_cov / 4:
                        low.append(e)
                    elif e.cov < max_cov * 5/4 and e.cov >= max_cov * 3 / 4:
                        high.append(e)
        if len(low) == 0 and len(high) == 0:
            cov = min(comp.covPerc(0.5, max_cov * 1.5), max_cov)
        elif len(low) == 0:
            cov = sum([edge.cov for edge in high]) / len(high)
        elif len(high) == 0:
            cov = sum([edge.cov for edge in low]) / len(low)
        else:
            cov = sum([edge.cov for edge in low]) / len(low)
        return cov

    def uniqueCondition(self, cov):
        return lambda edge: edge.cov < cov * 1.4 and edge.cov >= cov * 0.7

    def edgeColoring(self, cov):
        return lambda edge: "blue" if edge.cov < cov * 0.7 else ("black" if edge.cov < cov * 1.4 else "red")

    def isRepeat(self, edge, cov):
        return edge.cov >= cov * 1.4

class HaploidCalculator:
    def __init__(self, edge_length):
        self.edge_length = edge_length

    def calculateComponentCoverage(self, comp, max_cov):
        edges = []
        for vid in comp.v:
            for e in comp.v[vid].out:
                if e.len >= self.edge_length:
                    if e.cov > max_cov * 0.5 and e.cov < max_cov * 1.4:
                        edges.append(e)
        if len(edges) == 0:
            return max_cov
        else:
            return sum([edge.cov for edge in edges]) / len(edges)


    def uniqueCondition(self, cov):
        return lambda edge: edge.cov < cov * 1.4 and edge.cov >= cov * 0.5

    def edgeColoring(self, cov):
        return lambda edge: "blue" if edge.cov < cov * 0.5 else ("black" if edge.cov < cov * 1.4 else "red")

    def isRepeat(self, edge, cov):
        return edge.cov >= cov * 1.4


def SplitGraphByCondition(graph, unique_condition):
    # type: (SimpleGraph, Callable[[Edge], bool]) -> Generator[SimpleGraph]
    visited = set()
    for v in graph.v:
        if v in visited or graph.v[v].isIsolated() or not graph.v[v].core:
            continue
        queue = [v]
        comp = set()
        comp.add(v)
        comp.add(graph.rcVertex(graph.v[v]).id)
        while len(queue) > 0:
            next_id = queue.pop()
            if next_id in visited or len(graph.v[next_id].inc) + len(graph.v[next_id].out) > 50 or not graph.v[next_id].core:
                continue
            rcv = graph.rcVertex(graph.v[next_id]).id
            visited.add(next_id)
            visited.add(rcv)
            comp.add(next_id)
            comp.add(rcv)
            next = graph.v[next_id]
            for edge in next.inc + next.out:
                if not unique_condition(edge):
                    queue.append(edge.start)
                    queue.append(edge.end)
        if len(comp) > 1:
            yield graph.Component(comp)

def SplitGraph(graph, calculator):
    max_cov = graph.covPerc(0.5)
    for comp in SplitGraphByCondition(graph, lambda edge: edge.len >= calculator.edge_length and edge.cov < max_cov * 1.8):
        cov = calculator.calculateComponentCoverage(comp, max_cov)
        if cov == 0:
            print "Zero component"
        for comp1 in SplitGraphByCondition(comp, calculator.uniqueCondition(cov)):
            print comp1.v.__len__(), comp1.e.__len__()
            yield comp1, cov


if __name__ == "__main__":
    g = SimpleGraph()
    g.ReadDot(sys.argv[1])
    basic.ensure_dir_existance(sys.argv[2])
    args = sys.argv[3:]
    if "merge" in args:
        g = g.Merge()
    diploid = "--diploid" in args
    cnt = 0
    oppa = []
    simple = 0
    complex = 0
    max_cov = g.covPerc(0.5)
    if diploid:
        calculator = DipolidCalculator(150000)
    else:
        calculator = HaploidCalculator(150000)
    for comp, cov in SplitGraph(g, calculator):
        if len(comp) <= 2:
            continue
        if len(comp) <= 4:
            simple += 1
            oppa.extend(comp.v)
            continue
        if len(comp) > 100:
            print "Complex", comp.__len__()
            complex += len(comp)
            continue
        print cnt, len(comp), comp.__len__(), cov, max_cov
        f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w")
        coloring = lambda v: "white" if len(v.inc) + len(v.out) == len(g.v[v.id].inc) + len(g.v[v.id].out) else ("yellow" if len(g.v[v.id].inc) + len(g.v[v.id].out) < 50 else "red")
        comp.Print(f, coloring, calculator.edgeColoring(cov))
        f.close()
        cnt += 1
    f = open(os.path.join(sys.argv[2], "small.dot"), "w")
    g.Draw(set(oppa), f)
    f.close()
    print simple
