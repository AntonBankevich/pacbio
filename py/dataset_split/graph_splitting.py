import sys
import os
from typing import List, Callable, Generator

sys.path.append("py")

from common.SimpleGraph import SimpleGraph, Edge
from common import basic

def SplitGraphByCondition(graph, condition):
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
                if condition(edge):
                    queue.append(edge.start)
                    queue.append(edge.end)
        if len(comp) > 1:
            yield graph.Component(comp)

def SplitGraph(graph, edge_length = 150000):
    max_cov = graph.covPerc(0.5)
    for comp in SplitGraphByCondition(graph, lambda edge: edge.len < edge_length):
        cov = min(comp.covPerc(0.5, max_cov * 1.5), max_cov)
        low = []
        high = []
        for vid in comp.v:
            for e in comp.v[vid].out:
                if e.len >= edge_length:
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
        if cov == 0:
            print "Gopa", cov, str([edge.cov for edge in low]), str([edge.cov for edge in high])
        for comp1 in SplitGraphByCondition(comp, lambda edge: edge.cov > cov * 1.4 or edge.cov <= cov * 0.7):
            print comp1.v.__len__(), comp1.e.__len__()
            yield comp1, cov


if __name__ == "__main__":
    g = SimpleGraph()
    g.ReadDot(sys.argv[1])
    basic.ensure_dir_existance(sys.argv[2])
    args = sys.argv[3:]
    if "merge" in args:
        g = g.Merge()
    cnt = 0
    oppa = []
    simple = 0
    complex = 0
    max_cov = g.covPerc(0.5)
    print g.covPerc(0.1), g.covPerc(0.2), g.covPerc(0.5), g.covPerc(0.6), g.covPerc(0.7),g.covPerc(0.8)
    for comp, cov in SplitGraph(g, 150000):
        if len(comp) <= 1:
            continue
        if len(comp) < 3:
            simple += 1
            oppa.extend(comp.v)
            continue
        if len(comp) > 100:
            complex += len(comp)
            continue
        print cnt, len(comp), comp.__len__(), cov, max_cov
        f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w")
        coloring = lambda v: "white" if len(v.inc) + len(v.out) == len(g.v[v.id].inc) + len(g.v[v.id].out) else ("yellow" if len(g.v[v.id].inc) + len(g.v[v.id].out) < 50 else "red")
        comp.Print(f, cov, coloring)
        f.close()
        cnt += 1
    f = open(os.path.join(sys.argv[2], "small.dot"), "w")
    g.Draw(set(oppa), f)
    f.close()
    print simple
