import itertools
import sys

from alignment import align_tools
from common import basic
from dag_resolve import repeat_graph, line_tools, sequences, filters
from dag_resolve.line_tools import LineTail


class VertexResolver:
    def __init__(self, graph, polisher):
        # type: (repeat_graph.Graph, align_tools.Polisher) -> VertexResolver
        self.graph = graph
        self.lineStorage = None
        self.polisher = polisher

    def transitionSupport(self, prev, seg, next):
        # type: (repeat_graph.Edge, line_tools.LineSegment, repeat_graph.Edge) -> sequences.ReadCollection
        return seg.reads.filter(filters.EdgeTransitionFilter(prev, next))

    def resolveVertex(self, v):
        # type: (repeat_graph.Vertex) -> list[LineTail]
        print "Resolving vertex", v.id
        print "Incoming edges:", map(lambda e: e.id, v.inc)
        print "Outgoing edges:", map(lambda e: e.id, v.out)
        assert len(v.out) > 0, "Trying to resolve hanging vertex"
        res = []
        for edge in v.inc:
            print "Edge:", edge.id
            for lineSegment in self.lineStorage.resolved_edges[edge.id]:
                print "New line"
                support = []
                for next_candidate in v.out:
                    support.append(self.transitionSupport(edge, lineSegment, next_candidate))
                if len(v.out) == 1:
                    next = 0
                else:
                    supportWeight = map(len, support)
                    print "Number of edges that support transitions for the line from edge", edge.id, ":", ",".join(
                        map(str, supportWeight))
                    next, alternative = basic.best2(supportWeight, lambda x, y: x > y)
                    if supportWeight[alternative] * 5 > supportWeight[next]:
                        print "Support of alternative edge is too high. Aborting."
                        return None
                reads = support[next]# type: sequences.ReadCollection
                tail = self.polisher.polishAndAnalyse(reads, v.out[next])
                res.append(LineTail(lineSegment.line, v.out[next], tail, reads))
                print "Created tail on edge", edge.id, "of length", len(tail.cut())
        return res


class GraphResolver:
    def __init__(self, graph, vertexResolver, edgeResolver):
        # type: (repeat_graph.Graph, VertexResolver, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.lineStorage = line_tools.LineStorage(graph)
        self.vertexResolver = vertexResolver
        vertexResolver.lineStorage = self.lineStorage
        self.edgeResolver = edgeResolver
        edgeResolver.lineStorage = self.lineStorage

    def resolveVertexForward(self, v):
        # type: (repeat_graph.Vertex) -> None
        tails = sorted(self.vertexResolver.resolveVertex(v), key = lambda tail: tail.edge.id)
        for edge, tails_generator in itertools.groupby(tails, lambda tail: tail.edge):
            tails_list = list(tails_generator) # type: list[LineTail]
            if edge.id in self.lineStorage.resolved_edges:
                assert len(tails_list) == 1
                nextLine = self.lineStorage.resolved_edges[edge.id][0].line
                assert nextLine.leftSegment().edge.id == edge.id
                tails_list[0].line.merge(nextLine)
                print "Successfully connected line", tails_list[0].line.shortStr()
            else:
                while edge is not None:
                    finished, new_edge = self.edgeResolver.resolveEdge(edge, tails_list)
                    if finished:
                        print "Successfully resolved edge", edge.id, "into", len(tails_list), "lines"
                    edge = new_edge
                    if not finished:
                        print "Failed to resolve edge", edge.id
                        if edge is not None:
                            print "Graph modification detected. Trying to resolve again with edge:", new_edge.id, "of length", len(edge)
                            self.graph.printToFile(sys.stdout)
                            edge = new_edge

    def resolve(self):
        print self.graph.V
        visited_vertices = set()
        while True:
            cnt = 0
            for v in self.graph.V.values():
                v_id = v.id
                if v_id in [self.graph.source.id, self.graph.sink.id] or v_id in visited_vertices:
                    continue
                if self.lineStorage.isResolvableLeft(v):
                    self.resolveVertexForward(v)
                    visited_vertices.add(v.id)
                    cnt += 1
            if cnt == 0:
                break

    def printResults(self, handler):
        # type: (file) -> None
        nonterminal = set()
        for line in self.lineStorage.lines:
            if line.nextLine is not None:
                nonterminal.add(line.nextLine.id)
        printed = set()
        for line in self.lineStorage.lines:
            if line.id in nonterminal:
                continue
            tmp = line
            printed.add(tmp.id)
            handler.write(tmp.shortStr())
            while tmp.nextLine is not None:
                tmp = tmp.nextLine
                printed.add(tmp.id)
                handler.write("->")
                handler.write(tmp.shortStr())
            handler.write("\n")
        for line in self.lineStorage.lines:
            if line.id in printed:
                continue
            tmp = line
            printed.add(tmp.id)
            handler.write("->")
            handler.write(tmp.shortStr())
            while tmp.nextLine is not None and tmp.id not in printed:
                tmp = tmp.nextLine
                printed.add(tmp.id)
                handler.write("->")
                handler.write(tmp.shortStr())
            handler.write("->")
            handler.write("\n")


