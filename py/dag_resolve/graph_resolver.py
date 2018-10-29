import itertools
import sys

import os
from typing import Optional

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import sam_parser, basic
from dag_resolve.edge_resolver import EdgeResolver
from dag_resolve.knots import Knotter
from dag_resolve.line_tools import LineStorage
from dag_resolve.repeat_graph import Graph, Edge, Vertex
from dag_resolve.sequences import Contig, ReadCollection, Segment, ContigCollection
from dag_resolve.visualization import DotPrinter, FilterColoring

class GraphResolver:
    def __init__(self, graph, dir, lineStorage, edgeResolver):
        # type: (Graph, str, LineStorage, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.dir = dir
        self.lineStorage = lineStorage
        self.edgeResolver = edgeResolver
        self.printer = DotPrinter(self.graph)
        self.printer.edge_colorings.append(FilterColoring(lambda e: e.info.unique and self.lineStorage.getLine(e.id) is not None and self.lineStorage.getLine(e.id).knot is None, "black"))
        self.printer.edge_colorings.append(FilterColoring(lambda e: e.info.unique and self.lineStorage.getLine(e.id) is not None and self.lineStorage.getLine(e.id).knot is not None, "brown"))
        self.printer.edge_colorings.append(FilterColoring(lambda e: not e.info.unique and e.id in self.lineStorage.resolved_edges, "green"))
        self.printer.edge_colorings.append(FilterColoring(lambda e: not e.info.unique and e.id not in self.lineStorage.resolved_edges, "blue"))
        # self.printer.edge_colorings.append(FilterColoring(lambda e: abs(e.id) >= 6000, "red"))
        self.cur_picture = 1

    def printCurrentGraph(self, vertices, edges, message = ""):
        fn = os.path.join(self.dir, str(self.cur_picture) + "_" + "_".join(message.split()) + ".dot")
        f = open(fn, "w")
        self.printer.printToFile(f, [FilterColoring(lambda e: e in edges, "purple")], [FilterColoring(lambda v: v in vertices, "purple")])
        f.close()
        self.cur_picture += 1


    def resolveVertexForward(self, v):
        # type: (Vertex) -> None
        self.printCurrentGraph([v], [], "Resolving vertex " + str(v.id))
        print "Resolving vertex", v.id, "inc:", map(str, v.inc), "out:", map(str, v.out)
        for edge in v.out:
            lines_list = self.lineStorage.getEdgeLines(edge)
            if len(lines_list) == 0:
                print "No line entered edge " + str(edge.id) + ". Skipping vertex."
        for edge in v.out:
            print "Considering edge", edge
            if edge.id in self.lineStorage.resolved_edges:
                print "Encountered resolved edge. Skipping."
            elif edge.seq == basic.RC(edge.seq):
                print "Encountered self-rc edge. Skipping"
            else:
                lines_list = self.lineStorage.getEdgeLines(edge)
                while edge is not None:
                    finished, new_edge = self.edgeResolver.resolveEdge(edge, lines_list)
                    if finished:
                        self.lineStorage.resolved_edges.add(edge.id)
                        self.lineStorage.resolved_edges.add(edge.rc.id)
                        print "Successfully resolved edge", edge.id, "into", len(lines_list), "lines"
                        self.lineStorage.resolved_edges.add(edge.id)
                        for line in lines_list:
                            self.lineStorage.edgeLines[line.chain[-1].seg_to.contig.id].append(line)
                        self.printCurrentGraph([], [edge], "Resolved edge " + str(edge.id))
                    else:
                        print "Failed to resolve edge", edge.id
                        if new_edge is not None:
                            print "Graph modification detected. Trying to resolve again with edge:", new_edge.id, "of length", len(new_edge)
                            self.printCurrentGraph([], [edge])
                            self.graph.printToFile(sys.stdout)
                    edge = new_edge

    def resolve(self):
        self.graph.printToFile(sys.stdout)
        visited_vertices = set()
        print "Resolving unique edges"
        for edge in self.graph.E.values():
            if edge.info.unique:# and (len(edge.end.out) != 1 or len(edge.end.inc) != 1):
                res = self.edgeResolver.processUniqueEdge(edge, self.lineStorage.getEdgeLines(edge)[0])
                assert res is not None or len(edge.end.out) == 0
                if res is not None:
                    self.lineStorage.edgeLines[res.id].append(self.lineStorage.getEdgeLines(edge)[0])
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
        Knotter(self.lineStorage).knotGraph()
        self.printCurrentGraph([], [])

    # def printResults(self, handler):
    #     # type: (file) -> None
    #     nonterminal = set()
    #     for line in self.lineStorage.lines:
    #         if line.nextLine is not None:
    #             nonterminal.add(line.nextLine.id)
    #     printed = set()
    #     for line in self.lineStorage.lines:
    #         if line.id in nonterminal:
    #             continue
    #         tmp = line
    #         printed.add(tmp.id)
    #         handler.write(tmp.__str__())
    #         while tmp.nextLine is not None:
    #             tmp = tmp.nextLine
    #             printed.add(tmp.id)
    #             handler.write("->")
    #             handler.write(tmp.__str__())
    #         handler.write("\n")
    #     for line in self.lineStorage.lines:
    #         if line.id in printed:
    #             continue
    #         tmp = line
    #         printed.add(tmp.id)
    #         handler.write("->")
    #         handler.write(tmp.__str__())
    #         while tmp.nextLine is not None and tmp.id not in printed:
    #             tmp = tmp.nextLine
    #             printed.add(tmp.id)
    #             handler.write("->")
    #             handler.write(tmp.__str__())
    #         handler.write("->")
    #         handler.write("\n")

    def printResults(self, handler):
        printed = set()
        print "Linear line chains:"
        for line in self.lineStorage.lines:
            if line.id in printed or line.rc.knot is not None:
                continue
            while line is not None:
                assert line.id not in printed
                printed.add(line.id)
                printed.add(line.rc.id)
                handler.write(line.__str__())
                if line.knot is None:
                    line = None
                else:
                    handler.write("->")
                    line = line.knot.line2
            handler.write("\n")
        print "Circular line chains:"
        for line in self.lineStorage.lines:
            if line.id in printed:
                continue
            handler.write("->")
            while line is not None and line.id not in printed:
                printed.add(line.id)
                printed.add(line.rc.id)
                handler.write(line.__str__())
                if line.knot is None:
                    line = None
                else:
                    handler.write("->")
                    line = line.knot.line2
            handler.write("->\n")



