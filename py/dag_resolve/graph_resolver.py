import itertools
import sys

import os

from common import basic
from dag_resolve.edge_resolver import EdgeResolver
from dag_resolve.knots import Knotter
from dag_resolve.line_tools import LineStorage
from dag_resolve.repeat_graph import Graph, Vertex
from dag_resolve.visualization import DotPrinter, FilterColoring, HistoryPrinter


class GraphResolver:
    def __init__(self, graph, printer, lineStorage, edgeResolver):
        # type: (Graph, HistoryPrinter, LineStorage, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.lineStorage = lineStorage
        self.edgeResolver = edgeResolver
        self.printer = printer
        # self.printer.edge_colorings.append(FilterColoring(lambda e: abs(e.id) >= 6000, "red"))

    def resolveVertexForward(self, v):
        # type: (Vertex) -> None
        self.printer.printCurrentGraph([v], [], "Resolving vertex " + str(v.id))
        print "Resolving vertex", v.id, "inc:", map(str, v.inc), "out:", map(str, v.out)
        lines = list(itertools.chain(*[self.lineStorage.edgeLines[edge.id] for edge in v.inc]))
        resolved = self.edgeResolver.resolveVertex(v, lines)
        for edge, line in zip(resolved, lines):
            if edge is not None:
                self.lineStorage.edgeLines[edge.id].append(line)
        if None not in resolved:
            print "Successfully resolved vertex", v
            self.printer.printCurrentGraph([v], [], "Successfully resolved vertex " + str(v))
        else:
            print "Failed to resolve vertex", v
            print "WARNING: VERTEX WAS NOT RESOLVED"
            print map(str, zip(map(str, resolved), map(str, lines)))
            self.printer.printCurrentGraph([v], [], "Failed to resolve vertex " + str(v))
            return
        for edge in v.out:
            lines_list = self.lineStorage.getEdgeLines(edge)
            if len(lines_list) == 0:
                print "No line entered edge " + str(edge.id) + ". Skipping vertex."
        for edge in v.out:
            print "Considering edge", edge
            if edge.rc.id in self.lineStorage.resolved_edges:
                print "Encountered edge resolved from opposite side. Skipping."
            elif edge.seq == basic.RC(edge.seq):
                print "Encountered self-rc edge. Skipping"
            else:
                lines_list = self.lineStorage.getEdgeLines(edge)
                while edge is not None:
                    finished, new_edge = self.edgeResolver.resolveEdge(edge, lines_list)
                    if finished:
                        # self.lineStorage.resolved_edges.add(edge.id)
                        # self.lineStorage.resolved_edges.add(edge.rc.id)
                        print "Successfully resolved edge", edge.id, "into", len(lines_list), "lines"
                        self.lineStorage.resolved_edges.add(edge.id)
                        # for line in lines_list:
                        #     self.lineStorage.edgeLines[line.chain[-1].seg_to.contig.id].append(line)
                        self.printer.printCurrentGraph([], [edge], "Resolved edge " + str(edge.id))
                    else:
                        print "Failed to resolve edge", edge.id
                        if new_edge is not None:
                            print "Graph modification detected. Trying to resolve again with edge:", new_edge.id, "of length", len(new_edge)
                            self.printer.printCurrentGraph([], [edge])
                            self.graph.printToFile(sys.stdout)
                    edge = new_edge

    def resolve(self):
        self.graph.printToFile(sys.stdout)
        visited_vertices = set()
        print "Resolving unique edges"
        # for edge in self.graph.E.values():
        #     if edge.info.unique:# and (len(edge.end.out) != 1 or len(edge.end.inc) != 1):
        #         # res = self.edgeResolver.processUniqueEdge(edge, self.lineStorage.getEdgeLines(edge)[0])
        #         # assert res is not None or len(edge.end.out) == 0
        #         # if res is not None:
        #         self.lineStorage.edgeLines[res.id].append(self.lineStorage.getEdgeLines(edge)[0])
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
        for v in self.graph.V.values():# Aligning reads to remaining lines
            if v.id in visited_vertices:
                continue
            for edge in v.inc:
                if not edge.info.unique:
                    continue
                line = self.lineStorage.edgeLines[edge.id][0]
                relevant_reads = [self.lineStorage.reads[read.id] for read in edge.reads if read.id not in line.reads.reads]
                line.addReads(relevant_reads)
                line.invalidated_reads.extend(relevant_reads)
                line.fixLineAlignments()
        Knotter(self.lineStorage, self.edgeResolver.aligner).knotGraph()
        self.printer.printCurrentGraph([], [])

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



