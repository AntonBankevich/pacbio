import itertools
import sys

import os
from typing import Optional

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import sam_parser, basic
from dag_resolve.edge_resolver import EdgeResolver
from dag_resolve.filters import EdgeTransitionFilter
from dag_resolve.knots import Knotter
from dag_resolve.line_tools import LineTail, LineSegment, LineStorage
from dag_resolve.repeat_graph import Graph, Edge, Vertex
from dag_resolve.sequences import Contig, ReadCollection, Segment, ContigCollection
from dag_resolve.visualization import DotPrinter, FilterColoring


class VertexResolver:
    def __init__(self, graph, aligner, polisher):
        # type: (Graph, Aligner, Polisher) -> VertexResolver
        self.graph = graph
        self.lineStorage = None
        self.polisher = polisher
        self.aligner = aligner

    def transitionSupport(self, prev, seg, next):
        # type: (Edge, LineSegment, Edge) -> ReadCollection
        return seg.reads.filter(EdgeTransitionFilter(prev, next))

    def resolveVertex(self, v):
        # type: (Vertex) -> bool
        print "Resolving vertex", v.id
        print "Incoming edges:", map(lambda e: e.id, v.inc)
        print "Outgoing edges:", map(lambda e: e.id, v.out)
        assert len(v.out) > 0, "Trying to resolve hanging vertex " + str(v.id)
        # for edge in v.inc:
        #     for lineSegment in self.lineStorage.edgeLines(edge):# type:LineSegment
        #         print ""
        #         print "Reads from line", lineSegment.line.id
        #         lineSegment.reads.print_alignments(sys.stdout)
        for edge in v.out:
            if len(edge) < 500:
                print "One of outgoing edges is too short. Skipping vertex"
                return False
        for edge in v.inc:
            for lineSegment in self.lineStorage.edgeLines(edge):# type:LineSegment
                print "Creating a tail for line", lineSegment.line
                base_seq = lineSegment.seq[-min(5000, len(lineSegment.seq)):]
                consensus = self.polisher.polishQuiver(lineSegment.reads.inter(lineSegment.edge.suffix(-5000)), base_seq, len(base_seq), 3000)
                consensus = consensus.cut()
                print "Created consensus of length", len(consensus)
                first = None
                collection = ReadCollection(ContigCollection(v.out))
                collection.loadFromSam(self.aligner.align([Contig(consensus.seq, 0)], ContigCollection(v.out)))
                print collection.reads["0"].__str__()
                for rec in collection.reads["0"].alignments:#v.out):
                    if rec.seg_to.contig in v.out and len(rec.seg_from) > 300 and \
                            (len(rec.seg_from) > 500 or rec.seg_from.right > len(consensus.seq) - 100) and \
                            (first is None or rec.seg_from.left < first.seg_from.left):
                        first = rec
                if first is None:
                    print "Could not connect one of the tails. Aborting."
                    return False
                print "Found connection of consensus to position", first.seg_to.left, "of edge", first.seg_to.contig.id
                read_pos = first.seg_from.left
                new_edge = first.seg_to.contig
                filtered_reads = lineSegment.reads.inter(lineSegment.edge.suffix(-5000)).inter(Segment(new_edge, first.seg_to.left, first.seg_to.right))
                if first.seg_to.left < 500:
                    print first.seg_from.left, "nucleotedes were hidden inside a vertex"
                    lineSegment.line.tail = LineTail(lineSegment.line, new_edge, consensus.suffix(read_pos), filtered_reads)
                    lineSegment.line.tail.alignment = self.aligner.matchingAlignment([lineSegment.line.tail.tail_consensus.seq], new_edge)[0]
                else:
                    if read_pos < 500:
                        print "Late entrances not supported. Aborting."
                        return False
                    print "Reparing graph and restarting vertex resolution"
                    self.graph.addCuttingEdge(new_edge, 0, new_edge, first.seg_to.left, consensus.seq[:read_pos])
                    self.aligner.repairGraphAlignments(self.graph)
                    self.graph.printToFile(sys.stdout)
                    return self.resolveVertex(v)
                    # lineSegment.line.tail = LineTail(lineSegment.line, new_edge, consensus.cut(length=read_pos), filtered_reads)


                # support = []
                # for next_candidate in v.out:
                #     support.append(self.transitionSupport(edge, lineSegment, next_candidate))
                # if len(v.out) == 1:
                #     next = 0
                # else:
                #     supportWeight = map(len, support)
                #     print "Number of reads that support transitions for the line from edge", edge.id, ":", ",".join(
                #         map(str, supportWeight))
                #     next, alternative = basic.best2(supportWeight, lambda x, y: x > y)
                #     if supportWeight[alternative] * 5 > supportWeight[next]:
                #         print "Support of alternative edge is too high. Aborting."
                #         return False
                # reads = support[next] # type: sequences.ReadCollection
                # tail = self.polisher.polishAndAnalyse(reads, v.out[next])
                # lineSegment.line.tail = LineTail(lineSegment.line, v.out[next], tail, reads)
                # print "Created tail of line", lineSegment.line, "on edge", v.out[next].id, "of length", len(tail.cut())
        return True


class GraphResolver:
    def __init__(self, graph, dir, vertexResolver, edgeResolver):
        # type: (Graph, str, VertexResolver, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.dir = dir
        self.lineStorage = LineStorage(graph)
        self.vertexResolver = vertexResolver
        vertexResolver.lineStorage = self.lineStorage
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
        for edge in v.out:
            lines_list = self.lineStorage.getEdgeLines(edge)
            if len(lines_list) == 0:
                print "No line entered edge " + str(edge.id) + ". Skipping vertex."
        for edge in v.out:
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
        for line in self.lineStorage.lines:
            if line.id in printed or line.knot is not None:
                continue
            while line is not None:
                printed.add(line.id)
                printed.add(line.rc.id)
                handler.write(line.rcStr())
                handler.write("-")
                handler.write(line.rc.__str__())
                if line.rc.knot is None:
                    line = None
                else:
                    handler.write("->")
                    line = line.rc.knot.line2
            handler.write("\n")

        for line in self.lineStorage.lines:
            if line.id in printed:
                continue
            handler.write("->")
            while line is not None and line.id not in printed:
                printed.add(line.id)
                printed.add(line.rc.id)
                handler.write(line.rcStr())
                handler.write("-")
                handler.write(line.rc.__str__())
                if line.rc.knot is None:
                    line = None
                else:
                    handler.write("->")
                    line = line.rc.knot.line2
            handler.write("->\n")



