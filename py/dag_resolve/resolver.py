import sys

from common import basic
from dag_resolve import repeat_graph, line, sequences, align_tools
import itertools

class LineTail:
    def __init__(self, line, edge, tail_consensus):
        # type: (line.Line, repeat_graph.Edge, align_tools.Consensus) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus

class VertexResolver:
    def __init__(self, graph, alignerObject):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> VertexResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = alignerObject

    def resolveVertex(self, v):
        # type: (repeat_graph.Vertex) -> Generator[LineTail]
        if len(v.out) > 1:
            assert False, "Vertex resolution is not implemented"
        out_edge = v.out[0]
        for edge in v.inc:
            for lineSegment in self.lineStorage.resolved_edges[edge.id]:
                reads = lineSegment.reads.inter(edge.suffix(500)).inter(out_edge.prefix(500))# type: sequences.ReadCollection
                tail = self.aligner.polishAndAnalyse(reads, out_edge.prefix(3000).subcontig())
                tail.printQuality(sys.stderr)
                yield LineTail(lineSegment.line, out_edge, tail)



class Phasing:
    def __init__(self, states):
        # type: (list[line.DivergenceState]) -> Phasing
        self.states = states

class EdgeResolver:
    def __init__(self, graph, aligner):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> EdgeResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = aligner

    def findDivergence(self, tails):
        # type: (list[str]) -> list[line.Divergence]
        assert False, "Divergence search not implemented"

    def determinePhasing(self, seq, consensus, div_list):
        # type: (str, sequences.Contig, list[line.Divergence]) -> Generator[line.DivergenceState]
        alignment = self.aligner.matchingAlignment([seq], consensus)[0]
        for div in div_list:
            pos = alignment.findSeqPos(div.pos)
            nighborhood = seq[pos - 10: pos + 10 + 1]
            weights = [basic.diff(nighborhood, state.neighborhood) for state in div.states]
            sorted_weights = sorted(weights)
            if sorted_weights[0] > 1.2 * sorted_weights[1]:
                yield div.ambiguous
            else:
                for i, w in enumerate(weights):
                    if w == sorted_weights[0]:
                        yield div.states[i]
                        break




    def resolveEdge(self, e, tails):
        # type: (repeat_graph.Edge, Generator[LineTail]) -> None
        tails = list(tails)

        assert False, "Edge resolution is not implemented"
        pass



class GraphResolver:
    def __init__(self, graph, vertexResolver, edgeResolver):
        # type: (repeat_graph.Graph, VertexResolver, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.lineStorage = line.LineStorage(graph)
        self.vertexResolver = vertexResolver
        vertexResolver.lineStorage = self.lineStorage
        self.edgeResolver = edgeResolver
        edgeResolver.lineStorage = self.lineStorage

    def resolveVertexForward(self, v):
        tails = sorted(self.vertexResolver.resolveVertex(v), key = lambda tail: tail.edge.id)
        for edge, tails_generator in itertools.groupby(tails, lambda tail: tail.edge):
            print edge.id
            for tail in tails_generator:
                print tail.tail_consensus
            self.edgeResolver.resolveEdge(edge, list(tails_generator))

    def resolve(self):
        print self.graph.V
        potentially_resolvable = self.graph.V.values()
        while True:
            cnt = 0
            for v in potentially_resolvable:
                print v.id, "oppa"
                if self.lineStorage.isResolvableLeft(v):
                    self.resolveVertexForward(v)
                    potentially_resolvable.remove(v)
                    cnt += 1
            if cnt == 0:
                break


