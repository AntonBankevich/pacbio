from dag_resolve import repeat_graph, line
import itertools

class LineTail:
    def __init__(self, line, edge, tail_consensus):
        # type: (line.Line, repeat_graph.Edge, basestring) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus

class VertexResolver:
    def __init__(self, graph, lineStorage):
        # type: (repeat_graph.Graph, line.LineStorage) -> VertexResolver
        self.graph = graph
        self.lineStorage = lineStorage

    def resolveVertex(self, v):
        # type: (repeat_graph.Vertex) -> list[LineTail]
        assert False, "Vertex resolution is not implemented"
        pass

class EdgeResolver:
    def __init__(self, graph, lineStorage):
        # type: (repeat_graph.Graph, line.LineStorage) -> EdgeResolver
        self.graph = graph
        self.lineStorage = lineStorage

    def resolveEdge(self, e, tails):
        # type: (repeat_graph.Edge, list[LineTail]) -> None
        assert False, "Edge resolution is not implemented"
        pass



class GraphResolver:
    def __init__(self, graph, lineStorage, vertexResolver, edgeResolver):
        # type: (repeat_graph.Graph, line.LineStorage, VertexResolver, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.graph = repeat_graph.Graph()
        self.lineStorage = lineStorage
        self.vertexResolver = vertexResolver
        self.edgeResolver = edgeResolver

    def ResolveVertexForward(self, v):
        tails = sorted(self.vertexResolver.resolveVertex(v), lambda tail: tail.edge.id)
        for edge, tails_generator in itertools.groupby(tails, lambda tail: tail.edge):
            self.edgeResolver.resolveEdge(edge, list(tails_generator))

