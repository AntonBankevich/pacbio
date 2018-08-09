from dag_resolve import graph, line


class GraphResolver:
    def __init__(self, g, lineStorage):
        self.g = g
        self.g = graph.Graph()
        self.lineStorage = lineStorage
        self.lineStorage = line.LineStorage(g)

    def ResolveVertex(self, v):
        # type: (graph.Vertex) -> None
        lines = []
