from typing import Callable, Optional, Union

from common import basic
from dag_resolve.repeat_graph import Graph, Edge, Vertex


def FilterColoring(filter, color):
    # type: (Callable[[Union[Edge, Vertex]], bool], str) -> Callable[[Union[Edge, Vertex]], Optional[str]]
    return lambda edge: color if filter(edge) else None

class DotPrinter:
    def __init__(self, graph, edge_colorings = None, vertex_colorings = None):
        # type: (Graph, Optional[list[Callable[[Edge], Optional[str]]]]) -> DotPrinter
        if edge_colorings is None:
            edge_colorings = []
        if vertex_colorings is None:
            vertex_colorings = []
        self.edge_colorings = edge_colorings # type: list[Callable[[Edge], Optional[str]]]
        self.vertex_colorings = vertex_colorings # type: list[Callable[[Edge], Optional[str]]]
        self.graph = graph
        self.defaultEdgeColor = "black"
        self.defaultVertexColor = "white"

    def quoted(self, val):
        return "\"" + str(val) + "\""

    def getEdgeColor(self, e, additionalColorings = None):
        # type: (Edge, Optional[list[Callable[[Edge], Optional[str]]]]) -> str
        if additionalColorings is None:
            additionalColorings = []
        res = set()
        for coloring in self.edge_colorings + additionalColorings:
            color = coloring(e)
            if color is not None:
                res.add(color)
        if len(res) == 0:
            return self.defaultEdgeColor
        return ":".join(res)

    def getVertexColor(self, v, additionalColorings = None):
        # type: (Vertex, Optional[list[Callable[[Vertex], Optional[str]]]]) -> str
        if additionalColorings is None:
            additionalColorings = []
        res = set()
        for coloring in self.vertex_colorings + additionalColorings:
            color = coloring(v)
            if color is not None:
                res.add(color)
        if len(res) == 0:
            return self.defaultVertexColor
        return ":".join(res)

    def printToFile(self, handler, additionalEdgeColorings = None, additionalVertexColorings = None):
        # type: (file, Optional[list[Callable[[Edge], Optional[str]]]], Optional[list[Callable[[Vertex], Optional[str]]]]) -> None
        handler.write("digraph {\n")
        handler.write("nodesep = 0.5;\n")
        handler.write("node[shape = circle, label = \"\", height = 0.3];\n")
        for v in self.graph.V.values():
            handler.write(self.quoted(v.id) + " [style = \"filled\", fillcolor = \"" + self.getVertexColor(v, additionalVertexColorings) + "\"];\n")
        for e in self.graph.E.values():
            handler.write(self.quoted(e.start.id) + " -> " + self.quoted(e.end.id))
            handler.write(" [label = \"id " + str(e.id) + "\\l" + str((len(e) // 100) * 0.1) + "k ")
            handler.write(str(e.info.cov) + "x\", color = " + self.quoted(self.getEdgeColor(e, additionalEdgeColorings)))
            handler.write(" , penwidth = 3] ;\n")
        handler.write("}\n")
