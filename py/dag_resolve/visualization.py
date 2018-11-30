import os

from typing import Callable, Optional, Union, BinaryIO, List, Dict, Tuple

from common import basic
from common.basic import quoted
from common.sequences import UniqueList
from dag_resolve.line_tools import LineStorage
from dag_resolve.repeat_graph import Graph, Edge, Vertex
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.colors

def GetColorByNormalizedValue(cmap_name, norm_value):
    if norm_value < 0 or norm_value > 1:
        print "ERROR: value " + str(norm_value) + ' does not belong to [0, 1]'
    cmap = plt.cm.get_cmap(cmap_name)
    color = cmap(norm_value)
    return mplt.colors.rgb2hex(color[:3])

def FilterColoring(filter, color):
    # type: (Callable[[Union[Edge, Vertex]], bool], str) -> Callable[[Union[Edge, Vertex]], Optional[str]]
    return lambda edge: color if filter(edge) else None

def SimplePrintDot(graph, file):
    printer = DotPrinter(graph, [FilterColoring(lambda edge: not edge.info.unique, "red")], [])
    printer.printToFile(open(file, "w"))

class HistoryPrinter:
    def __init__(self, printer, dir):
        self.printer = printer
        basic.recreate(dir)
        self.dir = dir
        self.cur_picture = 1

    def printCurrentGraph(self, vertices = None, edges = None, message = ""):
        # type: (Optional[List[Vertex]], Optional[List[Edge]], str) -> None
        if vertices is None:
            vertices = []
        if edges is None:
            edges = []
        fn = os.path.join(self.dir, str(self.cur_picture) + "_" + "_".join(message.split()) + ".dot")
        f = open(fn, "w")
        self.printer.printToFile(f, [FilterColoring(lambda e: e in edges, "purple")], [FilterColoring(lambda v: v in vertices, "purple")])
        f.close()
        self.cur_picture += 1


class DotPrinter:
    def __init__(self, graph, edge_colorings = None, vertex_colorings = None):
        # type: (Graph, Optional[List[Callable[[Edge], Optional[str]]]], Optional[List[Callable[[Edge], Optional[str]]]]) -> DotPrinter
        if edge_colorings is None:
            edge_colorings = []
        if vertex_colorings is None:
            vertex_colorings = []
        self.edge_colorings = edge_colorings # type: List[Callable[[Edge], Optional[str]]]
        self.vertex_colorings = vertex_colorings # type: List[Callable[[Edge], Optional[str]]]
        self.graph = graph
        self.defaultEdgeColor = "black"
        self.defaultVertexColor = "white"

    def getEdgeColor(self, e, additionalColorings = None):
        # type: (Edge, Optional[List[Callable[[Edge], Optional[str]]]]) -> str
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
        # type: (Vertex, Optional[List[Callable[[Vertex], Optional[str]]]]) -> str
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
        # type: (BinaryIO, Optional[List[Callable[[Edge], Optional[str]]]], Optional[List[Callable[[Vertex], Optional[str]]]]) -> None
        handler.write("digraph {\n")
        handler.write("nodesep = 0.5;\n")
        handler.write("node[shape = circle, label = \"\", height = 0.3];\n")
        for v in self.graph.V.values():
            handler.write(quoted(v.id) + " [style = \"filled\", fillcolor = \"" + self.getVertexColor(v, additionalVertexColorings) + "\"];\n")
        for e in self.graph.E.values():
            handler.write(quoted(e.start.id) + " -> " + quoted(e.end.id))
            handler.write(" [label = \"id " + str(e.id) + "\\l" + str((len(e) // 100) * 0.1) + "k ")
            handler.write(str(e.info.cov) + "x\", color = " + quoted(self.getEdgeColor(e, additionalEdgeColorings)))
            handler.write(" , penwidth = 3] ;\n")
        handler.write("}\n")

class DotLinePrinter:
    def __init__(self, graph, storage):
        # type: (Graph, LineStorage) -> None
        self.graph = graph
        self.storage = storage
        self.step = 100

    def splitEdges(self):
        res = dict() # type: Dict[str, List[Tuple[str, int]]]
        for edge in UniqueList(self.graph.E.values()):
            if not edge.unique():
                n = min(1, len(edge) / self.step)
                step = len(edge) / n
                lens = [step] * n
                for i in range(len(edge) % n):
                    lens[i] += 1
                assert sum(lens) == len(edge)
                v = [(edge.start.id, 0)]
                v1 = [(edge.start.rc.id, len(edge))]
                cpos = 0
                for i in range(n - 1):
                    cpos += lens[i]
                    v.append((edge.id + "-" + str(cpos), cpos))
                    v1.append((edge.rc.id + "-" + str(len(edge) - cpos), len(edge) - cpos))
                v.append((edge.end.id, len(edge)))
                v1.append((edge.end.rc.id, 0))
                res[edge.id] = v
                res[edge.rc.id] = v1[::-1]
        return res

    def printToFile(self, handler, additionalEdgeColorings = None, additionalVertexColorings = None):
        vIndex = self.splitEdges()
        handler.write("digraph {\n")
        handler.write("nodesep = 0.5;\n")
        handler.write("node[shape = circle, label = \"\", height = 0.3];\n")
        for v in self.graph.V.values():
            handler.write(basic.quoted(v.id) + " [style = \"filled\", height = 0.6, fillcolor = " + basic.quoted("black") + "];\n")
        mpl.colo
        for edge in self.graph.E.values():
            for v in vIndex[edge.id][1:-1]:
                handler.write(basic.quoted(v[0]) + " [style = \"filled\", fillcolor = \"" + "black" + "\"];\n")

        edges = dict() # type: Dict[Tuple[str, str], Tuple[str, str]]



