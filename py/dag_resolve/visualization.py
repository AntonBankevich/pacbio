import os

from typing import Callable, Optional, Union, BinaryIO, List, Dict, Tuple

from common import basic
from common.basic import quoted
from common.sequences import UniqueList, Segment, MatchingSequence
from dag_resolve.line_align import Scorer
from dag_resolve.line_tools import LineStorage, Line
from dag_resolve.repeat_graph import Graph, Edge, Vertex
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.colors

class Palette:
    def __init__(self, num, cmap = "jet"):
        self.num = num
        self.cmap = cmap
        self.colorList = ["blue", "red", "green", "purple", "yellow", "violet", "orange", "pink", "brown", "navy", "honey"]
        if self.cmap == "basic" and num > len(self.colorList):
            self.cmap = "jet"
        if self.cmap != "basic":
            self.colorList = []
            for val in range(num):
                norm_value = float(val) / (self.num - 1)
                if norm_value < 0 or norm_value > 1:
                    print "ERROR: value " + str(norm_value) + ' does not belong to [0, 1]'
                cmap = plt.cm.get_cmap(self.cmap)
                color = cmap(norm_value)
                self.colorList.append(mplt.colors.rgb2hex(color[:3]))

    def getColor(self, val):
        return self.colorList[val]


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
        self.step = 40

    def splitEdges(self):
        good = self.goodEdges()
        print "Good edges:", map(str, list(good))
        res = dict() # type: Dict[str, List[Tuple[str, int]]]
        for edge in UniqueList(self.graph.E.values()):
            if edge.id not in good:
                res[edge.id] = [("V:" + str(edge.start.id), 0), ("V:" + str(edge.end.id), len(edge))]
                res[edge.rc.id] = [("V:" + str(edge.rc.start.id), 0), ("V:" + str(edge.rc.end.id), len(edge))]
            elif not edge.unique() or len(edge) < self.step * 30:
                n = max(1, len(edge) / self.step)
                step = len(edge) / n
                lens = [step] * n
                for i in range(len(edge) % n):
                    lens[i] += 1
                assert sum(lens) == len(edge)
                v = [("V:" + str(edge.start.id), 0)]
                v1 = [("V:" + str(edge.start.rc.id), len(edge))]
                cpos = 0
                for i in range(n - 1):
                    cpos += lens[i]
                    v.append((edge.id + "-" + str(cpos), cpos))
                    v1.append((edge.rc.id + "-" + str(len(edge) - cpos), len(edge) - cpos))
                v.append(("V:" + str(edge.end.id), len(edge)))
                v1.append(("V:" + str(edge.end.rc.id), 0))
                res[edge.id] = v
                res[edge.rc.id] = v1[::-1]
            else:
                v = [("V:" + str(edge.start.id), 0)]
                v1 = [("V:" + str(edge.start.rc.id), len(edge))]
                for i in range(1, 10):
                    cpos = i * self.step
                    v.append((edge.id + "-" + str(cpos), cpos))
                    v1.append((edge.rc.id + "-" + str(len(edge) - cpos), len(edge) - cpos))
                for i in range(1, 10)[::-1]:
                    cpos = len(edge) - i * self.step
                    v.append((edge.id + "-" + str(cpos), cpos))
                    v1.append((edge.rc.id + "-" + str(len(edge) - cpos), len(edge) - cpos))
                v.append(("V:" + str(edge.end.id), len(edge)))
                v1.append(("V:" + str(edge.end.rc.id), 0))
                res[edge.id] = v
                res[edge.rc.id] = v1[::-1]
        return res

    def goodEdges(self):
        visited = set()
        queue = [] #type: List[Vertex, bool]
        order = []
        for v in self.graph.V.values():
            if v.id in visited:
                continue
            visited.add(v.id)
            queue.append((v, True))
            for e in v.out:
                if not e.unique():
                    queue.append((e.end, False))
            while len(queue) != 0:
                v, finished = queue.pop()
                if finished:
                    order.append(v)
                else:
                    if v.id in visited:
                        continue
                    visited.add(v.id)
                    queue.append((v, True))
                    for e in v.out:
                        if not e.unique:
                            queue.append((v, False))
        visited = set()
        queue = [] #type: List[Vertex, bool]
        group = dict()
        for v in order:
            gid = v.id
            if v.id in visited:
                continue
            visited.add(v.id)
            queue.append((v, True))
            group[v.id] = gid
            for e in v.out:
                if not e.unique():
                    queue.append((e.end, False))
            while len(queue) != 0:
                v, finished = queue.pop()
                if finished:
                    order.append(v)
                else:
                    if v.id in visited:
                        continue
                    visited.add(v.id)
                    group[v.id] = gid
                    queue.append((v, True))
                    for e in v.out:
                        if not e.unique:
                            queue.append((v, False))
        res = set()
        for e in self.graph.E.values():
            if e.unique() or group[e.start.id] != group[e.end.id]:
                res.add(e.id)
        return res

    def printToFile(self, handler, additionalEdgeColorings = None, additionalVertexColorings = None):
        vIndex = self.splitEdges()
        handler.write("digraph {\n")
        handler.write("nodesep = 0.5;\n")
        handler.write("node[shape = circle, label = \"\", height = 0.3];\n")
        for v in self.graph.V.values():
            handler.write(basic.quoted("V:" + str(v.id)) + " [label=\"V:"+ str(v.id) + "\", height = 0.6, fillcolor = " + basic.quoted("black") + "];\n")
        for edge in self.graph.E.values():
            v_list = vIndex[edge.id]
            if len(vIndex[edge.id]) >=4:
                handler.write(basic.quoted(v_list[1][0]) + " [label = " + basic.quoted("E:" + str(edge.id)) + "];\n")
                for v in v_list[2:-2]:
                    handler.write(basic.quoted(v[0]) + " [style = \"filled\", fillcolor = \"" + "black" + "\"];\n")
                handler.write(basic.quoted(v_list[-2][0]) + " [label = " + basic.quoted("E:" + str(edge.id)) + "];\n")
            elif len(vIndex[edge.id]) == 3:
                handler.write(basic.quoted(v_list[1][0]) + " [label = " + basic.quoted("E:" + str(edge.id)) + "];\n")
        color_map = dict() # type: Dict[str, str]
        color_map["graph"] = "black"
        colors = Palette(len(self.storage.lines))
        for i, line in enumerate(self.storage.lines):
            color_map[line.id] = colors.getColor(i)

        edges = self.collectLineEdges(color_map, handler, vIndex)
        for edge in self.graph.E.values():
            v_list = vIndex[edge.id]
            for v1, v2 in zip(v_list[:-1], v_list[1:]):
                if v2[1] - v1[1] < 2 * self.step:
                    key = (v1[0], v2[0])
                    if key not in edges:
                        edges[key] = []
                    edges[key].append((Segment(edge, v1[1], v2[1]), "graph"))
        for v1, v2 in edges.keys():
            sequences = edges[(v1, v2)]
            scores = self.calculateScores(map(lambda x: x[0], sequences))
            group, group_lists = self.splitIntoGroups(scores, sequences)
            min_dists = self.countMinDists(group, group_lists, scores, sequences)
            tooltips = [sequence[0].Seq() for sequence in sequences]
            tooltips.append(";".join(str(sequence[1]) + ":" + color_map[sequence[1]] for sequence in sequences))
            for gList, min_dist in zip(group_lists, min_dists):
                if len(gList) == 0:
                    continue
                handler.write(quoted(v1) + " -> " + quoted(v2))
                handler.write("[color = " + quoted(":".join([color_map[sequences[i][1]] for i in gList])))
                style = "solid"
                if min_dist < 15:
                    style = "dotted"
                elif min_dist < 30:
                    style = "dashed"
                handler.write(" , style = " + basic.quoted(style))
                handler.write(" , tooltip = " + basic.quoted("\n".join(tooltips)))
                handler.write(" , penwidth = 3] ;\n")
        handler.write("}\n")

    def splitIntoGroups(self, scores, sequences):
        group = []
        n = len(sequences)
        for i in range(n):
            group.append(i)
        for i in range(n):
            for j in range(n):
                if scores[i][j] < 10:
                    group[i] = group[j]
        group_lists = []
        for i in range(n):
            group_lists.append([])
        for i in range(n):
            cur = i
            while group[cur] != cur:
                cur = group[cur]
            group_lists[cur].append(i)
            group[i] = cur
        return group, group_lists

    def countMinDists(self, group, group_lists, scores, sequences):
        n = len(sequences)
        min_dists = [100] * n
        black_group = -1
        for i in range(n):
            if len(group_lists[i]) == 1 and sequences[group_lists[i][0]][1] == "graph":
                black_group = i
        for i in range(n):
            if i == black_group:
                continue
            for j in range(n):
                if j == black_group:
                    continue
                if group[i] != group[j]:
                    min_dists[group[i]] = min(min_dists[group[i]], scores[i][j])
                    min_dists[group[j]] = min(min_dists[group[j]], scores[i][j])
        return min_dists

    def collectLineEdges(self, color_map, handler, vIndex):
        edges = dict()  # type: Dict[Tuple[str, str], List[Tuple[Segment, str]]]
        for line in self.storage.lines:
            line_sequence = [] # type: List[Tuple[str, int, Union[Line, Edge]]]
            edge_alignments = []
            cur = [line.chain[0]]
            for al in line.chain[1:]:
                if cur[-1].precedes(al):
                    cur.append(al)
                else:
                    edge_alignments.append(cur)
                    cur = [al]
            edge_alignments.append(cur)
            for als in edge_alignments:
                edge = als[0].seg_to.contig
                matchings = als[0].matchingSequence().concat(map(lambda al: al.matchingSequence(), als[1:]))
                cur_matching = 0
                for v, pos in vIndex[edge.id]:
                    if len(line_sequence) > 0 and line_sequence[-1][0] == v:
                        continue
                    while cur_matching < len(matchings.matches) and matchings.matches[cur_matching][1] < pos:
                        cur_matching += 1
                    if cur_matching < len(matchings.matches) and pos <= matchings.matches[cur_matching][1] < pos + self.step / 5:
                        line_sequence.append((v, matchings.matches[cur_matching][0], edge))
                        assert len(line_sequence) < 2 or line_sequence[-2][1] <= line_sequence[-1][1], ",".join(map(str, als)) + "\n" + str(line_sequence) + "\n" + v
                    if cur_matching == len(matchings.matches) and len(matchings.matches) >= 1 and matchings.matches[-1][0] >= pos - self.step / 5:
                        line_sequence.append((v, matchings.matches[-1][0] + 1, edge))
                        assert len(line_sequence) < 2 or line_sequence[-2][1] <= line_sequence[-1][1]
                        break
            if line_sequence[0][1] > self.step / 2:
                line_sequence.insert(0, ("L" + str(line.id) + "-" + str(0), 0, line))
                handler.write(basic.quoted(line_sequence[0][0]) + " [style = \"filled\", fillcolor = \"" + color_map[
                    line.id] + "\"];\n")
            if line_sequence[-1][1] < len(line) - self.step / 2:
                line_sequence.append(("L" + str(line.id) + "-" + str(len(line_sequence)), len(line), line))
                handler.write(basic.quoted(line_sequence[-1][0]) + " [style = \"filled\", fillcolor = \"" +
                              color_map[line.id] + "\"];\n")
            line_vertices = [line_sequence[0]]
            for v1, v2 in zip(line_sequence[:-1], line_sequence[1:]):
                if v1[2] == v2[2] and len(line_sequence) != 2 and v1[2].unique() and v2[1] - v1[1] > 2 * self.step:
                    continue
                if v2[1] > v1[1] + self.step * 2:
                    n = max(1, (v2[1] - v1[1]) / self.step)
                    step = (v2[1] - v1[1]) / n
                    lens = [step] * n
                    for i in range((v2[1] - v1[1]) % n):
                        lens[i] += 1
                    assert sum(lens) == (v2[1] - v1[1])
                    cpos = v1[1]
                    for i in range(n - 1):
                        cpos += lens[i]
                        v_name = "L" + str(line.id) + "-" + str(cpos)
                        line_vertices.append((v_name, cpos))
                        assert len(line_vertices) < 2 or line_vertices[-2][1] <= line_vertices[-1][1]
                        handler.write(basic.quoted(v_name) + " [style = \"filled\", fillcolor = \"" +
                                      color_map[line.id] + "\"];\n")
                line_vertices.append(v2)
                assert len(line_vertices) < 2 or line_vertices[-2][1] <= line_vertices[-1][1], str(v1) + " " +  str(v2)
            for v1, v2 in zip(line_vertices[:-1], line_vertices[1:]):
                if v2[1] - v1[1] > 2 * self.step:
                    continue
                key = (v1[0], v2[0])
                if key not in edges:
                    edges[key] = []
                edges[key].append((Segment(line, v1[1], v2[1]), line.id))
        return edges

    def calculateScores(self, segments):
        # type: (List[Segment]) -> List[List[int]]
        res = []
        scorer = Scorer()
        for seg1 in segments:
            res_line = []
            for seg2 in segments:
                if seg1.expand(10).Seq().find(seg2.Seq()) != -1 or seg2.expand(10).Seq().find(seg1.Seq()) != -1:
                    res_line.append(0)
                    continue
                tmp = MatchingSequence(seg1.contig.seq, seg2.contig.seq, [(seg1.left, seg2.left), (seg1.right - 1, seg2.right - 1)])
                res_line.append(scorer.accurateScore(tmp))


            res.append(res_line)
        return res

