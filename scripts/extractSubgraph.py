import os
import sys


sys.path.append("py")
from common.SimpleGraph import SimpleGraph
from common import basic
from common.parsing import parseUPaths


def extractSubgraph(dir, flye_dir, contigs):
    basic.ensure_dir_existance(dir)
    d = parseUPaths(flye_dir)
    edge_ids = []
    for contig in contigs:
        for s in d[contig]:
            edge_ids.append(s)
    graph = SimpleGraph().ReadDot(os.path.join(flye_dir, "20-repeat", "graph_after_rr.gv"))
    vertex_ids = set()
    len = 0
    for eid in edge_ids:
        len += graph.e[eid].len
        vertex_ids.add(graph.e[eid].start)
        vertex_ids.add(graph.e[eid].end)
        if len > 10000:
            break
    # print "{|}|" + "|".join(["id " + r + "\\\\" for r in edge_ids])
    print "{|}|" + "|".join(["\"" + str(r) + "\"" for r in vertex_ids])


if __name__ == "__main__":
    dir = sys.argv[1]
    flye_dir = sys.argv[2]
    contigs = sys.argv[3:]
    extractSubgraph(dir, flye_dir, contigs)