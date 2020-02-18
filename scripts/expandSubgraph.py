import os
import sys


sys.path.append("py")
from common.SimpleGraph import SimpleGraph
from common import basic
from common.parsing import parseUPaths


def extractSubgraph(gf):
#    graph = SimpleGraph().ReadDot(os.path.join(flye_dir, "20-repeat", "graph_after_rr.gv"))
    graph1 =SimpleGraph().ReadDot(os.path.join(gf))
    vertex_ids = graph1.v.keys()
    # print "{|}|" + "|".join(["id " + r + "\\\\" for r in edge_ids])
    print "{|}|" + "|".join(["\"" + str(r) + "\"" for r in vertex_ids])
    print " ".join(graph1.e.keys())


if __name__ == "__main__":
    gf = sys.argv[1]
#    flye_dir = sys.argv[2]
    extractSubgraph(gf)
