import sys

import os

from dag_resolve.repeat_graph import Graph
from dag_resolve.sequences import ContigCollection


def analyse(name, indir, outdir):
    # type: (str, str) -> None
    graph = Graph()
    contigs = os.path.join(indir, "2-repeat", "graph_final.fasta")
    dot = os.path.join(indir, "2-repeat", "graph_final.dot")
    edge_sequences = ContigCollection().loadFromFasta(open(contigs, "r"))
    graph.loadFromDot(ContigCollection())