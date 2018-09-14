import os
import sys


sys.path.append("py")

import alignment.polishing
import dag_resolve.edge_resolver
from common import basic
from dag_resolve import repeat_graph, sequences, graph_resolver, params
from alignment import align_tools

if __name__ == "__main__":
    sys.stdout.write("Started\n")
    edge_sequences = sys.argv[1]
    dot_file = sys.argv[2]
    reads = sys.argv[3]
    edges = dict()
    if sys.argv[4] == "":
        edges = None
    else:
        for s in sys.argv[4].split(";"):
            s = s.split(",")
            edges[int(s[0])] = s[1:]
    dir = sys.argv[5]
    basic.ensure_dir_existance(dir)
    if params.clean:
        basic.recreate(dir)
    dir_distributor = align_tools.DirDistributor(os.path.join(dir, "alignment"))
    log = open(os.path.join(dir, "log.info"), "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stdout.write("Collecting contig collection\n")
    edge_sequences = sequences.ContigCollection().loadFromFasta(open(edge_sequences, "r"))
    sys.stdout.write("Loading dot\n")
    dot = repeat_graph.DotParser(open(dot_file, "r")).parse(edges)
    graph = repeat_graph.Graph().loadFromDot(edge_sequences, dot)
    graph.printToFile(sys.stdout)
    sys.stdout.write("Aligning reads\n")
    al = align_tools.Aligner(dir_distributor)
    polisher = alignment.polishing.Polisher(al, dir_distributor)
    reads = sequences.ReadCollection().loadFromFasta(open(reads, "r"))
    alignment = al.align(reads, sequences.ContigCollection(graph.E.values()))
    sys.stdout.write("Filling alignments\n")
    graph.fillAlignments(reads.asSeqRecords(), alignment, False)
    sys.stdout.write("Resolving repeats\n")
    res = graph_resolver.GraphResolver(graph, graph_resolver.VertexResolver(graph, al, polisher), dag_resolve.edge_resolver.EdgeResolver(graph, al, polisher))
    res.resolve()
    sys.stdout.write("Printing results\n")
    res.printResults(sys.stdout)
    log.close()