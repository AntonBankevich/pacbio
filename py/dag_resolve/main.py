import sys
import os

sys.path.append("py")
from common import basic
from dag_resolve import repeat_graph, sequences, align_tools, resolver

if __name__ == "__main__":
    sys.stdout.write("Started\n")
    edge_sequences = sys.argv[1]
    dot_file = sys.argv[2]
    reads = sys.argv[3]
    edges = dict()
    for s in sys.argv[4].split(";"):
        s = s.split(",")
        edges[int(s[0])] = s[1:]
    dir = sys.argv[5]
    basic.ensure_dir_existance(dir)
    # sys.stderr = open(os.path.join(dir, "stderr.log"), "w")
    sys.stdout.write("Collecting contig collection\n")
    edge_sequences = sequences.ContigCollection().loadFromFasta(open(edge_sequences, "r"))
    sys.stdout.write("Loading dot\n")
    dot = repeat_graph.DotParser(open(dot_file, "r")).parse(edges)
    graph = repeat_graph.Graph().loadFromDot(edge_sequences, dot)
    sys.stdout.write("Aligning reads\n")
    al = align_tools.Aligner(os.path.join(dir, "alignment"))
    reads = sequences.ReadCollection().loadFromFasta(open(reads, "r"))
    alignment = al.align(reads, sequences.ContigCollection(graph.E.values()))
    sys.stdout.write("Filling alignments\n")
    graph.fillAlignments(reads.asSeqRecords(), alignment)
    sys.stdout.write("Resolving repeats\n")
    res = resolver.GraphResolver(graph, resolver.VertexResolver(graph, al), resolver.EdgeResolver(graph, al))
    res.resolve()
    sys.stdout.write("Printing results\n")
    res.lineStorage.printToFile(sys.stdout)