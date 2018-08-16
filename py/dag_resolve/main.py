import sys
import os

sys.path.append("py")
from common import basic
from dag_resolve import repeat_graph, sequences, align_tools, resolver

if __name__ == "__main__":
    sys.stderr.write("Started\n")
    edge_sequences = sys.argv[1]
    dot_file = sys.argv[2]
    reads = sys.argv[3]
    edges = map(int, sys.argv[4].split(";"))
    dir = sys.argv[5]
    basic.ensure_dir_existance(dir)
    sys.stderr.write("Collecting contig collection\n")
    edge_sequences = sequences.ContigCollection().loadFromFasta(open(edge_sequences, "r"))
    sys.stderr.write("Loading dot\n")
    dot = repeat_graph.DotParser(open(dot_file, "r")).parse(edges)
    graph = repeat_graph.Graph().loadFromDot(edge_sequences, dot)
    sys.stderr.write("Aligning reads\n")
    al = align_tools.Aligner(os.path.join(dir, "alignment"))
    reads = sequences.ReadCollection().loadFromFasta(open(reads, "r"))
    alignment = al.align(reads, sequences.ContigCollection(graph.E.values()))
    sys.stderr.write("Filling alignments\n")
    graph.fillAlignments(reads.asSeqRecords(), alignment)
    sys.stderr.write("Resolving repeats\n")
    res = resolver.GraphResolver(graph, resolver.VertexResolver(graph, al), resolver.EdgeResolver(graph))
    res.resolve()
    sys.stderr.write("Printing results\n")
    res.lineStorage.printToFile(sys.stdout)