import os
import sys


sys.path.append("py")

from alignment.polishing import Polisher
from alignment.align_tools import DirDistributor, Aligner
from dag_resolve import params
from dag_resolve.edge_resolver import EdgeResolver
from dag_resolve.graph_resolver import GraphResolver, VertexResolver
from dag_resolve.line_tools import LineStorage
from dag_resolve.repeat_graph import Graph, DotParser
from dag_resolve.sequences import ContigCollection, ReadCollection

from common import basic

def analyse(graph, storage):
    # type: (Graph, LineStorage) -> None
    components = dict()
    for edge in graph.E.values():
        if not edge.info.unique:
            components[edge.id] = edge.id
    for v in graph.V.values():
        vc = set()
        for e in v.inc:
            if e.id in components:
                vc.add(components[e.id])
        for e in v.out:
            if e.id in components:
                vc.add(components[e.id])
        if len(vc) > 1:
            nid = list(vc)[0]
            for eid in components:
                if components[eid] in vc:
                    components[eid] = nid
    comp_dict = dict()
    for eid in components:
        if components[eid] not in comp_dict:
            comp_dict[components[eid]] = []
        comp_dict[components[eid]].append(eid)
    unique = len(filter(lambda edge: edge.info.unique, graph.E.values()))
    print "Stat: Unique edges:", unique
    print "Stat: Repeat edges:", len(graph.E) - unique
    print "Stat: Repeat components:", len(comp_dict)
    print "Stat: 2in2out:", len(filter(lambda comp: len(comp) == 1, comp_dict.values()))
    print "Stat: Resolved edges:", len(storage.resolved_edges) - len(filter(lambda edge: edge.info.unique, graph.E.values()))
    print "Stat: Unique edge connections:", len(filter(lambda line: line.knot is not None, storage.lines))
    print "Stat: Loops:", len(filter(lambda line: line.isSimpleLoop(), storage.lines))

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
    dir_distributor = DirDistributor(os.path.join(dir, "alignment"))
    log = open(os.path.join(dir, "log.info"), "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stderr = sys.stdout
    print " ".join(sys.argv)
    sys.stdout.write("Collecting contig collection\n")
    edge_sequences = ContigCollection().loadFromFasta(open(edge_sequences, "r"))
    sys.stdout.write("Loading dot\n")
    dot = DotParser(open(dot_file, "r")).parse(edges)
    graph = Graph().loadFromDot(edge_sequences, dot)
    graph.printToFile(sys.stdout)
    sys.stdout.write("Aligning reads\n")
    al = Aligner(dir_distributor)
    polisher = Polisher(al, dir_distributor)
    reads = ReadCollection().loadFromFasta(open(reads, "r"))
    alignment = al.align(reads, ContigCollection(graph.E.values()))
    sys.stdout.write("Filling alignments\n")
    graph.fillAlignments(reads.asSeqRecords(), alignment, False)
    sys.stdout.write("Resolving repeats\n")
    res = GraphResolver(graph, VertexResolver(graph, al, polisher), EdgeResolver(graph, al, polisher))
    res.resolve()
    sys.stdout.write("Printing results\n")
    res.printResults(sys.stdout)
    analyse(graph, res.lineStorage)
    log.close()