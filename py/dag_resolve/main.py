import os
import sys
import time

import shutil
from typing import Dict

sys.path.append("py")

from alignment.polishing import Polisher
from alignment.align_tools import DirDistributor, Aligner
from dag_resolve import params
from dag_resolve.edge_resolver import EdgeResolver
from dag_resolve.graph_resolver import GraphResolver
from dag_resolve.line_tools import LineStorage
from dag_resolve.repeat_graph import Graph, DotParser
from dag_resolve.sequences import ContigCollection, ReadCollection, UniqueList
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
    comp_dict = dict() # type: Dict[int, list[int]]
    for eid in components:
        if components[eid] not in comp_dict:
            comp_dict[components[eid]] = []
        comp_dict[components[eid]].append(eid)
    unique = len(filter(lambda edge: edge.info.unique, graph.E.values()))
    print "Stat: Unique edges:", unique
    print "Stat: Repeat edges:", len(graph.E) - unique
    print "Stat: Repeat components:", len(comp_dict)
    print "Stat: 2in2out:", len(filter(lambda comp: len(comp) == 1 and len(graph.E[comp[0]].start.inc) == 2, comp_dict.values()))
    print "Stat: Resolved edges:", len(storage.resolved_edges) - len(filter(lambda edge: edge.info.unique, graph.E.values()))
    print "Stat: Unique edge connections:", len(filter(lambda line: line.knot is not None, storage.lines))
    print "Stat: Loops:", len(filter(lambda line: line.isSimpleLoop(), storage.lines))


def ParseEdges(e_str):
    if e_str == "":
        return None
    edges = dict()  # type: Dict[int, str]
    for s in e_str.split(";"):
        s = s.split(",")
        edges[int(s[0])] = s[1:]
    for eid in list(edges.keys()):
        if -eid not in edges:
            s = []
            for tmp in edges[eid]:
                if tmp == "source":
                    s.append("sink")
                if tmp == "sink":
                    s.append("source")
            edges[-eid] = s
    return edges


if __name__ == "__main__":
    start = time.time()
    sys.stdout.write("Started\n")
    indir = sys.argv[1]
    if indir.endswith("/"):
        indir = indir[:-1]
    edges = ParseEdges(sys.argv[2])
    dir = sys.argv[3]
    edge_sequences = os.path.join(indir, "2-repeat", "graph_final.fasta")
    dot_file = os.path.join(indir, "2-repeat", "graph_final.dot")
    if not os.path.isfile(dot_file):
        dot_file = os.path.join(indir, "2-repeat", "graph_final.gv")
    reads_file = os.path.join(indir, os.path.split(indir)[1] + ".fasta")
    relevant_reads = os.path.join(indir, "2-repeat", "repeats_dump.txt")
    basic.ensure_dir_existance(dir)
    if params.clean:
        basic.recreate(dir)
    dir_distributor = DirDistributor(os.path.join(dir, "alignment"))
    old_logs_dir = os.path.join(dir, "old")
    basic.ensure_dir_existance(old_logs_dir)
    log_file = os.path.join(dir, "log.info")
    if os.path.isfile(log_file):
        num = len(os.listdir(old_logs_dir))
        shutil.copy(log_file, os.path.join(old_logs_dir, str(num) + ".log"))
    log = open(log_file, "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stderr = sys.stdout
    print " ".join(sys.argv)
    print (time.strftime("%d.%m.%Y  %I:%M:%S"))
    sys.stdout.write("Collecting contig collection\n")
    edge_sequences = ContigCollection().loadFromFasta(open(edge_sequences, "r"))
    sys.stdout.write("Loading dot\n")
    dot = DotParser(open(dot_file, "r")).parse(edges)
    graph = Graph().loadFromDot(edge_sequences, dot)
    graph.printToFile(sys.stdout)
    sys.stdout.write("Aligning reads\n")
    al = Aligner(dir_distributor)
    polisher = Polisher(al, dir_distributor)
    reads = ReadCollection().loadFromFasta(open(reads_file, "r"))
    alignment = al.align(reads, ContigCollection(list(UniqueList(graph.E.values()))))
    sys.stdout.write("Filling alignments\n")
    graph.fillAlignments(reads.asSeqRecords(), alignment, False)
    relevant = graph.fillRelevant(relevant_reads, reads)
    sys.stdout.write("Added " + str(relevant) + " reads\n")
    sys.stdout.write("Resolving repeats\n")
    picture_dir = os.path.join(dir, "pictures")
    basic.recreate(picture_dir)
    lineStorage = LineStorage(graph, al)
    res = GraphResolver(graph, picture_dir, lineStorage, EdgeResolver(graph, al, polisher, lineStorage.reads))
    res.resolve()
    sys.stdout.write("Printing results\n")
    res.printResults(sys.stdout)
    analyse(graph, res.lineStorage)
    print "Finished in " + str(time.time() - start) + " seconds"
    log.close()