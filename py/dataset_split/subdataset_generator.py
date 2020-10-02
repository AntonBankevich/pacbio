import os
import sys


sys.path.append("py")
from typing import List, Dict, Generator, Tuple


from common.basic import CreateLog
from common import basic
from common.SimpleGraph import SimpleGraph
from common.save_load import TokenReader, TokenWriter
from dataset_split.graph_splitting import SplitGraph, DipolidCalculator, HaploidCalculator


class ComponentRecord:
    def __init__(self, component, cov):
        # type: (SimpleGraph, float) -> None
        self.component = component
        self.cov = cov
        self.reads = []
        self.unique = dict()#type: Dict[str, List[str]]
        self.red = 0
        self.zero = 0
        self.bad_border = 0
        self.unresolved_connections = 0
        self.inc = 0
        self.out = 0
        self.resolved_connections = []
        self.outside_connections = 0

    def half(self):
        res = set()
        tmp = [iter(self.component.v).next()]
        while tmp.__len__() > 0:
            vid = tmp.pop()
            if vid in res:
                continue
            res.add(vid)
            v = self.component.v[vid]
            for e in v.inc + v.out:
                tmp.append(e.start)
                tmp.append(e.end)
        return res


    def addUniqueEdge(self, eid):
        eid = basic.Normalize(eid)
        if eid not in self.unique:
            self.unique[eid] = []

    def addRead(self, rid, all_edges):
        # type: (str, List[str]) -> None
        edges = [eid for eid in all_edges if eid in self.component.e]
        if edges.__len__() == 1 and (basic.Normalize(edges[0]) in self.unique):
            return
        self.reads.append(rid)
        eset = set()
        for edge in edges:
            edge = basic.Normalize(edge)
            if edge in self.unique:
                eset.add(edge)
        for edge in eset:
            self.unique[edge].append(rid)

    def calcStats(self):
        for v in self.component.v:
            if self.component.isHub(v):
                self.red += 1
        for e in self.component.e.values():
            if (basic.Normalize(e.id) not in self.unique) and (self.component.isBorder(e.start) or self.component.isBorder(e.end)):
                self.bad_border += 1
            if e.cov == 0:
                self.zero += 1

    def printStats(self, fname):
        f = open(fname, "w")
        f.write("Cov " + str(self.cov) + "\n size " + str(self.component.e.__len__()) + "\n")
        f.write("Zero " + str(self.zero) + "\n bad_border " + str(self.bad_border) + "\n hubs " + str(self.red) + "\n")
        f.close()

    def printReads(self, fname):
        f = open(fname, "w")
        f.write(str(self.reads.__len__()) + " " + str(self.reads.__len__()) + "\n")
        for rid in self.reads:
            f.write(rid + "\n")
        f.close()

    def printContigs(self, fname):
        f = open(fname, "w")
        f.write(str(self.unique.__len__()) + "\n")
        for eid in self.unique:
            f.write(">" + eid + "\n")
            edge = self.component.e[eid]
            seq = edge.seq
            if len(seq) > 50000:
                if self.component.isBorder(edge.start):
                    seq = seq[-20000:]
                elif self.component.isBorder(edge.end):
                    seq = seq[:20000]
            f.write(seq)
            f.write("\n")
        f.close()

    def printInit(self, fname):
        file = open(fname, "w")
        f = TokenWriter(file)
        f.writeIntLine(self.unique.__len__())
        for eid in self.unique:
            f.writeTokenLine(eid)
            for rid in self.unique[eid]:
                f.writeToken(rid)
            f.newLine()
        file.close()

    def dump(self, dirname):
        basic.ensure_dir_existance(dirname)
        edge_file = os.path.join(dirname, "edges.txt")
        stats_file = os.path.join(dirname, "stats.txt")
        init_file = os.path.join(dirname, "init.txt")
        reads_file = os.path.join(dirname, "reads.txt")
        contigs_file = os.path.join(dirname, "contigs.fasta")
        self.printStats(stats_file)
        self.printEdges(edge_file)
        self.printInit(init_file)
        self.printReads(reads_file)
        self.printContigs(contigs_file)

    def printEdges(self, fname):
        file = open(fname, "w")
        f = TokenWriter(file)
        f.writeIntLine(self.component.e.__len__())
        for eid in self.component.e:
            f.writeToken(eid)
        f.newLine()
        f.writeIntLine(self.unique.__len__())
        for eid in self.unique:
            f.writeToken(eid)
        f.newLine()
        file.close()


class AlignmentDumpParser:
    def __init__(self, fname):
        self.fname = fname

    def parse(self):
        # type: () -> Generator[Tuple[str, List[str]]]
        edges = []
        rid = None
        for s in open(self.fname, "r").readlines():
            if s.startswith("Chain"):
                if rid is not None:
                    yield rid, edges
                    edges = []
                continue
            s = s.split()
            eid = s[6].split("_")[1]
            rid = s[2][1:]
            if s[6][0] == "-":
                eid = "-" + eid
            edges.append(eid)
        if rid is not None:
            yield rid, edges

def constructComponentRecords(graph, dir, calculator):
    basic.ensure_dir_existance(dir)
    oppa = []
    simple = 0
    max_cov = graph.covPerc(0.5)
    edgecomp = dict()  # type: Dict[str, List[int]]
    for eid in graph.e:
        edgecomp[eid] = []
    componentRecords = []  # type: List[ComponentRecord]
    for comp, cov in SplitGraph(graph, calculator):  # type: SimpleGraph, float
        if len(comp) <= 2:
            continue
        if len(comp) <= 4:
            simple += 1
            oppa.extend(comp.v)
            continue
        if len(comp) > 100:
            print "Complex", comp.__len__()
        print componentRecords.__len__(), len(comp), comp.__len__(), cov, max_cov
        if comp.__len__() <= 100:
            fig_file = os.path.join(dir, str(componentRecords.__len__()) + ".dot")
            f = open(fig_file, "w")
            coloring = lambda v: "white" if len(v.inc) + len(v.out) == len(graph.v[v.id].inc) + len(
                graph.v[v.id].out) else ("yellow" if len(graph.v[v.id].inc) + len(graph.v[v.id].out) < 50 else "red")
            comp.Print(f, coloring, calculator.edgeColoring(cov))
            f.close()
        rec = ComponentRecord(comp, cov)
        for e in comp.e.values():
            if calculator.uniqueCondition(cov)(e):
                rec.addUniqueEdge(e.id)
            else:
                edgecomp[e.id].append(componentRecords.__len__())
        rec.calcStats()
        componentRecords.append(rec)
    f = open(os.path.join(dir, "small.dot"), "w")
    graph.Draw(set(oppa), f)
    f.close()
    return componentRecords, edgecomp


def FillFlyeNext(componentRecords, log_file):
    flye_next = dict()
    for compRec in componentRecords:
        for eid in compRec.unique:
            flye_next[eid] = None
            flye_next[basic.Reverse(eid)] = None
    for s in open(log_file, "r").readlines():
        if s.find("UPath") == -1:
            continue
        s = s[s.find("UPath"):].split()[2::2]
        s = [eid for eid in s if eid in flye_next]
        for i in range(len(s) - 1):
            flye_next[s[i]] = s[i + 1]
            flye_next[basic.Reverse(s[i + 1])] = basic.Reverse(s[i])
    return flye_next


def main(flye_dir, output_dir, diploid):
    basic.ensure_dir_existance(output_dir)
    CreateLog(output_dir)
    graph_file = os.path.join(flye_dir, "20-repeat", "graph_before_rr.gv")
    edge_file = os.path.join(flye_dir, "20-repeat", "graph_before_rr.fasta")
    dump_file = os.path.join(flye_dir, "20-repeat", "read_alignment_dump")
    if diploid:
        calculator = DipolidCalculator(150000)
    else:
        calculator = HaploidCalculator(150000)
    print "Reading graph from", graph_file
    graph = SimpleGraph()
    graph.ReadDot(graph_file)
    print "Reading sequences from", edge_file
    graph.FillSeq(edge_file, True)
    print "Splitting graph", edge_file
    componentRecords, edgecomp = constructComponentRecords(graph, os.path.join(output_dir, "pics"), calculator)
    print "Reading alignment dump from", dump_file
    rcnt = 0
    for rid, eids in AlignmentDumpParser(dump_file).parse():
        compids = set()
        for eid in eids:
            if eid not in edgecomp:
                eid = basic.Normalize(eid)
            for compid in edgecomp[eid]:
                compids.add(compid)
        for compid in compids:
            comp_eids = [eid for eid in eids if eid in componentRecords[compid].component.e]
            if comp_eids.__len__() == 0:
                print "GOPA", compid, compids, rid, eids
            componentRecords[compid].addRead(rid, eids)
        rcnt += 1
        if rcnt % 100000 == 0:
            print "Processed", rcnt, "reads"
    print "Filling flye repeat resolution results"
    flye_next = FillFlyeNext(componentRecords, os.path.join(flye_dir, "flye.log"))
    for compRec in componentRecords:
        half = compRec.half()
        for norm_eid in compRec.unique:
            for eid in [norm_eid, basic.Reverse(norm_eid)]:
                if eid not in compRec.component.e:
                    assert not basic.isCanonocal(eid)
                    assert basic.Reverse(eid) in compRec.component.e
                    continue
                if compRec.component.e[eid].end in half:
                    if compRec.component.isBorder(compRec.component.e[eid].end):
                        compRec.out += 1
                    if compRec.component.isBorder(compRec.component.e[eid].start):
                        compRec.inc += 1
                if not compRec.component.isBorder(compRec.component.e[eid].end):
                    if flye_next[eid] is None:
                        compRec.unresolved_connections += 1
                    else:
                        compRec.resolved_connections.append((eid, flye_next[eid]))
                        if flye_next[eid] not in compRec.component.e:
                            compRec.outside_connections += 1
    basic.ensure_dir_existance(output_dir)
    print "Printing components to disk"
    for i, component in enumerate(componentRecords):
        component.dump(os.path.join(output_dir, str(i)))
    table_file = os.path.join(output_dir, "table.txt")
    print "Printing table to file", table_file
    f = open(table_file, "w")
    f.write("Id v e unique inc out unresolved resolved outside zero hub badborder\n")
    for i, compRec in enumerate(componentRecords):
        comp = compRec.component
        f.write(" ".join([str(i), str(comp.v.__len__()), str(comp.e.__len__()), str(compRec.unique.__len__() * 2),
                          str(compRec.inc), str(compRec.out),
                          str(compRec.unresolved_connections), str(compRec.resolved_connections.__len__()),
                          str(compRec.outside_connections), str(compRec.zero), str(compRec.red), str(compRec.bad_border)]) + "\n")
    f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "diploid" in sys.argv)