import sys

from common import basic
from dag_resolve import repeat_graph, line_tools, sequences, align_tools
import itertools
from typing import Generator

radius = 10

class LineTail:
    def __init__(self, line, edge, tail_consensus, read_collection):
        # type: (line_tools.Line, repeat_graph.Edge, align_tools.Consensus, sequences.ReadCollection) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus
        self.reads = read_collection

class VertexResolver:
    def __init__(self, graph, alignerObject):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> VertexResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = alignerObject

    def resolveVertex(self, v):
        # type: (repeat_graph.Vertex) -> Generator[LineTail]
        print "Resolving vertex", v.id
        print "Incoming edges:", map(lambda e: e.id, v.inc)
        print "Outgoing edges:", map(lambda e: e.id, v.out)
        if len(v.out) > 1:
            assert False, "Vertex resolution is not implemented"
        out_edge = v.out[0]
        res = []
        for edge in v.inc:
            for lineSegment in self.lineStorage.resolved_edges[edge.id]:
                reads = lineSegment.reads.inter(edge.suffix(500)).inter(out_edge.prefix(500))# type: sequences.ReadCollection
                tail = self.aligner.polishAndAnalyse(reads, out_edge.prefix(3000).subcontig())
                res.append(LineTail(lineSegment.line, out_edge, tail, reads))
                print "Created tail on edge", edge.id, "of length", len(tail.cut())
        return res.__iter__()



class EdgeResolver:
    def __init__(self, graph, aligner):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> EdgeResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = aligner

    def findDivergence(self, tails, edge):
        # type: (list[str], repeat_graph.Edge) -> tuple[list[line_tools.Divergence], list[line_tools.Phasing]]
        alignments_list = self.aligner.matchingAlignment(tails, edge)
        first = max(map(lambda al: al.first, alignments_list))
        last = min(map(lambda al: al.last, alignments_list))
        assert first + 20 < last
        bad_positions = list(basic.merge(*map(lambda alignment: alignment.divPositions(first, last), alignments_list)))
        div_list = []  # type: list[line_tools.Divergence]
        phasings = []  # type: list[line_tools.Phasing]
        for alignment in alignments_list:
            phasings.append(line_tools.Phasing())
        if len(bad_positions) == 0:
            return div_list, phasings
        div_pos = [[bad_positions[0], bad_positions[0], 1]]
        for pos in bad_positions[1:]:
            if pos > div_pos[-1][1] + radius:
                div_pos.append([pos, pos, 1])
            else:
                div_pos[-1][1] = pos
                div_pos[-1][2] += 1
        for rec in div_pos:
            if rec[1] - rec[0] < radius:
                divergence = line_tools.Divergence(edge, (rec[0], rec[1]))
                states = []
                for alignment in alignments_list:
                    seg = alignment.mapSegment(rec[0] - radius, rec[1] + radius)
                    seg_seq = alignment.seq_from[seg[0]: seg[1] + 1]
                    states.append(seg_seq)
                if len(set(states)) >= 2:
                    states_dict = dict()
                    for i, seg_seq in enumerate(states):
                        if seg_seq not in states_dict:
                            states_dict[seg_seq] = divergence.addState(seg_seq)
                        phasings[i].add(states_dict[seg_seq])
                    div_list.append(divergence)
            else:
                sys.stderr.write("Long divergence: " + str(rec[0]) + " " + str(rec[1]) + " " + str(rec[2]))
        print "Divergence rate:", float(sum([pos_rec[2] for pos_rec in div_pos])) / len(edge)
        for i, div in enumerate(div_list):
            print "Divergence", i, ":", div.pos, ",".join([state.seq for state in div.states if state.seq is not None])
        for phasing in phasings:
            phasing.printToFile(sys.stdout)
        return div_list, phasings

    def determinePhasing(self, read, consensus, div_list):
        # type: (sequences.Read, sequences.Contig, list[line_tools.Divergence]) -> line_tools.Phasing
        print "Determining phasing of read", read.id
        seq = read.seq
        res = line_tools.Phasing()
        alignment = self.aligner.matchingAlignment([seq], consensus)[0]
        seq = alignment.seq_from
        cnt = 0
        for div in div_list:
            cnt += 1
            pos = alignment.mapSegment(div.pos[0], div.pos[1])
            if pos is None:
                sys.stdout.write("*")
                continue
            nighborhood = seq[pos[0] - radius // 2: pos[1] + radius // 2 + 1]
            weights = [align_tools.AccurateAligner().align(nighborhood, state.seq) for state in div.states]
            sorted_weights = sorted(weights)
            if sorted_weights[0] * 1.3 > sorted_weights[1]:
                res.add(div.ambiguous)
            else:
                for i, w in enumerate(weights):
                    if w == sorted_weights[0]:
                        res.add(div.states[i])
                        break
            sys.stdout.write(res[-1].char())
        sys.stdout.write("\n")
        return res

    def callRead(self, read, consensus, div_list, phasings):
        # type: (sequences.Read, sequences.Contig, list[line_tools.Divergence], list[line_tools.Phasing]) -> int
        seq = read.seq
        print "Calling read", read.id
        read_phasing = self.determinePhasing(read, consensus, div_list)
        for phasing in phasings:
            phasing.printToFile(sys.stdout)
        print "Covered divergences:", read_phasing.__len__(), "of", div_list.__len__()
        dists = []
        phasings_diff = self.comparePhasings(read_phasing, phasings)
        for diff in phasings_diff:
            called_div_number = len(diff) - diff.count(-1)
            if called_div_number < 5:
                dists.append(1)
            else:
                dists.append(float(diff.count(0)) / (called_div_number))
        best = basic.smallest2(dists)
        print map(lambda diff: (diff.count(0), diff.count(1), diff.count(-1)), phasings_diff)
        if dists[best[0]] > 0.2 or dists[best[1]] < 0.2 or dists[best[1]] < dists[best[0]] * 1.5:
            print "Fail"
            return None
        print "Call success:", best[0]
        return best[0]


    def comparePhasings(self, query, targets):
        # type: (line_tools.Phasing, list[line_tools.Phasing]) -> list[list[int]]
        res = []
        for i in range(len(targets)):
            res.append([])
        cur = 0
        for query_num, query_state in enumerate(query.states):
            while cur < len(targets[0]) and targets[0][cur].divergence != query_state.divergence:
                cur += 1
            if cur == len(targets[0]):
                break
            for i, target in enumerate(targets):
                if query_state.isAmbiguous():
                    res[i].append(-1)
                elif query_state == target[cur]:
                    res[i].append(1)
                else:
                    res[i].append(0)
        return res

    def resolveEdge(self, e, tails):
        # type: (repeat_graph.Edge, list[LineTail]) -> None
        tails = list(tails)
        phased_reads = sequences.ReadCollection(sequences.ContigCollection([e]))
        for tail in tails:
            phased_reads.extend(tail.reads)
        cur_depth = min([len(lt.tail_consensus.cut()) for lt in tails])
        while True:
            line_reads = []  # type: list[sequences.ReadCollection]
            for tail in tails:
                line_reads.append(sequences.ReadCollection(sequences.ContigCollection([e])))
                tail.tail_consensus.printQuality(sys.stdout)
            print "Current depth:", cur_depth
            divergences, phasings = self.findDivergence([a.tail_consensus.seq[:cur_depth] for a in tails], e)
            for read in e.reads.inter(e.prefix(cur_depth)):
                if read.id in phased_reads.reads:
                    print "Read", read.id, "is already phased, skipping"
                call_result = self.callRead(read, e, divergences, phasings)
                if call_result is not None:
                    line_reads[call_result].add(read)
            print "Calling results for", len(tails), "tails"
            for tail, reads in zip(tails, line_reads):
                print tail.line.shortStr(), len(reads)
            for tail, read_collection in zip(tails, line_reads):
                read_collection.extend(tail.reads)
                tail.tail_consensus = self.aligner.polishAndAnalyse(read_collection, e.prefix(cur_depth + 3000).subcontig())
                tail.tail_consensus.printQuality(sys.stderr)
            new_depth = min([len(lt.tail_consensus.cut()) for lt in tails])
            cur_depth = new_depth
            if new_depth > len(e) - 500:
                for tail, line_reads in zip(tails, line_reads):
                    self.lineStorage.extendRight(tail.line, e, tail.tail_consensus, line_reads)
                return True
            if new_depth < cur_depth + 500:
                return False



class GraphResolver:
    def __init__(self, graph, vertexResolver, edgeResolver):
        # type: (repeat_graph.Graph, VertexResolver, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.lineStorage = line_tools.LineStorage(graph)
        self.vertexResolver = vertexResolver
        vertexResolver.lineStorage = self.lineStorage
        self.edgeResolver = edgeResolver
        edgeResolver.lineStorage = self.lineStorage

    def resolveVertexForward(self, v):
        tails = sorted(self.vertexResolver.resolveVertex(v), key = lambda tail: tail.edge.id)
        for edge, tails_generator in itertools.groupby(tails, lambda tail: tail.edge):
            tails_list = list(tails_generator)
            finished = self.edgeResolver.resolveEdge(edge, tails_list)
            if finished:
                print "Successfully resolved edge", edge.id, "into", len(tails_list), "lines"
            else:
                print "Failed to resolve edge", edge.id

    def resolve(self):
        print self.graph.V
        potentially_resolvable = set([v1.id for v1 in self.graph.V.values() if len(v1.out) == 1])
        while True:
            unresolvable = set()
            cnt = 0
            for v_id in potentially_resolvable:
                v = self.graph.V[v_id]
                if self.lineStorage.isResolvableLeft(v):
                    self.resolveVertexForward(v)
                    unresolvable.add(v.id)
                    cnt += 1
            potentially_resolvable = potentially_resolvable.difference(unresolvable)
            if cnt == 0:
                break


