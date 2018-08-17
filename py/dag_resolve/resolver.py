import sys

from common import basic
from dag_resolve import repeat_graph, line, sequences, align_tools
import itertools
from typing import Generator

radius = 10

class LineTail:
    def __init__(self, line, edge, tail_consensus):
        # type: (line.Line, repeat_graph.Edge, align_tools.Consensus) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus

class VertexResolver:
    def __init__(self, graph, alignerObject):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> VertexResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = alignerObject

    def resolveVertex(self, v):
        # type: (repeat_graph.Vertex) -> Generator[LineTail]
        if len(v.out) > 1:
            assert False, "Vertex resolution is not implemented"
        out_edge = v.out[0]
        for edge in v.inc:
            for lineSegment in self.lineStorage.resolved_edges[edge.id]:
                reads = lineSegment.reads.inter(edge.suffix(500)).inter(out_edge.prefix(500))# type: sequences.ReadCollection
                tail = self.aligner.polishAndAnalyse(reads, out_edge.prefix(3000).subcontig())
                tail.printQuality(sys.stderr)
                yield LineTail(lineSegment.line, out_edge, tail)



class EdgeResolver:
    def __init__(self, graph, aligner):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> EdgeResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = aligner

    def findDivergence(self, tails, edge):
        # type: (list[str], sequences.Contig) -> tuple[list[line.Divergence], list[line.Phasing]]
        alignments_list = self.aligner.matchingAlignment(tails, edge)
        bad_positions = []
        first = max(map(lambda al: al.first(), alignments_list))
        last = min(map(lambda al: al.last(), alignments_list))
        assert first + 20 < last
        for i in xrange(first + radius, last - radius):
            for alignment in alignments_list:
                if alignment.alignment[i] == -1 or \
                        (alignment.alignment[i - 1] != -1 and alignment.alignment[i] != alignment.alignment[i - 1] + 1) or \
                                alignment.seq_to[i] != alignment.seq_from[alignment.alignment[i]]:
                    bad_positions.append(i)
                    break
        div_list = []
        phasings = []  # type: list[line.Phasing]
        for alignment in alignments_list:
            phasings.append(line.Phasing())
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
                divergence = line.Divergence(edge, (rec[0], rec[1]))
                div_list.append(divergence)
                states = dict()
                for i, alignment in enumerate(alignments_list):
                    seg = alignment.mapSegment(rec[0] - radius, rec[1] + radius)
                    seg_seq = alignment.seq_from[seg[0]: seg[1] + 1]
                    if seg_seq not in states:
                        states[seg_seq] = divergence.addState(seg_seq)
                    phasings[i].add(states[seg_seq])
            else:
                sys.stderr.write("Long divergence: " + str(rec[0]) + " " + str(rec[1]) + " " + str(rec[2]))
        return div_list, phasings

    def determinePhasing(self, seq, consensus, div_list):
        # type: (str, sequences.Contig, list[line.Divergence]) -> Generator[line.DivergenceState]
        alignment = self.aligner.matchingAlignment([seq], consensus)[0]
        for div in div_list:
            pos = alignment.findSeqPos(div.pos)
            nighborhood = seq[pos - radius // 2: pos + radius // 2 + 1]
            weights = [align_tools.AccurateAligner().align(nighborhood, state.neighborhood) for state in div.states]
            sorted_weights = sorted(weights)
            if sorted_weights[0] > 1.2 * sorted_weights[1]:
                yield div.ambiguous
            else:
                for i, w in enumerate(weights):
                    if w == sorted_weights[0]:
                        yield div.states[i]
                        break

    def callRead(self, seq, consensus, div_list, phasings):
        # type: (str, sequences.Contig, list[line.Divergence], list[line.Phasing]) -> int
        read_phasing = self.determinePhasing(seq, consensus, div_list)
        dists = []
        phasings_diff = self.comparePhasings(read_phasing, phasings)
        for diff in phasings_diff:
            called_div_number = len(diff) - diff.count(-1)
            if called_div_number < 10:
                dists.append(1)
            else:
                dists.append(float(diff.count(0)) / (called_div_number))
        best = basic.smallest2(diff)
        if diff[best[0]] > 0.05 or diff[best[1]] < 0.1:
            return None
        return best[0]


    def comparePhasings(self, query, targets):
        # type: (line.Phasing, list[line.Phasing]) -> list[list[int]]
        res = []
        for i in range(len(targets)):
            res.append([])
        cur = 0
        for state_num, target_state in enumerate(targets[0].states):
            while cur < len(query) and query[cur].divergence != target_state.divergence:
                cur += 1
            if cur == len(query):
                break
            for i, target in enumerate(targets):
                if query[cur].isAmbiguous():
                    res[i].append(-1)
                elif query[cur] == target[state_num]:
                    res[i][0].append(1)
                else:
                    res[i][1].append(0)
        return res

    def resolveEdge(self, e, tails):
        # type: (repeat_graph.Edge, Generator[LineTail]) -> None
        tails_list = list(tails)
        while True:
            line_reads = []  # type: list[sequences.ReadCollection]
            for tail in tails:
                line_reads.append(sequences.ReadCollection(e))
                tail.tail_consensus.printQuality(sys.stdout)
            cur_depth = min([len(lt.tail_consensus.cut()) for lt in tails_list])
            print "Current depth:", cur_depth
            divergences, phasings = self.findDivergence(tails, e)
            for read in e.reads.inter(e.prefix(cur_depth)):
                call_result = self.callRead(read.seq, e, divergences, phasings)
                if call_result is not None:
                    line_reads[call_result].add(read)
            for tail, read_collection in zip(tails, line_reads):
                tail.tail_consensus = self.aligner.polishAndAnalyse(read_collection, e.prefix(cur_depth + 3000).subcontig())
                tail.printQuality(sys.stderr)


class GraphResolver:
    def __init__(self, graph, vertexResolver, edgeResolver):
        # type: (repeat_graph.Graph, VertexResolver, EdgeResolver) -> GraphResolver
        self.graph = graph
        self.lineStorage = line.LineStorage(graph)
        self.vertexResolver = vertexResolver
        vertexResolver.lineStorage = self.lineStorage
        self.edgeResolver = edgeResolver
        edgeResolver.lineStorage = self.lineStorage

    def resolveVertexForward(self, v):
        tails = sorted(self.vertexResolver.resolveVertex(v), key = lambda tail: tail.edge.id)
        for edge, tails_generator in itertools.groupby(tails, lambda tail: tail.edge):
            self.edgeResolver.resolveEdge(edge, list(tails_generator))

    def resolve(self):
        print self.graph.V
        potentially_resolvable = self.graph.V.values()
        while True:
            cnt = 0
            for v in potentially_resolvable:
                print v.id, "oppa"
                if self.lineStorage.isResolvableLeft(v):
                    self.resolveVertexForward(v)
                    potentially_resolvable.remove(v)
                    cnt += 1
            if cnt == 0:
                break


