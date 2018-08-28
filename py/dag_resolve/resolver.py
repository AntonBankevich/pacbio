import itertools
import sys

from typing import Optional

from common import basic
from dag_resolve import repeat_graph, line_tools, sequences, align_tools, params
from dag_resolve.params import radius


class LineTail:
    def __init__(self, line, edge, tail_consensus, read_collection):
        # type: (line_tools.Line, repeat_graph.Edge, align_tools.Consensus, sequences.ReadCollection) -> LineTail
        self.line = line
        self.edge = edge
        self.tail_consensus = tail_consensus
        self.reads = read_collection
        self.alignment = None # type: align_tools.AlignedSequences

    def __len__(self):
        # type: () -> int
        return self.tail_consensus.__len__()

class VertexResolver:
    def __init__(self, graph, alignerObject):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> VertexResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = alignerObject

    def transitionSupport(self, prev, seg, next):
        # type: (repeat_graph.Edge, line_tools.LineSegment, repeat_graph.Edge) -> sequences.ReadCollection
        return seg.reads.inter(prev.suffix(500)).inter(next.prefix(500))

    def resolveVertex(self, v):
        # type: (repeat_graph.Vertex) -> list[LineTail]
        print "Resolving vertex", v.id
        print "Incoming edges:", map(lambda e: e.id, v.inc)
        print "Outgoing edges:", map(lambda e: e.id, v.out)
        assert len(v.out) > 0, "Trying to resolve hanging vertex"
        res = []
        for edge in v.inc:
            for lineSegment in self.lineStorage.resolved_edges[edge.id]:
                support = []
                for next_candidate in v.out:
                    support.append(self.transitionSupport(edge, lineSegment, next_candidate))
                if len(v.out) == 1:
                    next = 0
                else:
                    supportWeight = map(len, support)
                    print "Number of edges that support transitions for the line from edge", edge.id, ":", ",".join(
                        map(str, supportWeight))
                    next, alternative = basic.best2(supportWeight, lambda x, y: x > y)
                    if supportWeight[alternative] * 5 > supportWeight[next]:
                        print "Support of alternative edge is too high. Aborting."
                        return None
                reads = support[next]# type: sequences.ReadCollection
                # tail = self.aligner.polishAndAnalyse(reads, out_edge.prefix(3000).subcontig())
                tail = self.aligner.polishAndAnalyse(reads, v.out[next])
                res.append(LineTail(lineSegment.line, v.out[next], tail, reads))
                print "Created tail on edge", edge.id, "of length", len(tail.cut())
        return res



class EdgeResolver:
    def __init__(self, graph, aligner):
        # type: (repeat_graph.Graph, align_tools.Aligner) -> EdgeResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = aligner

    def findDivergence(self, tails):
        # type: (list[LineTail]) -> tuple[list[line_tools.Divergence], list[line_tools.Phasing]]
        assert len(tails) > 0
        edge = tails[0].edge
        alignments_list = map(lambda t: t.alignment, tails)
        first = max(map(lambda al: al.first, alignments_list))
        last = min(map(lambda al: al.last, alignments_list))
        assert first + radius < last - radius
        bad_positions = list(basic.merge(*map(lambda alignment: alignment.divPositions(first + radius, last - radius), alignments_list)))
        div_list = []  # type: list[line_tools.Divergence]
        phasings = []  # type: list[line_tools.Phasing]
        for t in tails:
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
            if rec[1] - rec[0] < 2 * radius:
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
                sys.stdout.write("Long divergence: " + str(rec[0]) + " " + str(rec[1]) + " " + str(rec[2]) + "\n")
        print "Divergence rate:", float(sum([pos_rec[2] for pos_rec in div_pos])) / len(edge)
        for i, div in enumerate(div_list):
            print "Divergence", i, ":", div.pos, ",".join([state.seq for state in div.states if state.seq is not None])
        for phasing in phasings:
            phasing.printToFile(sys.stdout)
        return div_list, phasings

    def determinePhasing(self, read, consensus, div_list):
        # type: (sequences.Read, sequences.Contig, list[line_tools.Divergence]) -> line_tools.Phasing
        # print "Determining phasing of read", read.id
        seq = read.seq
        res = line_tools.Phasing()
        alignment = self.aligner.ReadToAlignedSequences(read, consensus)
        # alignment = self.aligner.matchingAlignment([read.seq], consensus)[0]
        # print "from-to:", read.__str__()
        if not alignment.aligned:
            print "WARNING: ALIGNED READ DID NOT ALIGN:", read.id
            return line_tools.Phasing()
        seq = alignment.seq_from
        cnt = 0
        toprint = []
        for div in div_list:
            cnt += 1
            pos = alignment.mapSegment(div.pos[0], div.pos[1])
            if pos is None:
                toprint.append("*")
                continue
            nighborhood = seq[pos[0] - radius // 2: pos[1] + radius // 2 + 1]
            weights = [align_tools.AccurateAligner().align(nighborhood, state.seq) for state in div.states]
            # print weights
            sorted_weights = sorted(weights)
            if sorted_weights[0] * 1.3 > sorted_weights[1]:
                res.add(div.ambiguous)
            else:
                for i, w in enumerate(weights):
                    if w == sorted_weights[0]:
                        res.add(div.states[i])
                        break
            toprint.append(res[-1].char())
        sys.stdout.write("".join(toprint[::-1]) + "\n")
        return res

    def callRead(self, read, consensus, div_list, phasings, min_div_number = 5):
        # type: (sequences.Read, sequences.Contig, list[line_tools.Divergence], list[line_tools.Phasing]) -> Tuple[Optional[int], bool]
        print "Calling read", read.id
        read_phasing = self.determinePhasing(read, consensus, div_list)
        # for phasing in phasings:
        #     phasing.printToFile(sys.stdout)
        dists = []
        phasings_diff = self.comparePhasings(read_phasing, phasings)
        for diff in phasings_diff:
            called_div_number = len(diff) - diff.count(-1)
            if called_div_number < min_div_number:
                dists.append(1)
            else:
                dists.append(float(diff.count(0)) / called_div_number)
        best = basic.best2(dists)
        # print map(lambda diff: (diff.count(0), diff.count(1), diff.count(-1)), phasings_diff)
        # print dists, best
        if dists[best[0]] > 0.2 or dists[best[1]] < 0.2 or dists[best[1]] < dists[best[0]] * 1.5:
            print "Fail"
            return None, False
        for phase in read_phasing:
            phase.divergence.statistics.called += 1
            if phase.isAmbiguous():
                phase.divergence.statistics.ambig += 1
        if read_phasing.called() >= 8 and dists[best[0]] < 0.1 and dists[best[0]] < 1. / len(phasings) and dists[best[1]] > 0.5:
            print "Call full success:", best[0]
            for i, val in enumerate(phasings_diff[best[0]]):
                read_phasing[i].divergence.statistics.called += 1
                if val == 1:
                    read_phasing[i].divergence.statistics.correct += 1
                elif val == 0:
                    read_phasing[i].divergence.statistics.wrong += 1
            return best[0], True
        else:
            print "Call conditional success:", best[0]
            for i, val in enumerate(phasings_diff[best[0]]):
                read_phasing[i].divergence.statistics.called += 1
                if val == 1:
                    read_phasing[i].divergence.statistics.correct += 1
                elif val == 0:
                    read_phasing[i].divergence.statistics.wrong += 1
            return best[0], False


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

    def alignAndCutTails(self, tails):
        # type: (list[LineTail]) -> int
        print "New tails lengths:", map(len, tails)
        for tail in tails:
            tail.tail_consensus = tail.tail_consensus.cut()
        print "Cut tails lengths:", map(len, tails)
        alignments_list = self.aligner.matchingAlignment(map(lambda tail: tail.tail_consensus.seq, tails), tails[0].edge)
        for tail, alignment in zip(tails, alignments_list):
            tail.alignment = alignment
        print "Tail alignments:", map(lambda tail: (tail.alignment.first, tail.alignment.last), tails)
        print "Tail aligned segments:", map(lambda tail: (tail.alignment[tail.alignment.first], tail.alignment[tail.alignment.last]), tails)
        cur_depth = min([tail.alignment.last for tail in tails])
        for tail in tails:
            tail_depth = tail.alignment.alignment[tail.alignment.findPreviousMatch(cur_depth)]
            tail.tail_consensus = tail.tail_consensus.cut(length = tail_depth)
        return cur_depth

    def constructPolishingBase(self, reads, tail):
        # type: (sequences.ReadCollection, LineTail) -> sequences.Contig
        # if len(tail.tail_consensus.full_seq) > len(tail) + 300:
        #     return sequences.Contig(tail.tail_consensus.full_seq, "old_consensus")

        best = tail.alignment.seq_from[:tail.alignment[tail.alignment.last]] + tail.alignment.seq_to[tail.alignment.last:]
        return sequences.Contig(best, "Auto")

        best_match = 0
        best = ""
        best_id = None
        for read in tail.reads:
            for alignment in read.alignments:
                if alignment.seg_to.contig.id != tail.edge.id or alignment.rc:
                    continue
                if alignment.seg_to.right > tail.tail_consensus.__len__() - 500 and \
                                                read.__len__() - alignment.seg_from.right + alignment.seg_to.right > tail.alignment.last + 3000:
                    pos = tail.alignment.findPreviousMatch(alignment.seg_to.right)
                    best = tail.alignment.seq_from[:tail.alignment[pos]] + \
                           read.seq[alignment.seg_from.right - alignment.seg_to.right + pos:]
                    best_match = read.__len__() - alignment.seg_from.right + alignment.seg_to.right
                    best_id = read.id
                    print "Best full:", read.__str__()
                    print "Final best:", best_id
                    return sequences.Contig(best, best_id)
        if best_id is not None and best_match > tail.alignment.last + 2000:
            print "Final best:", best_id
            return sequences.Contig(best, best_id)
        best_id = None
        best = ""
        best_match = 0
        for read in reads:
            for alignment in read.alignments:
                if alignment.seg_to.contig.id == tail.edge.id:
                    if not alignment.rc:
                        if alignment.seg_to.left > 500 and alignment.seg_from.left > 300:
                            continue
                        if read.__len__() - alignment.seg_from.right + alignment.seg_to.right - tail.alignment.last > 2000 and \
                                        alignment.seg_to.right > tail.alignment.last - 500 and \
                                        alignment.seg_from.right - alignment.seg_from.left > best_match:
                            pos = tail.alignment.findPreviousMatch(alignment.seg_to.right)
                            best = tail.alignment.seq_from[:tail.alignment[pos]] + \
                                   read.seq[alignment.seg_from.right - alignment.seg_to.right + pos:]
                            best_match = alignment.seg_from.right
                            best_id = read.id
                            print "Best:", read.__str__()
        if best_id is None:
            print sequences.Contig("No next read")
            return tail.edge.seq
        print "Final best:", best_id
        # best = tail.alignment.seq_from[:tail.alignment[tail.alignment.last]] + tail.alignment.seq_to[tail.alignment.last:]
        return sequences.Contig(best, best_id)

#This method is not used
    def AlignToTails(self, tails, reads):
        # type: (list[LineTail], sequences.ReadCollection) -> sequences.ReadCollection
        tail_contigs = sequences.ContigCollection()
        for i, tail in enumerate(tails):
            tail_contig = sequences.Contig(tail.tail_consensus.seq, i)
            tail_contigs.add(tail_contig)
        res = sequences.ReadCollection(tail_contigs)
        for tail_contig in tail_contigs:
            res.loadFromSam(self.aligner.align(reads, sequences.ContigCollection([tail_contig])))
        return res

    def resolveEdge(self, e, tails):
        # type: (repeat_graph.Edge, list[LineTail]) -> bool
        tails = list(tails)
        for tail in tails:
            tail.tail_consensus.printQuality(sys.stdout)
        phased_reads = sequences.ReadCollection(sequences.ContigCollection([e]))
        for tail in tails:
            phased_reads.extend(tail.reads)
        undecided = e.reads.minus(phased_reads)
        print len(e.reads), len(phased_reads), len(undecided)
        cur_depth = self.alignAndCutTails(tails)
        # edge_consensus = self.aligner.CalculateConsensusCoverage(e, e.reads)
        # edge_consensus.printCoverage(sys.stdout, 50)
        min_div_number = 5
        while True:
            print "Starting new iteration with current depth:", cur_depth
            print len(e.reads) - len(undecided), "reads already classified.", len(undecided), "reads left to classify"
            line_reads = []  # type: list[sequences.ReadCollection]
            for tail in tails:
                line_reads.append(sequences.ReadCollection(sequences.ContigCollection([e])))
            print "Tail lengths: ", map(lambda tail: len(tail.tail_consensus), tails)
            divergences, phasings = self.findDivergence(tails)
            cnt_end_div = 0
            for div in divergences:
                if div.pos[0] > cur_depth - 3000:
                    cnt_end_div += 1
            print "End-div:", cnt_end_div
            # if cnt_end_div < 5:
            #     print "Relaxing"
            #     min_div_number = 2
            #     params.reliable_coverage = 7
            if len(divergences) >= 5:
                interesting_reads = list(undecided.inter(e.prefix(cur_depth)).reads.values())
                interesting_reads = sorted(interesting_reads, key = lambda read: read.contigAlignment(e)[1])
                print "Attempting to classify", len(interesting_reads), "reads"
                # print "Aligning reads to tails"
                # tail_alignment = self.AlignToTails(tails, interesting_reads)
                # print "Alignment finished. Calling reads."
                if cur_depth > len(e) - 1000:
                    interesting_reads = e.reads
                for read in interesting_reads:
                    # if read.id not in tail_alignment.reads:
                    #     print "WARNING: Read", read.id, "was not aligned to tails. Skipping."
                    #     print read
                    # read_to_tail = tail_alignment.reads[read.id]
                    call_result, permanent = self.callRead(read, e, divergences, phasings, min_div_number)
                    if call_result is not None:
                        if permanent:
                            tails[call_result].reads.add(read)
                            line_reads[call_result].add(read)
                            undecided.remove(read)
                        else:
                            line_reads[call_result].add(read)
                    # if call_result is not None:
                    #     line_reads[call_result].add(read)
                for div in divergences:
                    print div.pos, div.statistics.__str__()
                print "Calling results for", len(tails), "tails"
                for tail, reads in zip(tails, line_reads):
                    print tail.line.shortStr(), len(reads), "of", len(interesting_reads)
            for tail, read_collection in zip(tails, line_reads):
                read_collection.extend(tail.reads)
                # tail.tail_consensus = self.aligner.polishAndAnalyse(read_collection, e.prefix(cur_depth + 3000).subcontig())
                tail.tail_consensus = self.aligner.polishAndAnalyse(read_collection, e, self.constructPolishingBase(read_collection, tail))
            for tail in tails:
                tail.tail_consensus.printQuality(sys.stdout)
            new_depth = self.alignAndCutTails(tails)
            print "New depth", new_depth
            if new_depth < cur_depth + 100:
                if new_depth > len(e) - 200:
                    for tail, read_collection, phasing in zip(tails, line_reads, phasings):
                        self.lineStorage.extendRight(tail.line, e, tail.tail_consensus, phasing, read_collection)
                    return True
                return False
            cur_depth = new_depth


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
        # type: (repeat_graph.Vertex) -> None
        tails = sorted(self.vertexResolver.resolveVertex(v), key = lambda tail: tail.edge.id)
        for edge, tails_generator in itertools.groupby(tails, lambda tail: tail.edge):
            tails_list = list(tails_generator) # type: list[LineTail]
            if edge.id in self.lineStorage.resolved_edges:
                assert len(tails_list) == 1
                tails_list[0].line.extendRight(edge, edge.seq, line_tools.Phasing(None), edge.reads)
                print "Successfully connected line", tails_list[0].line.shortStr()
            else:
                finished = self.edgeResolver.resolveEdge(edge, tails_list)
                if finished:
                    print "Successfully resolved edge", edge.id, "into", len(tails_list), "lines"
                else:
                    print "Failed to resolve edge", edge.id

    def resolve(self):
        print self.graph.V
        potentially_resolvable = set([v.id for v in self.graph.V.values() if v.id not in [self.graph.source.id, self.graph.sink.id]])
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


