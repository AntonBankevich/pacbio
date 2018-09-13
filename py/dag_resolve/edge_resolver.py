import sys

from typing import Optional, Tuple

import dag_resolve.phasing
from alignment import align_tools
from alignment.accurate_alignment import AccurateAligner
from common import basic, SeqIO
from dag_resolve import line_tools, sequences, params, repeat_graph
from dag_resolve.divergence import Divergence
from dag_resolve.line_tools import LineTail
from dag_resolve.params import radius


class EdgeResolver:
    def __init__(self, graph, aligner, polisher):
        # type: (repeat_graph.Graph, align_tools.Aligner, align_tools.Polisher) -> EdgeResolver
        self.graph = graph
        self.lineStorage = None
        self.aligner = aligner
        self.polisher = polisher

    def findDivergence(self, tails):
        # type: (list[LineTail]) -> tuple[list[dag_resolve.divergence.Divergence], list[dag_resolve.phasing.Phasing]]
        assert len(tails) > 0
        edge = tails[0].edge
        alignments_list = map(lambda t: t.alignment, tails)
        first = max(map(lambda al: al.first, alignments_list))
        last = min(map(lambda al: al.last, alignments_list))
        assert first + radius < last - radius
        bad_positions = list(basic.merge(*map(lambda alignment: alignment.divPositions(first + radius, last - radius), alignments_list)))
        div_list = []  # type: list[dag_resolve.divergence.Divergence]
        phasings = []  # type: list[dag_resolve.phasing.Phasing]
        for t in tails:
            phasings.append(dag_resolve.phasing.Phasing())
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
            if rec[1] - rec[0] < 3 * radius:
                divergence = Divergence(edge, (rec[0], rec[1]))
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
        # type: (sequences.AlignedRead, sequences.Contig, list[dag_resolve.divergence.Divergence]) -> line_tools.Phasing
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
                res.add(div.ambiguous)
                continue
            div.statistics.called += 1
            nighborhood = seq[pos[0] - radius // 2: pos[1] + radius // 2 + 1]
            weights = [AccurateAligner().align(nighborhood, state.seq) for state in div.states]
            # print weights
            sorted_weights = sorted(weights)
            if sorted_weights[0] * 1.3 > sorted_weights[1]:
                res.add(div.ambiguous)
                div.statistics.ambig += 1
            else:
                for i, w in enumerate(weights):
                    if w == sorted_weights[0]:
                        res.add(div.states[i])
                        break
            toprint.append(res[-1].char())
        sys.stdout.write("".join(toprint[::-1]) + "\n")
        return res

    def callRead(self, read, consensus, div_list, phasings, min_div_number = 5):
        # type: (sequences.AlignedRead, sequences.Contig, list[dag_resolve.divergence.Divergence], list[line_tools.Phasing]) -> tuple[Optional[int], bool]
        print "Calling read", read.__str__()
        read_phasing = self.determinePhasing(read, consensus, div_list)
        if len(phasings) > 2:
            for phasing in phasings:
                print phasing.__str__()[::-1]
        dists = []
        # phasings_diff = self.comparePhasings(read_phasing, phasings)
        called_div_number = read_phasing.called()
        if called_div_number < min_div_number:
            print "Fail. Not enough called positions"
            return None, False
        for phasing in phasings:
            dists.append(float(phasing.diff(read_phasing)))
        best = basic.best2(dists)
        if dists[best[0]] > 0.2 * called_div_number:
            print "Fail. No similar copy"
            return None, False
        d1, d2, diff = read_phasing.diff2(phasings[best[0]], phasings[best[1]])
        if diff < 3 or d1 > 0.2 * diff or d2 < 0.8* diff:
            print "Fail. Too similar copies: ", d1, d2, diff
            return None, False
        if read_phasing.called() >= 8 and diff >= 5 and dists[best[0]] < 0.1 * called_div_number \
                and dists[best[0]] < 1. / len(phasings) * called_div_number:
            print "Call full success:", best[0]
            for ph_read, ph_copy in zip(read_phasing.__iter__(), phasings[best[0]]):
                assert ph_read.divergence == ph_copy.divergence
                if ph_read.isAmbiguous():
                    continue
                if ph_read == ph_copy:
                    ph_read.divergence.statistics.correct += 1
                else:
                    ph_read.divergence.statistics.wrong += 1
            return best[0], True
        else:
            print "Call conditional success:", best[0]
            for ph_read, ph_copy in zip(read_phasing.__iter__(), phasings[best[0]]):
                assert ph_read.divergence == ph_copy.divergence
                if ph_read.isAmbiguous():
                    continue
                if ph_read == ph_copy:
                    ph_read.divergence.statistics.correct += 1
                else:
                    ph_read.divergence.statistics.wrong += 1
            return best[0], False


    def comparePhasings(self, query, targets):
        # type: (dag_resolve.phasing.Phasing, list[dag_resolve.phasing.Phasing]) -> list[list[int]]
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
        self.alignTails(tails)
        cur_depth = min([tail.alignment.last for tail in tails])
        for tail in tails:
            tail_depth = tail.alignment.alignment[tail.alignment.findPreviousMatch(cur_depth)]
            tail.tail_consensus = tail.tail_consensus.cut(length = tail_depth)
        return cur_depth

    def alignTails(self, tails):
        # type: (list[LineTail]) -> None
        print "New tails lengths:", map(len, tails)
        for tail in tails:
            tail.tail_consensus = tail.tail_consensus.cut()
        print "Cut tails lengths:", map(len, tails)
        alignments_list = self.aligner.matchingAlignment(map(lambda tail: tail.tail_consensus.seq, tails),
                                                         tails[0].edge)
        for tail, alignment in zip(tails, alignments_list):
            tail.alignment = alignment
        print "Tail alignments:", map(lambda tail: (tail.alignment.first, tail.alignment.last), tails)
        print "Tail aligned segments:", map(
            lambda tail: (tail.alignment[tail.alignment.first], tail.alignment[tail.alignment.last]), tails)

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
        # type: (repeat_graph.Edge, list[LineTail]) -> tuple[bool, Optional[repeat_graph.Edge]]
        tails = list(tails)
        if len(tails) == 1:
            tails[0].reads.extend(e.reads)
            tails[0].tail_consensus = self.polisher.polishQuiver(tails[0].reads, tails[0].tail_consensus.seq, 0)
        for tail in tails:
            tail.tail_consensus.printQuality(sys.stdout)
        phased_reads = sequences.ReadCollection(sequences.ContigCollection([e]))
        for tail in tails:
            phased_reads.extend(tail.reads)
        undecided = e.reads.minus(phased_reads)
        print len(e.reads), len(phased_reads), len(undecided)
        cur_depth = self.alignAndCutTails(tails)
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
            if len(divergences) >= 5:
                interesting_reads = list(undecided.inter(e.prefix(cur_depth)).reads.values())
                if params.full_stats and cur_depth > len(e) - 200:
                    interesting_reads = list(e.reads.reads.values())
                interesting_reads = sorted(interesting_reads, key = lambda read: read.contigAlignment(e)[1])
                print "Attempting to classify", len(interesting_reads), "reads"
                # print "Aligning reads to tails"
                # tail_alignment = self.AlignToTails(tails, interesting_reads)
                # print "Alignment finished. Calling reads."
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
                new_consensus = self.polisher.polishQuiver(read_collection, tail.tail_consensus.seq, 0)
                if new_consensus is not None:
                    print "Old consensus:", len(tail.tail_consensus), "New consensus:", len(new_consensus.cut())
                    tail.tail_consensus = new_consensus
                else:
                    print "Failed to construct longer consensus"
            new_depth = self.alignAndCutTails(tails)
            for tail in tails:
                print "New tail alignes to segment from ", tail.alignment.first, "to", tail.alignment.last, "with percent identity", tail.alignment.percentIdentity()
            for tail in tails:
                tail.tail_consensus.printQuality(sys.stdout)
            print "New depth:", new_depth, "Old depth:", cur_depth
            if new_depth < cur_depth + 100:
                for tail, read_collection in zip(tails, line_reads):
                    tail.reads = read_collection
                if new_depth > len(e) - 200:
                    for tail, phasing in zip(tails, phasings):
                        self.lineStorage.extendRight(tail.line, e, tail.tail_consensus, phasing, tail.reads)
                    return True, None
                repaired, new_edge = self.AttemptRepairGraph(tails)
                if repaired:
                    return False, new_edge
                repaired = min(map(lambda tail: tail.alignment.percentIdentity(), tails)) < 0.85 and self.AttemptRepairRecruitment(tails)
                if not repaired:
                    return False, None
            cur_depth = new_depth

    def AttemptRepairGraph(self, tails):
        # type: (list[LineTail]) -> Tuple[bool, Optional[repeat_graph.Edge]]
        print "Attempting to repair graph"
        for tail in tails:
            tail.tail_consensus = self.polisher.polishQuiver(tail.reads, tail.tail_consensus.seq, 0)
        self.alignTails(tails)
        largest_problem = tails[0]
        print "Tail alignments and lengths"
        for tail in tails[1:]:
            if tail.alignment.last < largest_problem.alignment.last:
                largest_problem = tail
        print "Largest problem is tail of line", largest_problem.line.id, "Alignment is broken at position", largest_problem.alignment.last
        largest_problem.tail_consensus = largest_problem.tail_consensus.full_seq.cut(5)
        if len(largest_problem) - largest_problem.alignment[largest_problem.alignment.last] > 700:
            next_alignment = None
            tmp = sequences.ReadCollection(self.graph.edgeCollection())
            tmp.loadFromSam(self.aligner.align([SeqIO.SeqRecord(largest_problem.tail_consensus.seq, "tail")], self.graph.edgeCollection()))
            read = tmp.reads["tail"]
            print "Largest problem alignment:", read.__str__()
            for al in read.alignments:
                if not al.rc and al.seg_from.right - al.seg_from.left > 500 \
                        and al.seg_from.left > largest_problem.alignment[largest_problem.alignment.last] - 50\
                        and (next_alignment is None or next_alignment.seg_from.left > al.seg_from.left):
                    next_alignment = al
            print "Closest alignment:", next_alignment.__str__()
            if next_alignment is None:
                print "Failed to repair graph. No alignment to known edges of the graph."
                return False, None
            new_edge = self.graph.addCuttingEdge(largest_problem.edge, largest_problem.alignment.last,
                                      self.graph.E[next_alignment.seg_to.contig.id], next_alignment.seg_to.left,
                                      largest_problem.tail_consensus.seq[largest_problem.alignment[largest_problem.alignment.last]:
                                      next_alignment.seg_from.left])
            self.aligner.repairGraphAlignments(self.graph)
            for tail in tails:
                tail.edge = new_edge.start.inc[0]
                tail.alignment = None
            return True, new_edge.start.inc[0]
        else:
            print "Failed to repair graph"
            return False, None

    def AttemptRepairRecruitment(self, tails):
        # type: (list[line_tools.LineTail]) -> bool
        print "Attempting to repair read recruitment"
        edge = tails[0].edge
        new_alignments = dict()
        for i, tail in enumerate(tails):
            for rec in self.aligner.align(self.graph.reads, [sequences.Contig(tail.tail_consensus.seq, 0)]):
                if rec.is_unmapped or rec.alen < 500:
                    continue
                if rec.query_name in edge.reads:
                    continue
                if rec.query_name not in new_alignments:
                    new_alignments[rec.query_name] = []
                new_alignments[rec.query_name].append(i)
        cnt = 0
        for read_id in new_alignments:
            if len(new_alignments[read_id]) == 1:
                tails[new_alignments[read_id][0]].reads.add(self.graph.reads[read_id])
                edge.reads.add(self.graph.reads[read_id])
                print "Added read", self.graph.reads[read_id]
                cnt += 1
        print "Recovered", cnt, "reads"
        return cnt > 0