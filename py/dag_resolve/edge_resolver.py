import sys
import itertools
from typing import Tuple, Optional, Dict

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common.seq_records import NamedSequence
from dag_resolve import params
from dag_resolve.line_align import Scorer
from dag_resolve.line_tools import Line, LinePosition
from dag_resolve.repeat_graph import Graph, Edge, Vertex
from common.sequences import ReadCollection, ContigCollection, Segment, Contig, AlignmentPiece


class EdgeResolver:
    def __init__(self, graph, aligner, polisher, reads):
        # type: (Graph, Aligner, Polisher, ReadCollection) -> None
        self.graph = graph
        self.aligner = aligner
        self.polisher = polisher
        self.scorer = Scorer()
        self.prolonger = Prolonger(graph, aligner, polisher)
        self.reads = reads

    def resolveVertex(self, v, lines):
        # type: (Vertex, list[Line]) -> list[Optional[Edge]]
        uncertain = ReadCollection()
        for edge in v.inc:
            uncertain.extend(map(lambda read: self.reads[read.id], edge.reads))
        uncertain = uncertain.minusAll([line.reads.inter(line.centerPos.suffix()) for line in lines])
        positions = []
        edges = []
        for line in lines:
            cur = len(line.chain) - 1
            while cur > 0 and line.chain[cur].seg_to.contig not in v.inc:
                cur -= 1
            assert line.chain[cur].seg_to.contig in v.inc
            edges.append(line.chain[cur].seg_to.contig)
            while cur > 0 and line.chain[cur].seg_to.contig == line.chain[cur-1].seg_to.contig and line.chain[-1].precedes(line.chain[cur]):
                cur -= 1
            if cur > 0:
                positions.append(LinePosition(line, max(line.centerPos.pos, line.chain[cur - 1].seg_from.left)))
            else:
                positions.append(line.centerPos)
        classifier = ReadClassifier(self.graph, self.aligner, lines, positions)
        classifier.classifyReads(lines, uncertain)
        res = []
        for edge, line in zip(edges, lines):
            if line.chain[-1].seg_to.contig not in v.out:
                self.prolonger.prolongConsensus(edge, line)
                res.append(self.attemptJump(edge, line))
            else:
                res.append(line.chain[-1].seg_to.contig)
        return res

    def resolveEdge(self, edge, lines):
        # type: (Edge, list[Line]) -> Tuple[bool, Optional[Edge]]
        if len(lines) == 0:
            print "WARNING: EDGE", edge, "HAD ZERO LINES PASSING THROUGH IT!!!"
            return True, None
        print "Resolving edge", edge, "into lines:", ", ".join(map(lambda line: str(line.id), lines))
        positions = []
        for line in lines:
            cur = len(line.chain) - 1
            while cur > 0 and line.chain[cur - 1].seg_to.contig == edge and line.chain[cur - 1].precedes(line.chain[cur]):
                cur -= 1
            if cur > 0:
                positions.append(LinePosition(line, line.chain[cur - 1].seg_from.right))
            else:
                positions.append(LinePosition(line, 0))
        classifier = ReadClassifier(self.graph, self.aligner, lines, positions)
        uncertain = ReadCollection(ContigCollection(lines), map(lambda read: self.reads[read.id], edge.reads))
        print "All reads:"
        edge.reads.print_alignments(sys.stdout)
        print "All on line:"
        uncertain.print_alignments(sys.stdout)
        uncertain = uncertain.minusAll([line.reads.inter(line.centerPos.suffix()) for line in lines])
        print "Uncertain:"
        uncertain.print_alignments(sys.stdout)
        passedLines = []
        for line in lines:
            if self.attemptJump(edge, line) is not None:
                passedLines.append(line)
        # self.prolongAll(edge, lines)
        while True:
            active_lines = [self.shortestLine(lines, passedLines)]
            classified = classifier.classifyReads(active_lines, uncertain)
            uncertain = uncertain.minusBoth(classified) #this should not be important !! Check!!
            total_extention = 0
            for line in active_lines:
                line_extention = self.prolonger.prolongConsensus(edge, line)
                if self.attemptJump(edge, line) is not None:
                    passedLines.append(line)
                total_extention += line_extention
            if len(passedLines) == len(lines):
                classifier.classifyReads(lines, uncertain)
                return True, None
            if len(classified) == 0 and total_extention < 100 and jumped == 0:
                return False, None

    def attemptJump(self, edge, line):
        # type: (Edge, Line) -> Optional[Edge]
        print "Jumping with line", line, "from edge", edge
        shift = line.chain[-1].seg_from.right
        seq = line.seq[shift:]
        if len(seq) < 300:
            print "Tail too short"
            return None
        alignments = ReadCollection(ContigCollection(edge.end.out))
        alignments.addNewRead(NamedSequence(seq, "tail"))
        self.aligner.alignReadCollection(alignments)
        best = None
        print alignments.reads["tail"]
        for al in alignments.reads["tail"].alignments:
            if al.seg_to.left > params.max_jump:
                continue
            if al.seg_to.left < 100 and len(al) > 100 and \
                    line.chain[-1].seg_to.contig == edge and \
                    line.chain[-1].seg_to.right > len(edge) - 100:
                best = al
                continue
            if len(al) < 300 or al.seg_to.contig not in edge.end.out:
                continue
            if al.seg_to.contig.unique():
                if al.seg_from.right < len(al.seg_from.contig) - 100 and al.seg_to.right < len(al.seg_to.contig) - 100:
                    continue
                if al.seg_to.left > 1000:
                    continue
            if best is None or al.seg_from.left < best.seg_from.left:
                best = al
        if best is None:
            print "No jump"
            alignments = ReadCollection(ContigCollection(self.graph.E.values()))
            alignments.addNewRead(NamedSequence(seq, "tail"))
            self.aligner.alignReadCollection(alignments)
            alignments.print_alignments(sys.stdout)
            return None
        line.addAlignment(AlignmentPiece(Segment(line, best.seg_from.left + shift, best.seg_from.right + shift), best.seg_to, best.cigar))
        print "Connected line", line, "to edge", best.seg_to.contig, "using alignment", best
        for al in alignments.reads["tail"].alignments:
            if al.seg_to.contig == best.seg_to.contig and best.precedes(al):
                print "Additional alignment:", al
                line.addAlignment(
                    AlignmentPiece(Segment(line, al.seg_from.left + shift, al.seg_from.right + shift), al.seg_to,
                                   al.cigar))
        return best.seg_to.contig

    # def attemptReattach(self, line):
    #     shift = line.chain[-1].seg_from.right
    #     seq = line.seq[shift:]
    #     if len(seq) < 1000:
    #         return None
    #     alignments = ReadCollection(ContigCollection(edge.end.out))
    #     alignments.addNewRead(NamedSequence(seq, "tail"))
    #     self.aligner.alignReadCollection(alignments)
    #     best = None
    #     print alignments.reads["tail"]
    #     for al in alignments.reads["tail"].alignments:
    #         if len(al) > 300 and (best is None or al.seg_from.left < best.seg_from.left) and al.seg_to.contig in edge.end.out:
    #             best = al
    #     if best is None:
    #         print "No jump"
    #         alignments = ReadCollection(ContigCollection(self.graph.E.values()))
    #         alignments.addNewRead(NamedSequence(seq, "tail"))
    #         self.aligner.alignReadCollection(alignments)
    #         alignments.print_alignments(sys.stdout)
    #         return None
    #     line.addAlignment(AlignmentPiece(Segment(line, best.seg_from.left + shift, best.seg_from.right + shift), best.seg_to, best.cigar))
    #     print "Connected line", line, "to edge", best.seg_to.contig, "using alignment", best
    #     return best.seg_to.contig

    def prolongAll(self, edge, lines):
        # type: (Edge, list[Line]) -> int
        res = 0
        for line in lines:
            tmp = self.prolonger.prolongConsensus(edge, line)
            res += tmp
        return res

    def shortestLine(self, lines, passedLines):
        shortest = None
        for line in lines:
            if line not in passedLines and (shortest is None or line.chain[-1].seg_to.right < shortest.chain[-1].seg_to.right):
                shortest = line
        return shortest


class Prolonger:
    def __init__(self, graph, aligner, polisher):
        # type: (Graph, Aligner, Polisher) -> Prolonger
        self.graph = graph
        self.aligner = aligner
        self.polisher = polisher

    def prolongConsensus(self, edge, line):
        # type: (Edge, Line) -> int
        print "Prolonging line", line
        step_back = min(5000, len(line))
        overlap = min(1000, step_back, len(line) - line.centerPos.pos)
        base_consensus = line.seq[-step_back:]
        # print "inter", line, len(line)
        # line.reads.inter(line.suffix(-step_back)).print_alignments(sys.stdout)
        reads = line.reads.inter(line.suffix(-overlap)).noncontradicting(line.asSegment())

        # print "Filtered reads:"
        # reads.print_alignments(sys.stdout)
        newConsensus = self.polisher.polishQuiver(reads, base_consensus, step_back - overlap).cut()
        if len(newConsensus) <= overlap * 1.2:
            print "Could not find good consensus"
            return 0
        if newConsensus.seq[:overlap] != line.suffix(-overlap).Seq():
            print "Inaccurate glue"
            print newConsensus.seq[:overlap]
            print line.suffix(-overlap).Seq()
            for i in range(overlap):
                if newConsensus.seq[i] != line[-overlap + i]:
                    print "First difference at position", i
                    break
        old_len = len(line)
        line.extendRight(newConsensus, -overlap)
        alignments = ReadCollection(ContigCollection([edge])) # alignments to previous edges may become corrupted!!!
        read = alignments.addNewRead(NamedSequence(newConsensus.seq, "tail"))
        self.aligner.alignReadCollection(alignments)
        read.sort()
        for al in read.alignments:
            if al.seg_to.contig == edge:
                seg_from = Segment(line, al.seg_from.left + old_len - overlap, al.seg_from.right + old_len - overlap)
                new_al = AlignmentPiece(seg_from, al.seg_to, al.cigar)
                print "Candidate alignment:", al, line.chain[-1], line.chain[-1].precedes(new_al)
                if line.chain[-1].seg_to.contig != edge or line.chain[-1].precedes(new_al):
                    line.addAlignment(new_al)
        print "Prolonged line", line, "by", len(line) - old_len, "nucleotides"
        return len(line) - old_len

class ReadClassifier:
    def __init__(self, graph, aligner, lines, positions):
        # type: (Graph, Aligner, list[Line], list[LinePosition]) -> ReadClassifier
        self.graph = graph
        self.aligner = aligner
        self.scorer = Scorer()
        self.lines = lines
        self.positions = dict()# type: Dict[int, LinePosition]
        for pos in positions:
            self.positions[pos.line.id] = pos

    def classifyReads(self, active, reads):
        # type: (list[Line], ReadCollection) -> ReadCollection
        print "Trying to classify reads to lines", map(str, active)
        print "Full list of lines:", map(str, self.lines)
        print "Active lines:", map(str, active)
        line_aligns = self.pairwiseAlign(self.lines)
        alignments = ReadCollection(ContigCollection(self.lines))
        for read in reads.reads.values():
            alignments.addNewRead(read)
        for line in self.lines:
            self.aligner.alignReadsToSegments(alignments, [line.centerPos.suffix()])
        classified = dict()
        for line in self.lines:
            classified[line.id] = []
        n_class = 0
        n_uncertain = 0
        n_irrelevant = 0
        for read in alignments:
            candidates = []
            for al in read.alignments:
                if al.seg_to.contig in self.lines and \
                        al.seg_from.left < 500 and \
                        len(al.seg_from) > min(700, len(al.seg_to.contig.centerPos.suffix()) * 0.8) and \
                        al.seg_to.inter(self.positions[al.seg_to.contig.id].suffix()):
                    candidates.append(al)
            # for al1 in candidates:
            #     for al2 in candidates:
            #         if al1 != al2 and al1.seg_to.contig == al2.seg_to.contig and len(al1.seg_from) < 1000 and len(al2.seg_from) < 1000:
            #             n_irrelevant += 1
            #             continue
            if len(candidates) == 0:
                n_irrelevant += 1
                continue
            print "Classifying read", read, "Initial:", reads[read.id]
            print "Candidates:", map(lambda al: str(al), candidates)
            has_active_choice = False
            for candidate in candidates:
                if candidate.seg_to.contig in active:
                    has_active_choice = True
            if not has_active_choice:
                n_class += 1
                print "No active candidates"
                continue
            res = self.tournament(candidates, line_aligns)
            if res is None:
                n_uncertain += 1
                print "Could not determine the champion"
                continue
            print "Champion:", res, "Active:", (res.seg_to.contig in active)
            n_class += 1
            # if res.seg_to.contig in active:
            if True:
                cread = reads[read.id]
                classified[res.seg_to.contig.id].append(cread)
                line = res.seg_to.contig #type: Line
                line.addRead(cread)
                if res.seg_to.left > res.seg_to.contig.centerPos.pos + 500: # if read does not intercect the center we use the selected alignment
                    seg_from = Segment(cread, res.seg_from.left, res.seg_from.right)
                    cread.addAlignment(AlignmentPiece(seg_from, res.seg_to, res.cigar))
                    print "New read alignments:", cread
                else: # otherwise we need to fix the alignment so that full read alignment is recorded
                    line.invalidateRead(cread)
                    print "Read intersects line center and alignments will be added later"
        for line in self.lines: # fixing read alignments for reads intersecting line center
            line.fixLineAlignments()
        for line in self.lines:
            print len(classified[line.id]), "reads were classified to line", line
            print map(str, classified[line.id])
        res = ReadCollection(reads.contigs).extend(list(itertools.chain(*classified.values())))
        print "Classified", len(res), "reads. Irrelevant:", n_irrelevant, ". Uncertain:", n_uncertain, ". Nonactive line:", n_class - len(res)
        return res

    def pairwiseAlign(self, lines):
        # type: (list[Line]) -> Dict[Tuple[int, int], list[AlignmentPiece]]
        print "Constructing pairwise line alignments"
        suffix_len = 10000
        for line in lines:
            suffix_len = min(len(line), suffix_len)
        assert suffix_len > 100
        shortened_lines = map(lambda line: Contig(line.seq[-suffix_len:], "short_" + str(line.id)), lines)
        lines_dict = dict()
        print "Aligning lines to lines"
        res = self.aligner.separateAlignments(shortened_lines, shortened_lines)
        res_map = dict() # type: Dict[Tuple[int, int], list[AlignmentPiece]]
        for l1 in lines:
            for l2 in lines:
                if l1 != l2:
                    res_map[(l1.id, l2.id)] = []
        for line in lines:
            lines_dict[line.id] = line
        print "Collecting line alignments"
        for read in res.reads.values():
            line_from = lines_dict[read.id[len("short_"):]]
            for al in read.alignments:
                print al
                line_to_id = read.id[len("short_"):]
                if line_to_id not in lines_dict:
                    continue
                line_to = lines_dict[line_to_id]
                if line_from == line_to:
                    continue
                seg_from = Segment(line_from, al.seg_from.left + len(line_from) - suffix_len, al.seg_from.right + len(line_from) - suffix_len)
                seg_to = Segment(line_to, al.seg_to.left + len(line_to) - suffix_len, al.seg_to.right + len(line_to) - suffix_len)
                # print seg_from, seg_to
                # print al.seg_from.Seq()
                # print seg_from.Seq()
                # print al.seg_to.Seq()
                # print al.cigar
                # print seg_to.Seq()
                new_piece = AlignmentPiece(seg_from, seg_to, al.cigar)
                print "Added line alignment:", new_piece
                res_map[(line_from.id, line_to.id)].append(new_piece)
        return res_map

    def fight(self, c1, c2, line_aligns):
        # type: (AlignmentPiece, AlignmentPiece, Dict[Tuple[int, int], list[AlignmentPiece]]) -> Optional[AlignmentPiece]
        assert c1.seg_from.contig == c2.seg_from.contig
        s1, s2, s12 = self.scorer.score3(c1, c2)
        if s12 is None:
            if s1 is None and s2 is not None:
                print "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "Winner:", c2
                return c2
            elif s1 is not None and s2 is None:
                print "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "Winner:", c1
                return c1
            elif s1 is None and s2 is None:
                print "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "No winner"
                return None
            assert False, "Strange comparison results"
        else:
            if s12 < 25 or (s12 < 100 and abs(s1 - s2) < s12 * 0.8) or abs(s1 - s2) < s12 * 0.65:
                print "Fight:", c1, c2, "Comparison results:", abs(s1 - s2), s12, s1, s2, "No winner"
                return None
            if s1 > s2:
                print "Fight:", c1, c2, "Comparison results:", abs(s1 - s2), s12, s1, s2, "Winner:", c2
                return c2
            else:
                print "Fight:", c1, c2, "Comparison results:", abs(s1 - s2), s12, s1, s2, "Winner:", c1
                return c1

    def tournament(self, candidates, line_aligns):
        # type: (list[AlignmentPiece], Dict[Tuple[int, int], list[AlignmentPiece]]) -> Optional[AlignmentPiece]
        best = None
        for candidate in candidates:
            if best is None:
                best = candidate
            else:
                best = self.fight(candidate, best, line_aligns)
        if best is None:
            return None
        if len(candidates) > 2:
            for candidate in candidates:
                if candidate == best:
                    continue
                fight_results = self.fight(candidate, best, line_aligns)
                if fight_results is None or fight_results != best:
                    return None
        return best

