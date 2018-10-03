import sys
from typing import Tuple, Optional, Dict

from alignment.align_tools import Aligner, AlignedSequences
from alignment.polishing import Polisher
from common.SeqIO import NamedSequence
from dag_resolve import params
from dag_resolve.line_align import Scorer
from dag_resolve.line_tools import Line
from dag_resolve.repeat_graph import Graph, Edge
from dag_resolve.sequences import ReadCollection, ContigCollection, Segment, AlignedRead, Contig, AlignmentPiece, \
    UniqueList


class EdgeResolver:
    def __init__(self, graph, aligner, polisher):
        # type: (Graph, Aligner, Polisher) -> EdgeResolver
        self.graph = graph
        self.aligner = aligner
        self.polisher = polisher
        self.scorer = Scorer()
        self.prolonger = Prolonger(graph, aligner, polisher)
        self.classifier = ReadClassifier(graph, aligner)

    def resolveEdge(self, edge, lines):
        # type: (Edge, list[Line]) -> Tuple[bool, Optional[Edge]]
        print "Resolving edge", edge, "into lines:"
        for line in lines:
            print line.__str__()
        uncertain = edge.reads.minusAll([line.reads for line in lines])
        passedLines = []
        for line in UniqueList(lines):
            self.aligner.realignCollection(line.reads)
        # self.prolongAll(edge, lines)
        while True:
            shortest_lines, classified = self.classifier.classifyReads([self.shortestLine(lines, passedLines)], lines, uncertain)
            uncertain = uncertain.minus(classified)
            total_extention = 0
            jumped = 0
            if len(classified) != 0:
                for line in shortest_lines:
                    line_extention = self.prolonger.prolongConsensus(edge, line)
                    if self.attemptJump(edge, line):
                        passedLines.append(line)
                        jumped += 1
                    print "Consensus of line", line.id, "was prolonged by", line_extention
                    total_extention += line_extention
            if len(classified) == 0 and total_extention < 100 and jumped == 0:
                if self.prolongAll(edge, lines) < 300:
                    if len(passedLines) == len(lines):
                        return True, None
                    return False, None

    def attemptJump(self, edge, line):
        # type: (Edge, Line) -> bool
        shift = line.chain[-1].seg_from.right
        seq = line.seq[shift:]
        if len(seq) < 200:
            return False
        alignments = ReadCollection(ContigCollection(edge.end.out))
        alignments.addNewRead(NamedSequence(seq, "tail"))
        self.aligner.alignReadCollection(alignments)
        best = None
        for al in alignments.reads["tail"].alignments:
            if (al.seg_to.right > len(al.seg_to.contig) - 200 or al.seg_from.right > len(al.seg_from.contig) - 500) and \
                    (best is None or al.seg_from.left < best.seg_from.left):
                best = al
        if best is None:
            return False
        line.addAlignment(AlignmentPiece(Segment(line, best.seg_from.left + shift, best.seg_from.right + shift), best.seg_to, best.cigar))
        print "Connected line", line, "to edge", best.seg_to.contig
        return True

    def prolongAll(self, edge, lines):
        # type: (Edge, list[Line]) -> int
        res = 0
        for line in lines:
            tmp = self.prolonger.prolongConsensus(edge, line)
            res += tmp
            print "Line", line, "was prolonged by", tmp
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
        base_consensus = line.seq[-5000:]
        # print "inter"
        # line.reads.inter(line.suffix(-5000)).print_alignments(sys.stdout)
        reads = line.reads.inter(line.suffix(-5000)).noncontradicting(line.asSegment())
        print "noncontradicting"
        reads.print_alignments(sys.stdout)
        newConsensus = self.polisher.polishQuiver(reads, base_consensus, 4900).cut()
        print "New length:", len(newConsensus)
        if len(newConsensus) <= 120:
            return 0
        if newConsensus.seq[:100] != line.suffix(-100).Seq():
            print "Inaccurate glue"
            print newConsensus.seq[:100]
            print line.suffix(-100).Seq()
        old_len = len(line)
        line.extendRight(newConsensus.suffix(100))
        alignments = ReadCollection(ContigCollection([edge]))
        read = alignments.addNewRead(NamedSequence(newConsensus.seq[-100:], "tail"))
        self.aligner.alignReadCollection(alignments)
        read.sort()
        for al in read.alignments:
            if al.seg_to.contig == edge:
                seg_from = Segment(line, al.seg_from.left + 5000, al.seg_from.right + 5000)
                line.addAlignment(AlignmentPiece(seg_from, al.seg_to, al.cigar))
        if len(line) - old_len > 100:
            self.updateAlignments(line)
        return len(line) - old_len

    def updateAlignments(self, line):
        # type: (Line) -> None
        print "Updating alignment for line:", line
        line.reads = line.reads.cleanCopy(ContigCollection([line,line.rc]))
        line.rc.reads = line.reads
        self.aligner.alignReadCollection(line.reads)

class ReadClassifier:
    def __init__(self, graph, aligner):
        # type: (Graph, Aligner) -> ReadClassifier
        self.graph = graph
        self.aligner = aligner
        self.scorer = Scorer()

    def classifyReads(self, active, lines, reads):
        # type: (list[Line], list[Line], ReadCollection) -> Tuple[list[Line], ReadCollection]
        print "Trying to classify reads to lines", map(str, active)
        print "Full list of lines:", map(str, lines)
        alignments = ReadCollection(ContigCollection(lines))
        for read in reads.reads.values():
            alignments.addNewRead(read)
        for line in lines:
            self.aligner.alignReadCollection(alignments, [line])
        line_aligns = self.pairwiseAlign(lines)
        classified = dict()
        for line in active:
            classified[line.id] = []
        for read in alignments:
            candidates = []
            for al in read.alignments:
                if al.seg_to.contig in lines and \
                        (al.seg_to.right > len(al.seg_to.contig) - 500 or al.seg_from.right > len(al.seg_from.contig) - 500) and \
                                al.seg_from.left < 500 and len(al.seg_from) > 500:
                    candidates.append(al)
            res = self.championship(candidates, line_aligns)
            if res is None:
                print "Could not determine the champion"
                continue
            if res.seg_to.contig in active:
                classified[res.seg_to.contig.id].append(read)
                new_read = res.seg_to.contig.reads.addNewRead(read) #type: AlignedRead
                seg_from = Segment(new_read, res.seg_from.left, res.seg_from.right)
                new_read.addAlignment(AlignmentPiece(seg_from, res.seg_to))
        if len(classified) > 0:
            for line in active:
                self.aligner.expandCollection(line.reads, classified[line.id])
        res = ReadCollection(reads.contigs).extend(sum(classified.values()))
        print "Classified", len(res), "reads"
        return active, res

    def pairwiseAlign(self, lines):
        # type: (list[Line]) -> Dict[Tuple[int, int], list[AlignmentPiece]]
        for line in lines:
            assert len(line) > 5000
        shortened_lines = map(lambda line: Contig(line.seq[-5000:], line.id), lines)
        lines_dict = dict()
        res = self.aligner.separateAlignments(shortened_lines, shortened_lines)
        res_map = dict() # type: Dict[Tuple[int, int], list[AlignmentPiece]]
        for l1 in lines:
            for l2 in lines:
                if l1 != l2:
                    res_map[(l1.id, l2.id)] = []
        for line in lines:
            lines_dict[line.id] = line
        for read in res.reads.values():
            line_from = lines_dict[int(read.id)]
            for al in read.alignments:
                line_to = lines_dict[int(al.seg_to.contig.id)]
                if line_from == line_to:
                    continue
                seg_from = Segment(line_from, al.seg_from.left + len(line_from) - 5000, al.seg_from.right + len(line_from) - 5000)
                seg_to = Segment(line_to, al.seg_to.left + len(line_to) - 5000, al.seg_to.right + len(line_to) - 5000)
                res_map[(int(read.id), al.seg_to.contig.id)].append(AlignmentPiece(seg_from, seg_to, al.cigar))
        return res_map

    def fight(self, c1, c2, line_aligns):
        # type: (AlignmentPiece, AlignmentPiece, Dict[Tuple[int, int], list[AlignmentPiece]]) -> Optional[AlignmentPiece]
        assert c1.seg_from.contig == c2.seg_from.contig
        s1, s2, s12 = self.scorer.score3(c1, c2, line_aligns[(c1.seg_to.contig.id, c2.seg_to.contig.id)])
        if s12 is None:
            if s1 is None and s2 is not None:
                return c2
            elif s1 is not None and s2 is None:
                return c1
            assert False, "Strange comparison results"
        else:
            print "Comparison results:", abs(s1 - s2), s12, s1, s2
            if s12 < 30 or abs(s1 - s2) < s12 * 0.8:
                print "No winner"
                return None
            if s1 > s2:
                return c1
            else:
                return c2

    def championship(self, candidates, line_aligns):
        # type: (list[AlignmentPiece], Dict[Tuple[int, int], list[AlignmentPiece]]) -> Optional[AlignmentPiece]
        best = None
        for candidate in candidates:
            if best is None:
                best = candidate
            else:
                best = self.fight(candidate, best, line_aligns)
        if best is None:
            return None
        for candidate in candidates:
            if candidate == best:
                continue
            if self.fight(candidate, best, line_aligns) != best:
                return None
        return best

