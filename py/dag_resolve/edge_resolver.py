from typing import Tuple, Optional

from alignment.align_tools import Aligner, AlignedSequences
from alignment.polishing import Polisher
from common.SeqIO import NamedSequence
from dag_resolve import params
from dag_resolve.line_tools import Line
from dag_resolve.repeat_graph import Graph, Edge
from dag_resolve.sequences import ReadCollection, ContigCollection, AlignmentPiece, Segment, AlignedRead


class EdgeResolver:
    def __init__(self, graph, aligner, polisher):
        # type: (Graph, Aligner, Polisher) -> EdgeResolver
        self.graph = graph
        self.aligner = aligner
        self.polisher = polisher

    def resolveEdge(self, edge, lines):
        # type: (Edge, list[Line]) -> Tuple[bool, Optional[Edge]]
        print "Resolving edge", edge, "into lines:"
        for line in lines:
            print line.__str__()
        uncertain = edge.reads.minusAll([line.reads for line in lines])
        passedLines = []
        while True:
            line, classified = self.classifyReads(self.shortestLine(lines, passedLines), lines, uncertain)
            print "Classified", len(classified), "reads to line", line
            if len(classified) == 0:
                return False, None
            extension_len = self.prolongConsensus(edge, line)
            self.updateAlignments(line)
            print "Consensus of line", line.id, "was prolonged by", extension_len
            if extension_len < 100:
                if self.prolongAll(edge, lines) < 300:
                    if len(passedLines) == len(lines):
                        return True, None
                    if self.attemptJump(edge, self.shortestLine(lines, passedLines)):
                        continue
                    else:
                        return False, None

    def attemptJump(self, edge, line):
        # type: (Edge, Line) -> bool
        shift = line.chain[-1].seg_from.right
        seq = line.seq[shift:]
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

    def updateAlignments(self, line):
        # type: (Line) -> None
        print "Updating alignment for line:", line
        line.reads = line.reads.cleanCopy(ContigCollection([line,line.rc]))
        line.rc.reads = line.reads
        self.aligner.alignReadCollection(line.reads)

    def prolongAll(self, edge, lines):
        # type: (Edge, list[Line]) -> int
        res = 0
        for line in lines:
            res += self.prolongConsensus(edge, line)
        return res

    def prolongConsensus(self, edge, line):
        # type: (Edge, Line) -> int
        base_consensus = line.seq[-5000:]
        reads = line.reads.inter(line.suffix(-5000))
        newConsensus = self.polisher.polishQuiver(reads, base_consensus, 4900).cut()
        if len(newConsensus) <= 100:
            return 0
        if newConsensus.seq[:100] != line.suffix(-100).Seq():
            print "Incaccurate glue"
            print newConsensus.seq[:100]
            print line.suffix(100).Seq()
        old_len = len(line)
        line.extendRight(newConsensus.suffix(100))
        alignments = ReadCollection(ContigCollection([edge]))
        read = alignments.addNewRead(NamedSequence(newConsensus.suffix(100), "tail"))
        self.aligner.alignReadCollection(alignments)
        read.sort()
        for al in read.alignments:
            if al.seg_to.contig == edge:
                seg_from = Segment(line, al.seg_from.left + 5000, al.seg_from.right + 5000)
                line.addAlignment(AlignmentPiece(seg_from, al.seg_to, al.cigar))
        self.updateAlignments(line)
        return len(line) - old_len

    def classifyReads(self, shortest, lines, reads):
        # type: (Line, list[Line], ReadCollection) -> Tuple[line, list[AlignedRead]]
        print "Trying to classify reads to line", shortest
        alignments = ReadCollection(ContigCollection(lines))
        for read in reads.reads.values():
            alignments.addNewRead(read)
        for line in lines:
            self.aligner.alignReadCollection(alignments, [line])

        classified = []
        for read in alignments:
            candidates = []
            for al in read.alignments:
                if al.seg_to.contig in lines and \
                        (al.seg_to.right > len(al.seg_to.contig) - 500 or al.seg_from.right > len(al.seg_from.contig) - 500) and \
                                al.seg_from.left < 500:
                    candidates.append(al)
            res = self.championship(candidates)
            if res is None:
                print "Could not determine the champion"
                continue
            if res == shortest:
                classified.append(read)
                shortest.reads.addNewRead(read)
        if len(classified) > 0:
            self.updateAlignments(shortest)
        return shortest, classified

    def fight(self, c1, c2):
        # type: (AlignmentPiece, AlignmentPiece) -> Optional[AlignmentPiece]
        assert c1.seg_from.contig == c2.seg_from.contig

        return None

    def championship(self, candidates):
        # type: (list[AlignmentPiece]) -> Optional[AlignmentPiece]
        best = None
        for candidate in candidates:
            if best is None:
                best = candidate
            else:
                best = self.fight(candidate, best)
        if best is None:
            return None
        for candidate in candidates:
            if candidate == best:
                continue
            if self.fight(candidate, best) != best:
                return None
        return best

    def shortestLine(self, lines, passedLines):
        shortest = None
        for line in lines:
            if line not in passedLines and (shortest is None or line.chain[-1].seg_to.right < shortest.chain[-1].seg_to.right):
                shortest = line
        return shortest
