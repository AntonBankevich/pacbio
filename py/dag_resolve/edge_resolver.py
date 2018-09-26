from typing import Tuple, Optional

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common.SeqIO import NamedSequence
from dag_resolve import params
from dag_resolve.line_tools import Line
from dag_resolve.repeat_graph import Graph, Edge
from dag_resolve.sequences import ReadCollection, ContigCollection, AlignmentPiece, Segment


class EdgeResolver:
    def __init__(self, graph, aligner, polisher):
        # type: (Graph, Aligner, Polisher) -> EdgeResolver
        self.graph = graph
        self.aligner = aligner
        self.polisher = polisher

    def resolveEdge(self, edge, lines):
        # type: (Edge, list[Line]) -> Tuple[bool, Optional[Edge]]
        uncertain = edge.reads.minusAll([line.reads for line in lines])
        while True:
            extension_len = 0
            for line in lines:
                tmp = self.prolongConsensus(edge, line)
                extension_len += tmp
                print "Consensus of line", line.id, "was prolonged by", tmp
            if extension_len < 200:
                return False, None
            self.classifyReads(edge, lines, uncertain)

        return False, None

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
        return len(line) - old_len

    def classifyReads(self, edge, lines, reads):
        # type: (Edge, list[Line], ReadCollection) -> int
        pass
