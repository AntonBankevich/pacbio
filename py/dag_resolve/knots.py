import sys

from dag_resolve import params
from dag_resolve.line_tools import LineStorage, Line, Knot
from dag_resolve.repeat_graph import Edge
from dag_resolve.sequences import AlignedRead, Segment


class Knotter:
    def __init__(self, storage):
        # type: (LineStorage) -> Knotter
        self.storage = storage
        self.graph = storage.g

    def tryKnot(self, line1, line2):
        # type: (Line, Line) -> bool
        print "Trying to knot lines:", line1, line2
        same_vertex = line1.chain[-1].edge.end.rc == line2.chain[-1].edge.end
        extreme_case = line1.rc.id == line2.id and len(line1.chain) == 1 and len(line2.chain) == 1
        if line1.tail is None or same_vertex:
            seg1 = line1.chain[-1].edge.asSegment()
            reads1 = line1.chain[-1].reads
        else:
            seg1 = line1.tail.edgeSegment()
            reads1 = line1.tail.reads
            extreme_case = False
        if line2.tail is None or same_vertex:
            seg2 = line2.chain[-1].edge.asSegment()
            reads2 = line2.chain[-1].reads
        else:
            seg2 = line2.tail.edgeSegment()
            reads2 = line2.tail.reads
            extreme_case = False
        reads = reads1.cap(reads2)
        if extreme_case:
            print "Extreme case."
            reads =reads.filter(lambda read: self.extremeConnect(read, line1.chain[0].edge))
        else:
            reads = reads.filter(lambda read: self.checkConnect(read, seg1, seg2))
        print len(reads), "supporting reads"
        if len(reads) >= params.min_reads_in_knot:
            line1.knot = Knot(line1, line2, "", reads)
            line2.knot = Knot(line2, line1, "", reads)
            line1.tail = None
            line2.tail = None
            print "Knotted lines", line1, "and", line2
            return True
        else:
            return False

    def knotGraph(self):
        print "Knotting lines"
        for line1 in self.storage.lines:
            for line2 in self.storage.lines:
                if line1.knot is None and line2.knot is None and line1.id < line2.id:
                    if self.tryKnot(line1, line2):
                        break

    def checkConnect(self, read, seg1, seg2):
        # type: (AlignedRead, Segment, Segment) -> bool
        for aln1 in read.alignments:
            if not aln1.seg_to.inter(seg1):
                continue
            for aln2 in read.alignments:
                if not aln2.seg_to.inter(seg2):
                    continue
                if aln1.rc == aln2.rc:
                    continue
                if aln1.rc == (aln1.seg_from.right > aln2.seg_from.right):
                    return True
                if aln1.seg_to.contig == aln2.seg_to.contig and aln1.seg_from.left == aln2.seg_from.right and aln2.seg_from.right == aln2.seg_from.left:
                    return True
        return False

    def extremeConnect(self, read, edge):
        # type: (AlignedRead, Edge) -> bool
        for aln1 in read.alignments:
            if not aln1.seg_to.inter(edge.suffix(1000)):
                continue
            for aln2 in read.alignments:
                if not aln2.seg_to.inter(edge.prefix(1000)):
                    continue
                if aln1.seg_from.right <= aln2.seg_from.left and aln1.seg_to.right > aln2.seg_to.left:
                    return True
        return False




