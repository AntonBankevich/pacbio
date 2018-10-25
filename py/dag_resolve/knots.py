import sys

from common import basic
from dag_resolve import params
from dag_resolve.line_tools import LineStorage, Line, Knot
from dag_resolve.repeat_graph import Edge
from dag_resolve.sequences import AlignedRead, Segment, ReadCollection, ContigCollection


class Knotter:
    def __init__(self, storage):
        # type: (LineStorage) -> Knotter
        self.storage = storage
        self.graph = storage.g

    def tryKnot(self, line1, line2):
        # type: (Line, Line) -> bool
        print "Trying to knot lines:", line1, line2
        # same_vertex = line1.chain[-1].seg_to.contig.end.rc == line2.chain[-1].seg_to.contig.end
        # extreme_case = line1.rc.id == line2.id and len(line1.chain) == 1 and len(line2.chain) == 1
        common_reads = line1.reads.cap(line2.reads).filter(
            lambda read: len(list(read.alignmentsTo(line1.centerPos.suffix()))) != 0 and
                         len(list(read.alignmentsTo(line2.centerPos.prefix()))) != 0)
        print len(common_reads), "supporting reads"
        for read in common_reads:
            print line1.reads[read.id], line1.reads[read.id] == line2.reads[read.id], self.graph.reads[read.id]
        if len(common_reads) >= params.min_reads_in_knot:
            line1.knot = Knot(line1, line2, "", common_reads)
            line2.rc.knot = Knot(line2.rc, line1.rc, "", common_reads.RC())
            print "Knotted lines", line1, "and", line2
            return True
        else:
            return False

    def knotGraph(self):
        print "Knotting lines"
        for line1 in self.storage.lines:
            for line2 in self.storage.lines:
                if line1 != line2.rc:
                    self.tryKnot(line1, line2)
                # if line1.knot is None and line2.knot is None and line1 != line2.rc:
                #     if self.tryKnot(line1, line2):
                #         break

    def extremeConnect(self, read, edge):
        # type: (AlignedRead, Edge) -> bool
        for aln1 in read.alignments:
            if not aln1.seg_to.common(edge.suffix(-1000)):
                continue
            for aln2 in read.alignments:
                if not aln2.seg_to.common(edge.prefix(1000)):
                    continue
                if aln1.seg_from.right <= aln2.seg_from.left and aln1.seg_to.right > aln2.seg_to.left:
                    return True
        return False




