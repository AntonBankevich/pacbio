import sys

from alignment.align_tools import Aligner
from common import basic
from dag_resolve import params
from dag_resolve.line_tools import LineStorage, Line, Knot
from dag_resolve.repeat_graph import Edge
from dag_resolve.sequences import AlignedRead, Segment, ReadCollection, ContigCollection, AlignmentPiece


class Knotter:
    def __init__(self, storage, aligner):
        # type: (LineStorage, Aligner) -> Knotter
        self.storage = storage
        self.graph = storage.g
        self.aligner = aligner

    def tryKnot(self, line1, line2):
        # type: (Line, Line) -> bool
        print "Trying to knot lines:", line1, line2
        # same_vertex = line1.chain[-1].seg_to.contig.end.rc == line2.chain[-1].seg_to.contig.end
        # extreme_case = line1.rc.id == line2.id and len(line1.chain) == 1 and len(line2.chain) == 1
        common_reads = ReadCollection(ContigCollection([line1, line2]))
        no_inter = 0
        aligned = 0
        unaligned = 0
        for read in line1.reads.cap(line2.reads):
            point = None
            for al1 in read.alignmentsTo(line1.centerPos.suffix()): # type: AlignmentPiece
                for al2 in read.alignmentsTo(line2.centerPos.prefix()): # type: AlignmentPiece
                    if al1 != al2:
                        common_reads.add(read)
                        point = (al1.seg_to.left - al1.seg_from.left, al2.seg_to.left - al2.seg_from.left)
            if point is None or aligned >= params.min_reads_in_knot:
                continue
            left1 = -point[0]
            right1 = len(line1) - point[0]
            left2 = -point[1]
            right2 = len(line2) - point[1]
            if left1 > left2 or right1 > right2:
                print "WARNING: STRANGE LINE KNOTTING"
            if right1 < left2 + 500:
                no_inter += 1
            else:
                seg1 = line1.suffix(left2 - right1)
                seg2 = line2.prefix(right1 - left2)
                segment_read = AlignedRead(seg1.subcontig())
                alignment = ReadCollection(ContigCollection([seg2.subcontig()]), [segment_read])
                self.aligner.alignReadCollection(alignment)
                found = False
                for al in segment_read.alignments:
                    if len(al.seg_from) > 0.8 * len(segment_read) and al.percentIdentity() > 0.8:
                        found = True
                if found:
                    aligned += 1
                else:
                    unaligned += 1

        print len(common_reads), "supporting reads.", "Aligned:", aligned, "Unaligned:", unaligned, "No inter:", no_inter
        for read in common_reads:
            print line1.reads[read.id], line1.reads[read.id] == line2.reads[read.id], self.graph.reads[read.id]
        if len(common_reads) >= params.min_reads_in_knot and (aligned >= params.min_reads_in_knot or no_inter >= params.min_reads_in_knot):
            line1.knot = Knot(line1, line2, "", common_reads)
            line2.rc.knot = Knot(line2.rc, line1.rc, "", common_reads.RC())
            print "Knotted lines", line1, "and", line2
            return True
        else:
            return False

    def knotGraph(self):
        for line in self.storage.lines:
            print "Line center:", line, len(line), line.centerPos.pos, line.rc.centerPos.pos, len(line.reads.inter(line.centerPos.suffix())), len(line.reads.inter(line.centerPos.prefix()))
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




