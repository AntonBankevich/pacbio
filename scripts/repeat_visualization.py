import itertools
import os
import sys


sys.path.append("py")
from common.basic import CreateLog
from typing import List
from common.alignment_storage import AlignmentPiece
from common.sequences import ContigStorage, Segment
from disjointig_resolve.dot_plot import DotPlot

from alignment.align_tools import Aligner, DirDistributor

class DisjointSet:
    def __init__(self):
        self.vals = dict()

    def add(self, item):
        self.vals[item] = item

    def get(self, item):
        if self.vals[item] != item:
            self.vals[item] = self.get(self.vals[item])
        return self.vals[item]

    def union(self, val1, val2):
        self.vals[self.get(val1)] = self.get(val2)

    def listComponenets(self):
        res = list(self.vals.keys())
        res = sorted(res, key = lambda item: self.get(item))
        return itertools.groupby(res, lambda item: self.get(item))

    def items(self):
        return list(self.vals.keys())

class Component:
    def __init__(self, segments, dot_plot):
        # type: (List[Segment], DotPlot) -> None
        self.segments = segments # type: List[Segment]
        self.alignments = [] # type: List[AlignmentPiece]
        for seg in self.segments:
            self.alignments.extend(filter(lambda al: not al.isIdentical(), dot_plot.allInter(seg)))

def ExtractRepeatComponents(contigs, aligner):
    # type: (ContigStorage, Aligner) -> List[Component]
    dot_plot = DotPlot(contigs)
    dot_plot.construct(aligner)
    segs = DisjointSet()
    for contig in contigs:
        for al in dot_plot.allInter(contig.asSegment()):
            if not al.isIdentical():
                segs.add(al.seg_to)
                segs.add(al.seg_from)
                segs.union(al.seg_from, al.seg_to)
    cur_seg = None #type: Segment
    for seg in sorted(segs.items(), key = lambda seg: (seg.contig.id, seg.left, seg.right)):
        if cur_seg is None or cur_seg.contig != seg.contig or cur_seg.dist(seg) > 1000:
            cur_seg = seg
        else:
            segs.union(cur_seg, seg)
            if seg.right > cur_seg.right:
                cur_seg = seg
    res = []
    for id, group in segs.listComponenets():
        comp_segs = []
        cur_seg = None
        for seg in sorted(group, key=lambda seg: (seg.contig.id, seg.left, seg.right)):
            if cur_seg is None or cur_seg.contig != seg.contig or cur_seg.dist(seg) > 500:
                if cur_seg is not None:
                    comp_segs.append(cur_seg)
                cur_seg = seg
            else:
                cur_seg = cur_seg.merge(seg)
        if cur_seg is not None:
            comp_segs.append(cur_seg)
        res.append(Component(comp_segs, dot_plot))
    return res


class Block:
    id_cnt = 0
    def __init__(self, segs):
        self.segs = segs
        self.x = None
        self.y = None
        self.out = [] # type: List[Block]
        self.inc = [] # type: List[Block]
        self.id = Block.id_cnt
        Block.id_cnt += 1

    def addOut(self, other):
        # type: (Block) -> None
        if other not in self.out:
            self.out.append(other)
            other.inc.append(self)




def CreateBlocks(comp):
    # type: (Component) -> List[Block]
    points = []
    for al in comp.alignments:
        if al.seg_from.left > 50:
            points.append(al.seg_to.prefix(length=1))
        if al.seg_to.left > 50:
            points.append(al.seg_from.prefix(length=1))
        if al.seg_from.right > 50:
            points.append(al.seg_to.suffix(length=1))
        if al.seg_to.right > 50:
            points.append(al.seg_from.suffix(length=1))
    new_points = []
    for al in comp.alignments:
        p = filter(lambda seg: al.seg_from.contains(seg), points)
        p = al.matchingSequence().mapPositionsDown([seg.left for seg in p], True)
        p = filter(lambda pos: pos is not None, p)
        new_points.extend([al.seg_to.contig.segment(pos, pos + 1) for pos in p])
    points.extend(new_points)
    seg_set = DisjointSet()
    seg_map = dict()
    for main_seg in comp.segments:
        seg_points = filter(lambda seg: main_seg.contains(seg), points)
        seg_points.append(main_seg.prefix(length=1))
        seg_points.append(main_seg.suffix(length=1))
        seg_points = sorted(seg_points, key=lambda seg: seg.left)
        # print "Points:", seg_points
        culsters = []
        prev = main_seg.prefix(length=1)
        for pos in seg_points[1:]:
            if pos.left > prev.right + 100:
                culsters.append((prev.right + prev.left) / 2)
                if len(prev) > 100:
                    print "Warning: large cluster.", prev
                prev = pos
            else:
                prev = prev.merge(pos)
        culsters.append((prev.right + prev.left) / 2)
        if len(prev) > 100:
            print "Warning: large cluster.", prev
        segments = []
        for p1, p2 in zip(culsters[:-1], culsters[1:]):
            segments.append(main_seg.contig.segment(p1, p2))
        segments[0].left = main_seg.left
        segments[-1].right = main_seg.right
        # print "Segments:", segments
        seg_map[main_seg] = segments
        for seg in segments:
            seg_set.add(seg)
    for al in comp.alignments:
        for seg_from in comp.segments:
            if not seg_from.inter(al.seg_from):
                continue
            for seg_to in comp.segments:
                if not seg_to.inter(al.seg_to):
                    continue
                for seg1 in seg_map[seg_from]:
                    if seg1.interSize(al.seg_from) > len(seg1) / 2:
                        al2 = al.reduce(query=seg1)
                        seg2 = al2.seg_to
                        for seg3 in seg_map[seg_to]:
                            if seg3.interSize(seg2) > len(seg3) * 3 / 4:
                                if seg_set.get(seg1) != seg_set.get(seg3):
                                    seg_set.union(seg1, seg3)
                                    # print "United blocks", seg1, seg3
                                    # print al2
    blocks = dict()
    res = []
    for key, iter in seg_set.listComponenets():
        block = Block(list(iter))
        res.append(block)
        for seg in block.segs:
            blocks[seg] = block
    for seg_list in seg_map.values():
        for seg1, seg2 in zip(seg_list[:-1], seg_list[1:]):
            blocks[seg1].addOut(blocks[seg2])
    return res



def draw(contigs_file, output_dir, k):
    aligner = Aligner(DirDistributor(os.path.join(output_dir, "alignments")))
    CreateLog(output_dir)
    print "Reading contigs"
    contigs = ContigStorage().loadFromFasta(open(contigs_file, "r"))
    print "Constructing components"
    componenets = ExtractRepeatComponents(contigs, aligner)
    print "Components:"
    for comp in componenets:
        print comp.segments
        print comp.alignments
    for cnt, comp in enumerate(componenets):
        print "Processing component", cnt
        print comp.segments
        # print comp.alignments
        print "Forming blocks"
        blocks = CreateBlocks(comp)
        for block in blocks:
            print "Block", block.id, ":", block.segs
        for block in blocks:
            for other in block.out:
                print block.id, "->", other.id

        print "Placing blocks on X axis"
        print "Placing blocks on Y axis"




if __name__ == "__main__":
    contig_file = sys.argv[1]
    output_dir = sys.argv[2]
    k = int(sys.argv[3])
    draw(contig_file, output_dir, k)
