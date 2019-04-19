import heapq
import sys
from typing import List, Tuple


sys.path.append("py")
from common import SeqIO
from common.sam_parser import Samfile
from common.sequences import ContigCollection, Segment
from common.alignment_storage import ReadCollection
from dag_resolve.repeat_graph import DotParser, Graph, Vertex


def ParseVertices(s):
    return s.split(";")

class PriorityQueue:
    def __init__(self):
        self.vals = [] #type: List[Tuple[int,Vertex]]

    def empty(self):
        # type: () -> bool
        return len(self.vals) == 0

    def push(self, item):
        # type: (Tuple[int, Vertex]) -> None
        heapq.heappush(self.vals, item)

    def pop(self):
        # type: () -> Tuple[int, Vertex]
        return heapq.heappop(self.vals)

dist = 5000

def main(argv):
    sys.stdout.write("Started\n")
    dot_file = argv[1]
    edge_sequences = argv[2]
    reference_file = argv[3]
    alignment_file = argv[4]
    edges = ParseVertices(argv[5])
    output_file = argv[6]
    sys.stdout.write("Loading dot\n")
    dot = DotParser(open(dot_file, "r")).parse()
    edge_collection = ContigCollection().loadFromFasta(open(edge_sequences, "r"), True)
    graph = Graph().loadFromDot(edge_collection, dot)
    vertices = [graph.E[id].start.id for id in edges]
    graph.printToFile(sys.stdout)
    print vertices
    ref = ContigCollection().loadFromFasta(open(reference_file, "r"), False)

    print "Looking for relevant"
    pq = PriorityQueue()
    for v in graph.V.values():
        if v.id in vertices:
            pq.push((0, v))
    visited = []
    while not pq.empty():
        d, v = pq.pop()
        if v in visited:
            continue
        visited.append(v)
        for e in v.inc:
            print e.id, e.start.id, e.end.id
            if d + len(e) < dist:
                pq.push((d + len(e), e.start))
        for e in v.out:
            print e.id, e.start.id, e.end.id
            if d + len(e) < dist:
                pq.push((d + len(e), e.end))
    print "Visited", len(visited)
    print map(str, list(visited))
    relevant = []
    edge_alignments = ReadCollection().loadFromFasta(open(edge_sequences, "r")).addAllRC()
    for edge in graph.E.values():
        if edge.start in visited or edge.start.rc in visited:
            relevant.append(edge_alignments[edge.id])
    print "Loading sam"
    edge_alignments.fillFromSam(Samfile(open(alignment_file, "r")), ref)
    for rel in relevant:
        print rel.__str__()
    print "Collecting segments"
    segments = []
    chr1 = ref["chr1"]
    for edge in relevant:
        for al in edge.alignments:
            print al
            if al.seg_from.inter(edge.prefix(dist)):
                l = dist - al.seg_from.left
                contig = al.seg_to.contig
                start = al.seg_to.left
                segments.append(Segment(contig, start, min(start + l, len(contig))))
                print segments[-1]
    tmp = []
    print "Rotating"
    for seg in segments:
        if seg.contig != chr1:
            seg = seg.RC()
        if seg.contig != chr1:
            print "WARNING", seg
        tmp.append(seg)
    segments = sorted(tmp, key = lambda seg: seg.left)
    print "All relevant segments"
    print "\n".join(map(str, segments))
    cur_seg = None
    interesting_segments = []
    print "Gluing"
    for seg in segments:
        if cur_seg is None:
            cur_seg = seg.copy()
            continue
        if cur_seg.right + 20000 < seg.left:
            interesting_segments.append(cur_seg.copy())
            cur_seg = seg.copy()
        else:
            cur_seg.right = max(cur_seg.right, seg.right)
    if cur_seg is not None:
        interesting_segments.append(cur_seg.copy())

    alignments = []
    for edge in edge_alignments:
        for al in edge.alignments:
            ok = False
            for seg in interesting_segments:
                if al.seg_to.inter(seg):
                    alignments.append(al)
    alignments = sorted(alignments, key = lambda al: al.seg_to.left)
    print "All relevant alignments"
    print "\n".join(map(str, alignments))

    print "Interesting segments:", len(interesting_segments), sum(map(len, interesting_segments))
    for seg in interesting_segments:
        print seg
    f = open(output_file, "w")
    tmp = []
    for seg in interesting_segments:
        SeqIO.write(SeqIO.SeqRecord(seg.Seq(), seg.__str__()), f, "fasta")
        tmp.append(seg.Seq())
    f.close()
    f1 = open(output_file + "1", "w")
    SeqIO.write(SeqIO.SeqRecord(("N" * 20000).join(tmp), "concat"), f1, "fasta")

if __name__ == "__main__":
    main(sys.argv)