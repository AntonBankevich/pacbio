import itertools
import os
import sys


sys.path.append("py")

from common.basic import CreateLog
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.smart_storage import SegmentStorage
from alignment.align_tools import Aligner, DirDistributor
from common import basic
from common.sequences import ContigStorage
from common.line_align import Scorer, events


def evaluatePI(dir, contigs_file, initial_file, ref_file):
    basic.ensure_dir_existance(dir)
    CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    contigs = ContigStorage().loadFromFasta(open(contigs_file, "r"), False)
    initial = ContigStorage().loadFromFasta(open(initial_file, "r"), False)
    ref = ContigStorage().loadFromFasta(open(ref_file, "r"), False)
    segs = []
    for al in aligner.overlapAlign(initial.unique(), contigs):
        if basic.isCanonocal(al.seg_to.contig.id):
            segs.append(al.seg_to)
        else:
            segs.append(al.rc.seg_to)
    segs = sorted(segs, key = lambda seg: basic.Normalize(seg.contig.id))
    interesting = dict()
    print "Interesting segments:"
    for contig, segit in itertools.groupby(segs, lambda seg: seg.contig):
        csegs = SegmentStorage().addAll(segit)
        csegs.mergeSegments()
        csegs = csegs.reverse(contig)
        interesting[contig.id] = list(csegs)
        print list(csegs)
    print "Analysis of contigs"
    scorer = Scorer()
    for al in aligner.localAlign(contigs.unique(), ref):
        print al
        for seg in interesting[al.seg_from.contig.id]:
            if al.seg_from.expand(500).contains(seg):
                tmp_al = al.reduce(query=al.seg_from.cap(seg))
                scorer.polyshMatching(tmp_al.matchingSequence())
                print tmp_al.seg_from, tmp_al.seg_to, str(events)
    print ""
    print "Analysis of initial"
    for al in aligner.overlapAlign(initial, ref):
        scorer.polyshMatching(al.matchingSequence())
        print al.seg_from, al.seg_to, str(events)




if __name__ == "__main__":
    dir = sys.argv[1]
    contigs = sys.argv[2]
    initial = sys.argv[3]
    ref = sys.argv[4]
    evaluatePI(dir, contigs, initial, ref)