import sys

sys.path.append("py")

from common import basic, params
from common.basic import CreateLog

from alignment.align_tools import Aligner, DirDistributor
from common.line_align import Scorer
from common.sequences import ContigStorage

if __name__ == "__main__":
    basic.ensure_dir_existance(sys.argv[1])
    CreateLog(sys.argv[1])
    reads = ContigStorage().loadFromFile(sys.argv[2])
    contigs = ContigStorage().loadFromFile(sys.argv[3])
    scorer = Scorer()
    dd = DirDistributor(sys.argv[1])
    aligner = Aligner(dd)
    for read in reads.unique():
        print "Processing read", read
        als = [scorer.polyshAlignment(al, params.alignment_correction_radius) for al in aligner.localAlign([read], contigs)]
        for al1 in als:
            for al2 in als:
                if al1.seg_to.contig == al2.seg_to.contig:
                    continue
                print al1, "vs", al2
                scorer.scoreInCorrectSegments(al1, al1.seg_to.contig.asSegment(), al2, al2.seg_to.contig.asSegment())

