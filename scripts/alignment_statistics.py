import os
import sys

sys.path.append("py")
from alignment.align_tools import DirDistributor, Aligner
from common.basic import CreateLog
from common import params, basic
from common.scoring_model import ComplexScores
from common.sam_parser import Samfile
from common.sequences import ContigCollection, ContigStorage, Contig
from common.alignment_storage import ReadCollection, AlignedRead
from common.line_align import Scorer


def main(model_file, k, dir, contigs_file, reads_file):
    # type: (str, int, str, str, str) -> None
    basic.ensure_dir_existance(dir)
    CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
    params.scores = ComplexScores()
    params.scores.load(open(model, "r"))
    params.k = k
    print "Loading contigs"
    tmp = sorted(ContigStorage().loadFromFasta(open(contigs_file, "r"), False).unique(), key = lambda contig: len(contig))
    cnt = 1
    contigs = ContigStorage()
    for c1, c2 in zip(tmp[::2], tmp[1::2]):
        # if c1.seq == c2.rc.seq:
        contigs.add(Contig(c1.seq, str(cnt)))
        print cnt, c1.id, c2.id
        cnt += 1
        # else:
        #     contigs.add(Contig(c1.seq, str(cnt)))
        #     print cnt, c1.id
        #     cnt += 1
        #     contigs.add(Contig(c2.seq, str(cnt)))
        #     print cnt, c2.id
        #     cnt += 1
    print "Loading reads"
    reads = ReadCollection().loadFromFasta(open(reads_file, "r"))
    print "Aligning reads"
    for al in aligner.localAlign(reads, contigs):
        if len(al) > k:
            read = al.seg_from.contig # type:AlignedRead
            read.addAlignment(al)
    for read in reads:
        if not basic.isCanonocal(read.id):
            continue
        cnt = 0
        al0 = None
        others = []
        for al in read.alignments:
            if al.contradictingRTC():
                cnt += 1
                al0 = al
            else:
                others.append(al)
        if cnt != 1 or len(others) == 0:
            continue
        print al0
        print others
        seg = al0.seg_from
        for al in others:
            if al.seg_from.interSize(seg) < k:
                seg = None
                break
            else:
                seg = al.seg_from.cap(seg)
        print seg
        if seg == None:
            continue
        al0 = al0.reduce(query=seg)
        others = [al.reduce(query=seg) for al in others]
        scorer = Scorer(params.scores)
        print "Winner"
        for al in others:
            scorer.scoreCommon(al0, al)
        print "No winner"
        for al1 in others:
            for al2 in others:
                if al1 == al2:
                    continue
                scorer.scoreCommon(al1, al2)

if __name__ == "__main__":
    model = sys.argv[1]
    k = int(sys.argv[2])
    dir = sys.argv[3]
    contigs_file = sys.argv[4]
    reads_file = sys.argv[5]
    main(model, k, dir, contigs_file, reads_file)