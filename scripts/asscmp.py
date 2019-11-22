import itertools
import os
import sys
sys.path.append("py")
from alignment.align_tools import Aligner, DirDistributor
from common.basic import CreateLog
from common.sequences import ContigStorage


def main(args):
    contigsf = sys.argv[1]
    reff = sys.argv[2]
    dir = sys.argv[3]
    CreateLog(dir)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    print "Reading contigs and reference"
    contigs = ContigStorage().loadFromFasta(open(contigsf, "r"), False)
    ref = ContigStorage().loadFromFasta(open(reff, "r"), False)
    print "Aligning reads"
    filtered = filter(lambda al: len(al) > 50000, aligner.referenceAlign(contigs, ref))
    print "Sorting reads along reference"
    als = sorted(filtered, key = lambda al: al.seg_to.contig.id)
    filtered = []
    for contig, it in itertools.groupby(als, lambda al: al.seg_to.contig):
        print contig
        cals = sorted(it, key = lambda al: -len(al))
        clen = 0
        long_filtered = []
        for al in cals:
            if clen  + 10000 > len(contig):
                return
            filtered.append(al)
            filtered.append(al.rc)
            long_filtered.append(al)
        res = []
        for read, it in itertools.groupby(sorted(long_filtered, key = lambda al: al.seg_from.contig.id), lambda al: al.seg_from.contig):
            als = sorted(it, key = lambda al : al.seg_to.left)
            cseg1 = als[0].seg_from
            cseg2 = als[0].seg_to
            for al1, al2 in zip(als[:-1], als[1:]):
                if al1.__le__(al2) and al1.seg_from.dist(al2.seg_from) < 30000 and al1.seg_to.dist(al2.seg_to) < 30000:
                    cseg1 = cseg1.merge(al2.seg_from)
                    cseg2 = cseg2.merge(al2.seg_to)
                else:
                    res.append((cseg1, cseg2))
                    cseg1 = al2.seg_from
                    cseg2 = al2.seg_to
        res_sorted = sorted(res, key = lambda al: al[1].left)
        print ", ".join("(%s[%d:%d]->%s[%d:%d])" % (f.contig.id, f.left, f.right, t.contig.id, t.left, t.right) for f, t in res_sorted)

    # print "Sorting reference along reads"
    # als = sorted(filtered, key = lambda al: al.seg_from.contig.id)
    # for contig, it in itertools.groupby(als, lambda al: al.seg_from.contig):
    #     print contig
    #     cals = sorted(it, key = lambda al: al.seg_from.left)
    #     print cals

if __name__ == "__main__":
    main(sys.argv)
