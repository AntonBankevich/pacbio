import sys

from common import basic

sys.path.append("py")

from alignment.align_tools import Aligner, DirDistributor
from common.sequences import ContigCollection


def main(contigs_file, parts_file, dir):
    contigs = ContigCollection().loadFromFasta(open(contigs_file, "r"))
    parts = ContigCollection().loadFromFasta(open(parts_file, "r"))
    basic.CreateLog(dir)
    aligner = Aligner(DirDistributor(dir))
    res = dict()
    for al in aligner.localAlign(parts, contigs):
        if al.seg_to.contig.id not in res:
            res[al.seg_to.contig.id] = []
            res[al.seg_to.contig.rc.id] = []
        res[al.seg_to.contig.id].append(al)
        res[al.seg_to.contig.rc.id].append(al.rc)
    for cname, arr in res.items():
        print cname
        arr = filter(lambda al: len(al.seg_to) > min(len(al.seg_to.contig) - 1000, 5000), arr)
        arr = sorted(arr, key = lambda al: al.seg_to.left)
        print arr


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])