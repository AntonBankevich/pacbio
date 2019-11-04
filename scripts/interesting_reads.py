import os
import sys

sys.path.append("py")
from alignment.align_tools import DirDistributor, Aligner
from common.basic import CreateLog
from common import params, basic, SeqIO
from common.scoring_model import ComplexScores
from common.sam_parser import Samfile
from common.sequences import ContigCollection, ContigStorage, Contig
from common.alignment_storage import ReadCollection, AlignedRead
from common.line_align import Scorer


def main(k, dir, contigs_file, reads_file):
    # type: (int, str, str, str) -> None
    basic.ensure_dir_existance(dir)
    CreateLog(dir)
    dd = DirDistributor(os.path.join(dir, "alignments"))
    aligner = Aligner(dd)
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
    res = open(os.path.join(dir, "reads.fasta"), "w")
    for read in reads:
        if not basic.isCanonocal(read.id):
            continue
        if len(read.alignments) > 1:
            SeqIO.write(read, res, "fasta")
    res.close()

if __name__ == "__main__":
    k = int(sys.argv[2])
    dir = sys.argv[3]
    contigs_file = sys.argv[4]
    reads_file = sys.argv[5]
    main(k, dir, contigs_file, reads_file)