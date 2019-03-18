import itertools
import sys
import random


sys.path.append("py")
from common.line_align import Scorer
from alignment.align_tools import DirDistributor, Aligner
from common import basic
from common.sequences import ContigCollection, Contig, ReadCollection
import common.SeqIO as SeqIO


def main(reads_file, ref_file, dir, error_rate):
    sys.stderr.write("Reading reference" + "\n")
    ref = sorted(list(SeqIO.parse_fasta(open(ref_file, "r"))), key = lambda rec: len(rec))[-1]
    ref = Contig(ref.seq, ref.id)
    refs = ContigCollection()
    for i in range(0, len(ref) - 500, 500):
        if random.random() > 0.95:
            tmp = list(ref.segment(i, i + 500).Seq())
            for j in range(error_rate * 500 / 100):
                pos = random.randint(0, 499)
                tmp[pos] = basic.rc[tmp[pos]]
            refs.add(Contig("".join(tmp), ref.id + "(" + str(i) + "," + str(i + 500) + ")"))
    refs.print_names(sys.stderr)
    sys.stderr.write("Reading reads" + "\n")
    reads = ReadCollection()
    reads.loadFromFasta(open(reads_file, "r"))

    sys.stderr.write("Aligning reads" + "\n")
    basic.ensure_dir_existance(dir)
    aligner = Aligner(DirDistributor(dir))
    aligner.alignReadCollection(reads, refs)
    sys.stderr.write("Analysing alignments" + "\n")
    alignments = []
    for read in reads:
        alignments.extend(read.alignments)
    alignments = filter(lambda al: len(al) > 450, alignments)
    alignments = sorted(alignments, key = lambda al: (al.seg_to.contig.id, al.seg_from.contig.id))
    scorer = Scorer()
    scorer.scores.homo_score = 3
    scorer.scores.ins_score = 5
    scorer.scores.del_score = 5
    cnt = 0
    for contig, iter in itertools.groupby(alignments, key=lambda al:al.seg_to.contig):
        iter = list(iter)
        sys.stderr.write(str(contig) + " " + str(len(iter)) + "\n")
        if len(iter) < 150:
            for al in iter:
                print scorer.accurateScore(al.matchingSequence())
                cnt += 1
                if cnt >= 5000:
                    break
        if cnt >= 5000:
            break




if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))