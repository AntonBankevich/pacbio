import itertools
import os
import shutil
import sys
import time


from typing import List, Tuple, Dict

sys.path.append("py")

from common.log_params import LogPriority
from alignment.align_tools import Aligner, DirDistributor
from common import basic
from common.sequences import ContigCollection, Contig
from common.alignment_storage import AlignmentPiece


def CreateLog(dir):
    old_logs_dir = os.path.join(dir, "old")
    basic.ensure_dir_existance(old_logs_dir)
    log_file = os.path.join(dir, "log.info")
    if os.path.isfile(log_file):
        num = len(os.listdir(old_logs_dir))
        shutil.copy(log_file, os.path.join(old_logs_dir, str(num) + ".log"))
    log = open(log_file, "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stdout.prefix = lambda s: time.strftime("%I:%M:%S") + "  "
    sys.stderr = sys.stdout

def extract_transfers(contigs, als):
    # type: (ContigCollection, List[AlignmentPiece]) -> Tuple[Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]], Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]], Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]]]
    als = als + map(lambda al: al.rc, als)
    als = sorted(als, key= lambda al: (len(al.seg_to.contig), al.seg_to.contig.id, al.seg_to.left, al.seg_to.right))
    res = dict()
    possible_term = dict()
    for id, it in itertools.groupby(als, lambda al: al.seg_to.contig.id):
        contig_als = list(it)
        tmp = list(contig_als)
        starts = set([al.seg_from.contig.id for al in contig_als if al.seg_to.left < 50 and al.seg_from.left > 500])
        ends = [al for al in contig_als if al.rc.seg_to.left < 50 and al.rc.seg_from.left > 500]
        contig_als = [al for al in contig_als if not (al.rc.seg_to.left < 50 and al.rc.seg_from.left > 500)]
        for al in ends:
            if al.seg_from.contig.id not in starts:
                contig_als.append(al)
        # print contig_als[0].seg_to.contig
        # print [al.seg_from.contig.id for al in contig_als]
        # print map(str, contig_als)
        print id, ":", contig_als
        for al1, al2 in zip(contig_als[:-1], contig_als[1:]):
            if al1.seg_from.contig.id in possible_term:
                possible_term.pop(al1.seg_from.contig.id)
            if al2.rc.seg_from.contig.id in possible_term:
                possible_term.pop(al2.rc.seg_from.contig.id)
            res[al1.seg_from.contig.id] = (al2.seg_from.contig, al1, al2)
            res[al2.rc.seg_from.contig.id] = (al1.rc.seg_from.contig, al2.rc, al1.rc)
        if len(contig_als) == 0:
            continue
        al1 = contig_als[-1]
        al2 = contig_als[0]
        if len(contig_als) == 1 or al1.seg_from.contig != al2.seg_from.contig:
            if al1.seg_from.contig.id in res:
                res.pop(al1.seg_from.contig.id)
            if al2.rc.seg_from.contig.id in res:
                res.pop(al2.rc.seg_from.contig.id)
            possible_term[al1.seg_from.contig.id] = (al2.seg_from.contig, al1, al2)
            possible_term[al2.rc.seg_from.contig.id] = (al1.rc.seg_from.contig, al2.rc, al1.rc)
    all = dict(list(res.items()) + list(possible_term.items()))
    return res, possible_term, all

def main(dir, contigs_file1, contigs_file2, unique_contigs_file):
    CreateLog(dir)
    sys.stdout.level = LogPriority.warning
    unique = ContigCollection().loadFromFasta(open(unique_contigs_file, "r"), False).filter(lambda contig: len(contig) > 5000)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    contigs1 = ContigCollection().loadFromFasta(open(contigs_file1, "r"), False)
    cals1= list(aligner.overlapAlign(unique.unique(), contigs1))
    transfers1, term1, all1 = extract_transfers(contigs1, cals1)
    contigs2 = ContigCollection().loadFromFasta(open(contigs_file2, "r"), False)
    cals2 = list(aligner.overlapAlign(unique.unique(), contigs2))
    transfers2, term2, all2 = extract_transfers(contigs2, cals2)
    missing1 = []
    missing2 = []
    different = dict()
    unresolved1 = []
    unresolved2 = []
    same = []
    for ucontig in list(unique) + [contig.rc for contig in unique]:
        uid = ucontig.id
        in1 = uid in all1
        in2 = uid in all2
        if not in1 and not in2:
            continue
        if not in1:
            missing1.append(uid)
        elif not in2:
            missing2.append(uid)
        else:
            if all1[uid][0] == all2[uid][0]:
                same.append(uid)
            elif uid in transfers1 and uid in transfers2:
                different[uid] = (all1[uid][0], all2[uid][0])
            elif uid in transfers1:
                unresolved2.append(uid)
            elif uid in transfers2:
                unresolved1.append(uid)
    out = open(os.path.join(dir, "contigs.txt"), "w")
    out.write("Different: " + str(different) + "\n")
    out.write("Unresolved1: " + str(unresolved1) + "\n")
    out.write("Unresolved2: " + str(unresolved2) + "\n")
    out.write("Same: " + str(same) + "\n")
    out.write("Missing1: " + str(missing1) + "\n")
    out.write("Missing2: " + str(missing2) + "\n")
    out.write("Contig1 transfers: " + str(transfers1) + "\n")
    out.write("Contig1 term: " + str(term1) + "\n")
    out.write("Contig2 transfers: " + str(transfers2) + "\n")
    out.write("Contig2 term: " + str(term2) + "\n")
    out.close()
    print contigs_file1, contigs_file2
    print len(different), len(unresolved1), len(unresolved2), len(missing1), len(missing2), len(same)

if __name__ == "__main__":
    dir = sys.argv[1]
    unique_contigs_file = sys.argv[2]
    contigs_file1 = sys.argv[3]
    contigs_file2 = sys.argv[4]
    main(dir, contigs_file1, contigs_file2, unique_contigs_file)
