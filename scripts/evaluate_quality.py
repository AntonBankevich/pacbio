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
    # type: (ContigCollection, List[AlignmentPiece]) -> Tuple[Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]], Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]], List[List[AlignmentPiece]]]
    als = als + map(lambda al: al.rc, als)
    als = sorted(als, key= lambda al: (len(al.seg_to.contig), al.seg_to.contig.id, al.seg_to.left, al.seg_to.right))
    res = dict()
    possible_term = dict()
    all_als = []
    for id, it in itertools.groupby(als, lambda al: al.seg_to.contig.id):
        contig_als = list(it)
        tmp = list(contig_als)
        starts = set([al.seg_from.contig.id for al in contig_als if al.seg_to.left < 50 and al.seg_from.left > 500])
        ends = [al for al in contig_als if al.rc.seg_to.left < 50 and al.rc.seg_from.left > 500]
        contig_als = [al for al in contig_als if not (al.rc.seg_to.left < 50 and al.rc.seg_from.left > 500)]
        for al in ends:
            if al.seg_from.contig.id not in starts:
                contig_als.append(al)
        all_als.append(contig_als)
        # print contig_als[0].seg_to.contig
        # print [al.seg_from.contig.id for al in contig_als]
        # print map(str, contig_als)
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

    return res, possible_term, all_als

def main(dir, contigs_files, reference_file, unique_contigs_file):
    CreateLog(dir)
    sys.stdout.level = LogPriority.warning
    ref = ContigCollection().loadFromFasta(open(reference_file, "r"), False)
    unique = ContigCollection().loadFromFasta(open(unique_contigs_file, "r"), False).filter(lambda contig: len(contig) > 5000)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    ref_als= list(aligner.overlapAlign(unique.unique(), ref))
    ref_transfers, ref_term, all_ref_als = extract_transfers(ref, ref_als)
    for uid in ref_term:
        ref_transfers[uid] = ref_term[uid]
    print "#", "file", "wrong", "unresolved", "correct", "missing"
    for i, contigs_file in enumerate(contigs_files):
        contigs = ContigCollection().loadFromFasta(open(contigs_file, "r"), False)
        contig_als = list(aligner.overlapAlign(unique.unique(), contigs))
        contig_transfers, contig_term, all_contig_als = extract_transfers(contigs, contig_als)
        missing = []
        wrong = dict()
        unresolved = []
        correct = []
        for uid in ref_transfers:
            if uid not in contig_transfers and uid not in contig_term:
                # print uid, "missing"
                missing.append(uid)
            elif uid in contig_transfers:
                if ref_transfers[uid][0] == contig_transfers[uid][0]:
                    # print uid, "correct"
                    correct.append(uid)
                else:
                    # print uid, "wrong", ref_transfers[uid][0].id, contig_transfers[uid][0].id
                    wrong[uid] = (ref_transfers[uid][0], contig_transfers[uid][0])
            else:
                if ref_transfers[uid][0] == contig_term[uid][0]:
                    # print uid, "correct"
                    correct.append(uid)
                else:
                    # print uid, "unresolved"
                    unresolved.append(uid)
        out = open(os.path.join(dir, "contigs_" + str(i) +".txt"), "w")
        out.write("Wrong: " + str(wrong) + "\n")
        out.write("Unresolved: " + str(unresolved) + "\n")
        out.write("Correct: " + str(correct) + "\n")
        out.write("Missing: " + str(missing) + "\n")
        out.write("Contig transfers: " + str(contig_transfers) + "\n")
        out.write("Contig term: " + str(contig_term) + "\n")
        out.write("Ref transfers: " + str(ref_transfers) + "\n")
        out.write("Ref als:\n")
        for c in all_ref_als:
            out.write(str(c) + "\n")
        out.write("Contig als:\n")
        for c in all_contig_als:
            out.write(str(c) + "\n")
        out.close()
        print "result", i, contigs_file, len(wrong), len(unresolved), len(correct), len(missing)

if __name__ == "__main__":
    dir = sys.argv[1]
    reference_file = sys.argv[2]
    unique_contigs_file = sys.argv[3]
    contigs_files = sys.argv[4:]
    main(dir, contigs_files, reference_file, unique_contigs_file)