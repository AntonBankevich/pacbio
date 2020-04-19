import itertools
import os
import shutil
import sys
import time
sys.path.append("py")

from typing import List, Tuple, Dict

from disjointig_resolve import cl_params


from common.dot_parser import DotParser
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
    # type: (ContigCollection, List[AlignmentPiece]) -> Tuple[Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]], Dict[str, Tuple[Contig, AlignmentPiece, AlignmentPiece]]]
    als = als + map(lambda al: al.rc, als)
    als = sorted(als, key= lambda al: (al.seg_to.contig.id, al.seg_to.left, al.seg_to.right))
    res = dict()
    possible_term = dict()
    for id, it in itertools.groupby(als, lambda al: al.seg_to.contig.id):
        contig_als = list(it)
        if len(contig_als) == 0:
            continue
        print contig_als[0].seg_to.contig
        print [al.seg_from.contig.id for al in contig_als]
        print map(str, contig_als)
        for al1, al2 in zip(contig_als[:-1], contig_als[1:]):
            res[al1.seg_from.contig.id] = (al2.seg_from.contig, al1, al2)
            res[al2.rc.seg_from.contig.id] = (al1.rc.seg_from.contig, al2.rc, al1.rc)
        al1 = contig_als[-1]
        al2 = contig_als[0]
        if al1.seg_from.contig != al2.seg_from.contig:
            possible_term[al1.seg_from.contig.id] = (al2.seg_from.contig, al1, al2)
            possible_term[al2.rc.seg_from.contig.id] = (al1.rc.seg_from.contig, al2.rc, al1.rc)

    return res, possible_term

def main(dir, contigs_file, reference_file, unique_contigs_file):
    CreateLog(dir)
    contigs = ContigCollection().loadFromFasta(open(contigs_file, "r"), False)
    ref = ContigCollection().loadFromFasta(open(reference_file, "r"), False)
    unique = ContigCollection().loadFromFasta(open(unique_contigs_file, "r"), False).filter(lambda contig: len(contig) > 5000)
    aligner = Aligner(DirDistributor(os.path.join(dir, "alignments")))
    ref_als= list(aligner.overlapAlign(unique.unique(), ref))
    contig_als = list(aligner.overlapAlign(unique.unique(), contigs))
    ref_transfers, ref_term = extract_transfers(ref, ref_als)
    contig_transfers, contig_term = extract_transfers(contigs, contig_als)
    for uid in ref_term:
        ref_transfers[uid] = ref_term[uid]

    missing = 0
    wrong = 0
    unresolved = 0
    correct = 0
    for uid in ref_transfers:
        if uid not in contig_transfers and uid not in ref_term:
            print uid, "missing"
            missing += 1
        elif uid in contig_transfers:
            if ref_transfers[uid][0] == contig_transfers[uid][0]:
                print uid, "correct"
                correct += 1
            else:
                print uid, "wrong", ref_transfers[uid][0].id, contig_transfers[uid][0].id
                wrong += 1
        else:
            if ref_transfers[uid][0] == contig_term[uid][0]:
                print uid, "correct"
                correct += 1
            else:
                print uid, "unresolved"
                unresolved += 1
    print "Wrong:", wrong
    print "Unresolved:", unresolved
    print "Correct:", correct
    print "Missing:", missing

if __name__ == "__main__":
    dir = sys.argv[1]
    contigs_file = sys.argv[2]
    reference_file = sys.argv[3]
    unique_contigs_file = sys.argv[4]
    main(dir, contigs_file, reference_file, unique_contigs_file)