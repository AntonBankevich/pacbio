import itertools
import os
import shutil
import sys
import time

from typing import List, Tuple, Dict

sys.path.append("py")

from common.dot_parser import DotParser
from alignment.align_tools import Aligner, DirDistributor
from common import basic
from common.sequences import ContigCollection
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
    # type: (ContigCollection, List[AlignmentPiece]) -> Dict[Tuple[str, str], Tuple[AlignmentPiece, AlignmentPiece]]
    als = als + map(lambda al: al.rc, als)
    als = sorted(als, key= lambda al: (al.seg_to.contig.id, al.seg_to.left))
    res = dict()
    for id, it in itertools.groupby(als, lambda al: al.seg_to.contig.id):
        contig_als = list(it)
        print contig_als[0].seg_to.contig
        print [al.seg_from.contig.id for al in contig_als]
        print map(str, contig_als)
        for al1, al2 in zip(contig_als[:-1], contig_als[1:]):
            res[(al1.seg_from.contig.id, al2.seg_from.contig.id)] = (al1, al2)
            res[(al2.rc.seg_from.contig.id, al1.rc.seg_from.contig.id)] = (al2.rc, al1.rc)
    return res

def main(args):
    unique_file = sys.argv[1]
    graph_file = sys.argv[2]
    our_file = sys.argv[3]
    their_file = sys.argv[4]
    dir = sys.argv[5]
    min_cov = int(sys.argv[6])
    CreateLog(dir)
    unique_ids = [str(val[0]) for val in DotParser(open(graph_file, "r")).parse() if val[4].unique and val[4].cov >= min_cov]
    unique = ContigCollection()
    unique.loadFromFasta(open(unique_file, "r"), num_names=True)
    unique = unique.filter(lambda contig: contig.id in unique_ids)
    our = ContigCollection().loadFromFasta(open(our_file, "r"), False)
    their = ContigCollection().loadFromFasta(open(their_file, "r"), False)
    aligner = Aligner(DirDistributor(dir))
    our_als = list(aligner.overlapAlign(unique.unique(), our))
    their_als = list(aligner.overlapAlign(unique.unique(), their))
    our_transfers = extract_transfers(our, our_als)
    their_transfers = extract_transfers(their, their_als)

    out = open(os.path.join(dir, "res.txt"), "w")
    for o1, o2 in our_transfers:#type: str, str
        alo1, alo2 = our_transfers[(o1,o2)]# type: AlignmentPiece, AlignmentPiece
        if (o1,o2) in their_transfers:
            print "Common transition:", o1, o2
            alt1, alt2 = their_transfers[(o1,o2)] # type: AlignmentPiece, AlignmentPiece
            alot1 = alo1.composeTargetDifference(alt1)
            alot2 = alo2.composeTargetDifference(alt2)
            if alot1.seg_from.right < alot2.seg_from.left:
                copyo = alot1.seg_from.contig.segment(alot1.seg_from.right - 1, alot2.seg_from.left + 1).asContig()
                copyt = alot1.seg_to.contig.segment(alot1.seg_to.right - 1, alot2.seg_to.left + 1).asContig()
            else:
                copyo = alot1.seg_from.contig.segment(alot2.seg_from.left, alot1.seg_from.right).asContig()
                copyt = alot1.seg_to.contig.segment(alot2.seg_to.left, alot2.seg_to.right).asContig()
            diff = list(aligner.overlapAlign(ContigCollection([copyo]), ContigCollection([copyt])))
            if len(diff) != 1:
                print "Bad copy alignments"
                print map(str, diff)
            else:
                diff = diff[0]
                out.write(str(len(copyo)) + " " + str(len(copyt)) + " " + str(diff.percentIdentity()) + "\n")
                print [diff.seg_from.left, diff.seg_from.right], diff.seg_to, len(diff), diff.percentIdentity()
                print "\n".join(diff.asMatchingStrings())
        else:
            print "Our new transition:", o1, o2
    for o1, o2 in their_transfers:#type: str, str
        if (o1,o2) not in our_transfers:
            print "Missing transition:", o1, o2
    out.close()

if __name__ == "__main__":
    main(sys.argv)