import sys
sys.path.append("py")

from common import sam_parser
from common.sequences import ContigStorage
from common.alignment_storage import AlignmentPiece


def printAlignments(sam_handler, reference_handler, reads_handler):
    print "Loading reference"
    cc = ContigStorage(add_rc=False).loadFromFasta(reference_handler, False)
    print "Loading query"
    reads = ContigStorage().loadFromFasta(reads_handler, False)
    print "Loading result"
    res = []
    for rec in sam_parser.Samfile(sam_handler):
        if rec.query_name in reads.items and cc[rec.tname] is not None:
            al = AlignmentPiece.FromSamRecord(reads[rec.query_name], cc[rec.tname], rec)
            if al is None:
                print rec.query_name, rec.tname
                continue
            if al.seg_to.contig not in cc:
                al = al.rc
            res.append(al)
    print "Printing result", len(res)
    res = sorted(res, key = lambda al: al.seg_to.left)
#    res = sorted(res, key = lambda al: len(al))[::-1]
    up = 0
    down = 0
    for al in res:
        print al
        print list(al.splitRead())
        s1, s2 = al.asMatchingStrings()
        up += s1.count("-")
        down += s2.count("-")
        s = []
        if len(list(al.splitRead())) > 1:
            nums = []
            for al1 in al.splitRead():
                nums.append(al1.seg_from.left)
                nums.append(al1.seg_from.right - 1)
            cur_num = 0
            cur = al.seg_from.left

            for c in s1:
                if cur == nums[cur_num] and c != "-":
                    if cur_num % 2 == 0:
                        s.append("[")
                    else:
                        s.append("]")
                    cur_num += 1
                else:
                    if cur_num % 2 == 0:
                        s.append("-")
                    else:
                        s.append("+")
                if c != "-":
                    cur += 1
            print "".join(s)
        print s1
        print s2
    print up, down



if __name__ == "__main__":
    sam_handler = open(sys.argv[1], "r")
    ref_handler = open(sys.argv[2], "r")
    reads_handler = open(sys.argv[3], "r")
    printAlignments(sam_handler, ref_handler, reads_handler)
