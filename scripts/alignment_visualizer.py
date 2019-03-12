import sys
sys.path.append("py")

from common import sam_parser, basic
from common.seq_records import NamedSequence
from common.sequences import ContigCollection
from common.alignment_storage import AlignmentPiece, AlignedRead


def printAlignments(sam_handler, reference_handler, handler):
    cc = ContigCollection().loadFromFasta(reference_handler)
    # print cc.contigs.keys()
    for rec in sam_parser.Samfile(sam_handler):
        if rec.secondary or rec.is_unmapped or "H" in rec.cigar:
            continue
        if rec.rc:
            read = AlignedRead(NamedSequence(basic.RC(rec.seq), rec.query_name))
        else:
            read = AlignedRead(NamedSequence(rec.seq, rec.query_name))
        print read, rec.rc, rec.tname
        ap = read.AddSamAlignment(rec, cc[rec.tname])
        s1, s2 = ap.asMatchingStrings()
        handler.write(str(read) + "\n")
        handler.write(s1 + "\n" + s2 + "\n")


if __name__ == "__main__":
    sam_handler = open(sys.argv[1], "r")
    ref_handler = open(sys.argv[2], "r")
    out_handler = open(sys.argv[3], "w")
    printAlignments(sam_handler, ref_handler, out_handler)
    out_handler.close()
