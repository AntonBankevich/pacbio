import os
import sys


sys.path.append("py")
from alignment.align_tools import Aligner, DirDistributor
from common import basic
from common.sequences import ContigCollection, ReadCollection, Segment

def main(ref_file, segment, dir):
    ref = ContigCollection().loadFromFasta(open(ref_file, "r"), False)
    chr1 = ref["chr1"]
    if segment[0] < 0:
        segment = (-segment[0], -segment[1])
        chr1 = chr1.rc
    reads = ReadCollection(ref)
    reads_list = []
    for i in range(segment[0], segment[1], 500):
        read = reads.addNewRead(Segment(chr1, i, i + 500).asNamedSequence())
        reads_list.append(read)
    chr1.seq = chr1.seq[:segment[0]] + "N" * (segment[1] - segment[0]) + chr1.seq[segment[1]:]
    chr1.rc.seq = basic.RC(chr1.seq)
    basic.ensure_dir_existance(dir)
    aligner = Aligner(DirDistributor(dir))
    aligner.alignReadCollection(reads)
    out = sys.stdout
    for read in reads_list:
        # print read
        out.write(str(len(read.alignments)) + " " + str(max([0] + map(lambda al: al.percentIdentity(), read.alignments))) + "\n")
    out.close()

if __name__ == "__main__":
    main(sys.argv[1], (int(sys.argv[2]), int(sys.argv[3])), sys.argv[4])
