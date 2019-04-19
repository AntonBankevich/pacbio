import sys
sys.path.append("py")

from common.sam_parser import Samfile
from common.sequences import ContigCollection, loadFromSam
from common.alignment_storage import ReadCollection
from common.line_align import Scorer


# def main(contig_file, sam_file):
#     # type: (str, str) -> None
#     contigs = ContigCollection()
#     contigs.loadFromFasta(open(contig_file, "r"))
#     reads = loadFromSam(ReadCollection(), Samfile(open(sam_file, "r")), contigs)
#     scores = []
#     full_count = dict()
#     for read in reads:
#         full = 0
#         for al in read.alignments:
#             if len(al) > len(read) - 100:
#                 full += 1
#                 scores.append(float(Scorer().accurateScore(al.matchingSequence())) / len(read))
#         key = (full, len(read.alignments))
#         if key not in full_count:
#             full_count[key] = 0
#         full_count[key] += 1
#     full_count_list = []
#     for key in full_count:
#         full_count_list.append((full_count[key], key))
#     full_count_list = sorted(full_count_list)[::-1]
#     print "Full alignments:"
#     for num, key in full_count_list:
#         if num > 50:
#             print key, float(num) / len(reads)
#     print "Scores:"
#     scores = sorted(scores)
#     for score in scores[::1000]:
#         print score
#
# if __name__ == "__main__":
#     main(sys.argv[1], sys.argv[2])