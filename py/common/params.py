from common.scoring_model import SimpleScores

MINIMAP_BIN = "bin/flye-minimap2"

clean = False
reliable_coverage = 10
radius = 8
alignment_correction_radius = 30
alignment_smoothing_radius = 60
score_counting_radius = 10
full_stats = False
num_iters = 2
min_reads_in_knot = 2
min_expected_Pacbio_PI = 0.78
min_allowed_Pacbio_PI = 0.7
max_jump = 6000
assert_pi = False
min_pi = 0.35
max_allowed_unaligned = 200
k = 1500
l = 2500
bad_end_length = 500
# k = 500
# l = 1500
threads = 16
min_k_mer_cov = 10
min_contra_for_break = 5 # minimal nunber of unaligned reads needed to break resolved segment
max_read_length = 80000
min_alignment_size = 100 # minimal size of alignment that will be considered. Important for composition of alignments.
technology = "pacbio"
expected_size = 5000000


scores = SimpleScores(scoreIns = 8, scoreDel = 7, scoreMM = 10, scoreHomo = 4)
# class Scores:
#     ins_score = 8
#     del_score = 7
#     sub_score = 10
#     homo_score = 4
#     switch_core = 1
#     center_score = 20
#     inf = 10000000


downsample = 100000000


flanking_size = 500
window_size = 1500
