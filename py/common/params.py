MINIMAP_BIN = "bin/flye-minimap2"

clean = False
reliable_coverage = 10
redo_alignments = False
radius = 8
alignment_correction_radius = 10
full_stats = False
num_iters = 2
min_reads_in_knot = 2
min_expected_Pacbio_PI = 0.78
min_allowed_Pacbio_PI = 0.7
max_jump = 6000
assert_pi = False
min_pi = 0.35
max_allowed_unaligned = 100
k = 1500
l = 2500
bad_end_length = 500
# k = 500
# l = 1500
threads = 16
min_k_mer_cov = 10
max_read_length = 80000

class Scores:
    ins_score = 8
    del_score = 7
    sub_score = 10
    homo_score = 4
    switch_core = 1
    center_score = 20
    inf = 10000000

class LogPriority:
    warning = 2
    log_level = 5
    alignment_files = 5
    main_parts = 0
    common = 1


downsample = 100000000


flanking_size = 500
window_size = 1500
