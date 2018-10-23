clean = False
reliable_coverage = 10
radius = 8
alignment_correction_radius = 4
full_stats = False
num_iters = 2
min_reads_in_knot = 3
min_expected_Pacbio_PI = 0.8

class Scores:
    ins_score = 8
    del_score = 7
    sub_score = 10
    homo_score = 4
    switch_core = 1
    center_score = 20
    inf = 10000000