import os
import sys

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print "flye_folder contigs"
    flye_dir = sys.argv[1]
    contigs = sys.argv[2]
    flye_contigs = os.path.join(flye_dir, "assembly.fasta")
    if "20-repeat" in os.listdir(flye_dir):
        edges = os.path.join(flye_dir, "20-repeat", "graph_before_rr.fasta")
        graph = os.path.join(flye_dir, "20-repeat", "graph_before_rr.gv")
    else:
        edges = os.path.join(flye_dir, "2-repeat", "graph_before_rr.fasta")
        graph = os.path.join(flye_dir, "2-repeat", "graph_before_rr.gv")
    eval_polish(contigs, flye_contigs, edges, graph)