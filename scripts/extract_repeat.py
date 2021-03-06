import os
import sys

sys.path.append("py")

from common import SeqIO, basic
from common.basic import CreateLog
from common.seq_records import NamedSequence
from common.sequences import ContigStorage
from disjointig_resolve import cl_params


def parse(f):
    repeats = []
    starts = []
    ends = []
    for s in open(f, "r").readlines():
        s = s.split()
        if len(s) == 0:
            continue
        if s[0] == "Repeat":
            repeats.append(s[1])
        elif s[0] == "Start":
            starts.append(s[1:])
        elif s[0] == "End":
            ends.append(s[1:])
        else:
            assert False, "Incorrect dataset file"
    return repeats, starts, ends


def main(args):
    flye_dir = sys.argv[1]
    repeats, starts, ends = parse(sys.argv[2])
    graph_file, unique_file, disjointigs_file, rep_dir, tmp, their_file = cl_params.parseFlyeDir(flye_dir)
    dump = os.path.join(rep_dir, "read_alignment_dump")
    reads_file = sys.argv[3]
    dir = sys.argv[4]
    CreateLog(dir)
    print " ".join(args)
    print "Printing contigs"
    edges_file = os.path.join(rep_dir, "graph_before_rr.fasta")
    edges = ContigStorage().loadFromFasta(open(edges_file, "r"))
    unique = open(os.path.join(dir, "contigs"), "w")
    for l in starts:
        seq = "".join(map(lambda eid: edges[eid].seq, l))
        if len(seq) > 15000:
            seq = seq[-15000:]
        SeqIO.write(NamedSequence(seq, "(" + "_".join(l) + ")"), unique, "fasta")
    for l in ends:
        seq = "".join(map(lambda eid: edges[eid].seq, l))
        if len(seq) > 15000:
            seq = seq[:15000]
        SeqIO.write(NamedSequence(basic.RC(seq), "(" + "_".join(l) + ")"), unique, "fasta")
    unique.close()
    print "Selecting reads"
    reads = set()
    cur_read = None
    als = []
    for s in open(dump).readlines():
        if s.startswith("Chain"):
            if cur_read == "06dfd536-e258-446f-a019-8d8340ef831e":
                print als
            for al in als:
                if al in repeats:
                    if cur_read == "06dfd536-e258-446f-a019-8d8340ef831e":
                        print "oppa"
                    reads.add(cur_read)
                    break
            als = []
            
        else:
            s = s.split()
            cur_read = s[2][1:]
            eid = s[6].split("_")[1]
            if s[6][0] == "-":
                eid = "-" + eid
            if cur_read == "06dfd536-e258-446f-a019-8d8340ef831e":
                print eid
            als.append(eid)
    print "Selected", len(reads), "reads"
    print "\n".join(reads)
    print "Reading and printing reads"
    freads = open(os.path.join(dir, "reads.fasta"), "w")
    cnt = 0
    for read in SeqIO.parse_by_name(reads_file):
        cnt += 1
        if cnt % 10000 == 0:
            print cnt
        if read.id in reads:
            SeqIO.write(read, freads, "fasta")
    freads.close()








if __name__ == "__main__":
    main(sys.argv)
