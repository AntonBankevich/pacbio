import json
import sys
import os

import pickle

sys.path.append("py")
from common import sam_parser, basic
from common.basic import RC
from flye.alignment import make_alignment
from flye import polysh_job
import sequences

import common.SeqIO as SeqIO

def MakeAlignment(reads, consensus, dir):
    contigs_file = os.path.join(dir, "contigs.fasta")
    reads_file = os.path.join(dir, "reads.fasta")
    alignment_dir = os.path.join(dir, "alignment")
    alignment_file = os.path.join(dir, "alignment.sam")
    basic.ensure_dir_existance(dir)
    basic.ensure_dir_existance(alignment_dir)
    reads_handler = open(reads_file, "w")
    reads.print_fasta(reads_handler)
    reads_handler.close()
    contigs_handler = open(contigs_file, "w")
    consensus.print_fasta(contigs_handler)
    contigs_handler.close()
    make_alignment(contigs_file, [reads_file], 8, alignment_dir, "pacbio", alignment_file)
    return alignment_file

def CountCoverage(contig, alignment_file):
    res = [0] * (len(contig) + 1)
    for rec in sam_parser.Samfile(open(alignment_file, "r")):
        res[rec.pos] += 1
        res[min(rec.pos + rec.alen, len(contig))] -= 1
    res_list = [[0, 0, res[0]]]
    cur = res[0]
    for i in range(1, len(res)):
        cur += res[i]
        if res[i] == 0:
            res_list[-1][1] = i
        else:
            res_list.append([i,i,cur])
    return res_list


def Polish(reads, consensus, dir):
    basic.ensure_dir_existance(dir)
    consensus_file_name = os.path.join(dir, "ref.fasta")
    consensus_file = open(consensus_file_name, "w")
    SeqIO.write(consensus, consensus_file, "fasta")
    consensus_file.close()
    reads_file_name = os.path.join(dir, "reads.fasta")
    reads_file = open(reads_file_name, "w")
    reads.print_fasta(reads_file)
    reads_file.close()
    res = polysh_job.polish_from_disk(dir, consensus_file_name, reads_file_name)
    return list(SeqIO.parse_fasta(open(res, "r")))[0]



# 1 - contigs names 2 - contigs file 3- reads file 4 - dir
def construct_alignments(contigs_names, initial_contigs_file, reads_file, dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

    sys.stderr.write("Collecting contigs\n")
    contig_names = dict()
    for rec in contigs_names.split(";"):
        tmp = rec.split(",")
        contig_names[tmp[0]] = tmp[1:]
    contigs = sequences.ContigCollection()
    contigs_file = os.path.join(dir, "contigs.fasta")
    contigs_handler = open(contigs_file, "w")
    for rec in SeqIO.parse_fasta(open(initial_contigs_file, "r")):
        if rec.id in contig_names:
            mods = contig_names[rec.id]
            if "RC" in mods:
                rec.seq = RC(rec.seq)
            # rec.id += ";" + ";".join(mods)
            SeqIO.write(rec, contigs_handler, "fasta")
            contigs.add(sequences.Contig(rec.seq, rec.id, mods))
    contigs_handler.close()
    contigs.print_names(sys.stderr)
    alignment_dir = os.path.join(dir, "alignment")
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)

    sys.stderr.write("Aligning\n")
    alignment = os.path.join(dir, "alignment.sam")
    make_alignment(contigs_file, [reads_file], 8, alignment_dir, "pacbio", alignment)
    sys.stderr.write("Filtering out irrelevant\n")
    aligned = set()
    for rec in sam_parser.Samfile(open(alignment, "r")):
        if not rec.is_unmapped and "main" in contigs[rec.tname].info.misc:
            aligned.add(rec.query_name)

    sys.stderr.write("Collecting reads\n")

    output = open(os.path.join(dir, "filtered_reads.fasta"), "w")
    reads = sequences.ReadCollection(contigs)
    for rec in SeqIO.parse_fasta(open(reads_file, "r")):
        rname = rec.id.split()[0]
        rec.id = rec.id.split()[0]
        if rname in aligned:
            SeqIO.write(rec, output, "fasta")
            reads.addNewRead(rec)
    output.close()

    sys.stderr.write("Collecting alignments\n")
    for rec in sam_parser.Samfile(open(alignment, "r")):
        reads.addNewAlignment(rec)

    dump = open(os.path.join(dir, "aln.txt"), "w")
    reads.print_alignments(dump)
    dump.close()

    polished = open(os.path.join(dir, "polished.fasta"), "w")
    cnt = 0
    for contig in contigs.incoming():
        sys.stderr.write("Polishing incoming contig " + contig.id + " of length " + str(len(contig)) + "\n")
        sorted_reads = reads.inter(contig.asSegment())
        dump = open(os.path.join(dir, "aln" + contig.id+ ".txt"), "w")
        sorted_reads.print_alignments(dump)
        dump.close()
        rec = Polish(sorted_reads, contigs.main, os.path.join(dir, "tmp"))
        rec.id += "_polished_copy_" + contig.id
        SeqIO.write(rec, polished, "fasta")
        polished_contig = sequences.Contig(rec.seq, rec.id, [])
        cov_dir = os.path.join(dir, "coverage_" + polished_contig.id)
        print cov_dir
        coverage_list = CountCoverage(polished_contig, MakeAlignment(sorted_reads, polished_contig, cov_dir))
        cov_file = open(os.path.join(cov_dir, "cov.txt"), "w")
        cov_file.write("\n".join(map(str, coverage_list)))
        cov_file.close()
        cnt += 1

    polished.close()

# 1 - contigs names 2 - contigs file 3- reads file 4 - dir
if __name__ == "__main__":
    construct_alignments(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
