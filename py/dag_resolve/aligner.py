import os

from common import basic, sam_parser
from flye.alignment import make_alignment


class Aligner:
    def __init__(self, dir, threads = 16):
        self.dir = dir
        self.cur_alignment = 0
        self.threads = threads

    def next_dir(self):
        self.cur_alignment += 1
        name = os.path.join(self.dir, str(self.cur_alignment - 1))
        basic.ensure_dir_existance(name)
        return name

    def align(self, reads, consensus, threads = None):
        if threads is None:
            threads = self.threads
        dir = self.next_dir()
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
        make_alignment(contigs_file, [reads_file], threads, alignment_dir, "pacbio", alignment_file)
        return sam_parser.Samfile(alignment_file)