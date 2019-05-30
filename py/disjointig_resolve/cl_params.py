import getopt
import os
import sys

from typing import List

from common.save_load import TokenWriter, TokenReader


class Params:
    def __init__(self):
        self.reads_file = None
        self.contigs_file = None
        self.disjointigs_file_list = []
        self.disjointigs_file = None
        self.load_from = None
        self.graph_file = None
        self.flye_dir = None
        self.dir = None
        self.args = None
        self.threads = 8
        self.test = False
        self.long_params = "test stats flye-dir= graph= output-dir= reads= contigs= disjointigs= load= help".split(" ")
        self.short_params = "o:t:h"
        self.stats = False

    def check(self):
        if self.dir is None:
            print "Define output dir"
            self.print_usage_and_exit(1)
        if os.path.exists(self.dir) and os.path.isfile(self.dir):
            print "Incorrect output dir"
            self.print_usage_and_exit(1)

    def parse(self, argv):
        # type: (List[str]) -> Params
        self.args = argv
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], self.short_params, self.long_params)
            if len(tmp) != 0:
                self.print_usage_and_exit(1)
        except getopt.GetoptError:
            _, exc, _ = sys.exc_info()
            sys.stderr.write(str(exc) + "\n")
            self.print_usage_and_exit(1)
            return
        for (key, value) in options_list:
            if key == "--output-dir" or key == "-o":
                self.dir = value
                self.save_dir = os.path.join(self.dir, "saves")
                self.disjointigs_file = os.path.join(self.dir, "disjointigs.fasta")
            elif key == "--test":
                self.test = True
            elif key == "--flye-dir":
                self.flye_dir = value
                if self.graph_file is None:
                    self.graph_file = os.path.join(self.flye_dir, "assembly_graph.gv")
                if self.contigs_file is None:
                    self.contigs_file = os.path.join(self.flye_dir, "3-polishing", "polished_edges.fasta")
                self.disjointigs_file_list.append(os.path.join(self.flye_dir, "1-consensus", "consensus.fasta"))
                # self.disjointigs_file = os.path.join(self.flye_dir, "0-assembly", "draft_assembly.fasta")
            elif key == "--stats":
                self.stats = True
            elif key == "--graph":
                self.graph_file = value
            elif key == "--reads":
                self.reads_file = value
            elif key == "--contigs":
                self.contigs_file = value
            elif key == "--disjointigs":
                self.disjointigs_file_list.append(value)
            elif key == "--load":
                self.load_from = value
            elif key == "-t":
                self.threads = int(value)
            elif key == "--help" or key == "-h":
                self.print_usage_and_exit(0)
            else:
                self.print_usage_and_exit(1)
        return self

    def alignmentDir(self):
        return os.path.join(self.dir, "alignment")

    def saveDir(self):
        return os.path.join(self.dir, "saves")

    def print_usage_and_exit(self, code):
        # TODO: write usage
        print "Error in params"
        print self.long_params
        print self.short_params
        sys.exit(code)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.reads_file)
        handler.writeTokenLine(self.contigs_file)
        handler.writeTokenLine(self.disjointigs_file)
        handler.writeTokenLine(self.load_from)
        handler.writeTokenLine(self.dir)
        handler.writeTokenLine(self.save_dir)

    def load(self, handler):
        # type: (TokenReader) -> None
        self.reads_file = handler.readToken()
        self.contigs_file = handler.readToken()
        self.disjointigs_file = handler.readToken()
        self.load_from = handler.readToken()
        self.dir = handler.readToken()
        self.save_dir = handler.readToken()
