import getopt
import os
import sys

from typing import List

from common.save_load import TokenWriter, TokenReader


class Params:
    def __init__(self):
        self.reads_file = None
        self.contigs_file = None
        self.disjointigs_file = None
        self.load_from = None
        self.dir = None
        self.args = None
        self.threads = 8

    def parse(self, argv):
        # type: (List[str]) -> Params
        self.args = argv
        long_params = "output-dir= reads= contigs= disjointigs= load= help".split(" ")
        short_params = "o:t:h"
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], short_params, long_params)
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
            elif key == "--reads":
                self.reads_file = value
            elif key == "--contigs":
                self.contigs_file = value
            elif key == "--disjointigs":
                self.disjointigs_file = value
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
