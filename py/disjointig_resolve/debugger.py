import os
import sys

from alignment.align_tools import Aligner
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.line_storage import NewLineStorage
from common.alignment_storage import ReadCollection
from common import basic


class Debugger:
    def __init__(self, dir, lines, dot_plot, reads, aligner):
        # type: (str, NewLineStorage, LineDotPlot, ReadCollection, Aligner) -> None
        self.dir = dir
        basic.ensure_dir_existance(dir)
        self.lines = lines
        self.dot_plot = dot_plot
        self.reads = reads
        self.aligner = aligner
        self.cnt = 0

    def dump(self):
        fname = os.path.join(self.dir, str(self.cnt) + ".txt")
        self.cnt += 1
        print "Dumping state to ", fname
        of = open(fname, "w")
        for line in self.lines:
            of.write(" ".join(map(str, [line.id, line, line.correct_segments, line.completely_resolved, line.initial])) + "\n")
        for line in self.lines:
            of.write(" ".join(map(str, [line.id, line, line.correct_segments, line.completely_resolved, line.initial])) + "\n")
            for al in line.read_alignments:
                of.write(str(al) + "\n")
        self.dot_plot.printAll(of)
        for al in self.aligner.localAlign(self.lines, self.lines):
            of.write(str(al) + "\n")
        for al in self.aligner.localAlign(self.reads, self.lines):
            of.write(str(al) + "\n")
        of.close()

debugger = None
