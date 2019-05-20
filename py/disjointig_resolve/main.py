import itertools
import os
import shutil
import sys
import time

sys.path.append("py")

from disjointig_resolve.unique_marker import UniqueMarker
from alignment.align_tools import Aligner, DirDistributor
from alignment.polishing import Polisher
from common import basic
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.tests import Tester
from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.sequences import ContigCollection
from common.alignment_storage import ReadCollection
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.saves_io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.cl_params import Params


def CreateLog(params):
    old_logs_dir = os.path.join(params.dir, "old")
    basic.ensure_dir_existance(old_logs_dir)
    log_file = os.path.join(params.dir, "log.info")
    if os.path.isfile(log_file):
        num = len(os.listdir(old_logs_dir))
        shutil.copy(log_file, os.path.join(old_logs_dir, str(num) + ".log"))
    log = open(log_file, "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stderr = sys.stdout
    print " ".join(params.args)


def main(args):
    params = Params().parse(args)
    params.check()
    CreateLog(params)
    sys.stdout.write("Started\n")
    print (time.strftime("%d.%m.%Y  %I:%M:%S"))
    if params.test:
        aligner = Aligner(DirDistributor(params.alignmentDir()))
        Tester(aligner).testAll("tests/cases.txt")
        sys.stdout.write("Finished\n")
        print (time.strftime("%d.%m.%Y  %I:%M:%S"))
        return
    print "Preparing initial state"
    if params.load_from is not None:
        print "Loading initial state from saves"
        params, aligner, contigs, reads, disjointigs, lines, dot_plot = loadAll(TokenReader(open(params.load_from, "r")))
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
    else:
        aligner = Aligner(DirDistributor(params.alignmentDir()))

        print "Creating disjointig collection"
        disjointigs = DisjointigCollection()
        disjointigs.loadFromFasta(open(params.disjointigs_file, "r"))

        print "Creating read collection"
        reads = ReadCollection()
        reads.loadFromFasta(open(params.reads_file, "r"))

        print "Aligning reads to disjointigs"
        disjointigs.addAlignments(aligner.alignClean(reads, disjointigs))

        print "Creating contig collection"
        contigs = ContigCollection()
        contigs.loadFromFasta(open(params.contigs_file, "r"), num_names=True)

        print "Creating line collection"
        lines = NewLineStorage(disjointigs, aligner)
        lines.fillFromContigs(contigs)
        lines.alignDisjointigs()
        # lines.fillFromDisjointigs()

        dot_plot = LineDotPlot(lines, aligner)
        dot_plot.construct(aligner)

        UniqueMarker().markAllUnique(lines, dot_plot)
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
        extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))

    print "Resolving"
    save_handler = SaveHandler(params.save_dir)
    cnt = 0
    while True:
        stop = True
        for line_id in list(lines.items.keys()):
            if line_id not in lines.items:
                continue
            line = lines[line_id]
            extended = extender.tryExtend(line)
            if extended:
                cnt += 1
                stop = False
            if cnt > 20:
                cnt = 0
                saveAll(save_handler.getWriter(), params, aligner, contigs, reads, disjointigs, lines, dot_plot)
        if stop:
            break

    lines.printToFile(sys.stdout)
    lines.printToFasta(open(os.path.join(params.dir, "lines.fasta"), "w"))
    sys.stdout.write("Finished\n")
    print (time.strftime("%d.%m.%Y  %I:%M:%S"))


if __name__ == "__main__":
    main(sys.argv)