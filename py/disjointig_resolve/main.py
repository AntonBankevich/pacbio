import os
import sys

from alignment.align_tools import Aligner, DirDistributor
from disjointig_resolve.dot_plot import LineDotPlot

sys.path.append("py")

from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.sequences import ContigCollection, ReadCollection
from disjointig_resolve.accurate_line import NewLineStorage
from disjointig_resolve.io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.cl_params import Params


def main(args):
    params = Params().parse(args)
    # IMPLEMENT add logger from old main
    print "Preparing initial state"
    if params.load_from is not None:
        print "Loading initial state from saves"
        params, aligner, contigs, reads, disjointigs, lines, dot_plot = loadAll(TokenReader(open(params.load_from, "r")))
    else:
        aligner = Aligner(DirDistributor(params.alignmentDir()))

        print "Creating disjointig collection"
        disjointigs = DisjointigCollection()
        disjointigs.loadFromFasta(open(params.disjointigs_file, "r"))

        print "Creating read collection"
        reads = ReadCollection()
        reads.loadFromFasta(open(params.reads_file, "r"))

        print "Aligning reads to disjointigs"
        disjointigs.addAll(aligner.alignClean(reads, disjointigs))

        print "Creating contig collection"
        contigs = ContigCollection()
        contigs.loadFromFasta(open(params.contigs_file, "r"), num_names=True)

        print "Creating line collection"
        lines = NewLineStorage(disjointigs)
        lines.fillFromContigs(contigs)
        lines.fillFromDisjointigs()

        dot_plot = LineDotPlot(lines)
        dot_plot.construct(aligner)

    save_handler = SaveHandler(params.save_dir)
    print "Resolving"
    extender = LineExtender(disjointigs)
    knotter = LineKnotter(lines)#IMPLEMENT
    cnt = 0
    while True:
        stop = True
        for line in lines:
            extended = extender.tryExtend(line)
            if extended:
                cnt += 1
                stop = False
                knotter.tryKnot(line)
            if cnt > 20:
                cnt = 0
                saveAll(save_handler.getWriter(), params, aligner, contigs, reads, disjointigs, lines, dot_plot)
        if stop:
            break

    lines.printToFile(sys.stdout)
    lines.printToFasta(open(os.path.join(params.dir, "lines.fasta"), "w"))




if __name__ == "__main__":
    main(sys.argv)