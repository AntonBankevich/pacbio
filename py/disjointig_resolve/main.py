import os
import sys


sys.path.append("py")

from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.sequences import ContigCollection, ReadCollection
from disjointig_resolve.accurate_line import NewLineStorage
from disjointig_resolve.io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.params import Params


def main(args):
    params = Params().parse(args)
    # IMPLEMENT add logger from old main
    print "Preparing initial state"
    if params.load_from is not None:
        print "Loading initial state from saves"
        params, contigs, reads, disjointigs, lines = loadAll(TokenReader(open(params.load_from, "r")))
    else:
        print "Creating disjointig collection"
        disjointigs = DisjointigCollection()
        disjointigs.loadFromFasta(open(params.disjointigs_file, "r"))
        disjointigs.calculateDotPlot()

        print "Creating read collection"
        reads = ReadCollection()
        reads.loadFromFasta(open(params.reads_file, "r"))

        print "Creating contig collection"
        contigs = ContigCollection()
        contigs.loadFromFasta(open(params.contigs_file, "r"), num_names=True)

        print "Creating line collection"
        lines = NewLineStorage(disjointigs)
        lines.fillFromContigs(contigs)
        lines.fillFromDisjointigs()

    save_handler = SaveHandler(params.save_dir)
    print "Resolving"
    extender = LineExtender(disjointigs)
    knotter = LineKnotter(lines)
    cnt = 0
    while True:
        stop = True
        for line in lines:
            extended = extender.tryExtend(line)
            if extended:
                cnt += 1
                stop = False
                knotter.tryKnot(line)
            if cnt > 50:
                cnt = 0
                saveAll(save_handler.getWriter(), params, contigs, reads, disjointigs, lines)
        if stop:
            break

    lines.printToFile(sys.stdout)
    lines.printToFasta(open(os.path.join(params.dir, "lines.fasta"), "w"))




if __name__ == "__main__":
    main(sys.argv)