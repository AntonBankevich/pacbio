import os
import sys

from typing import Tuple

from common.save_load import TokenReader, TokenWriter, SaveHandler
from common.sequences import ContigCollection, UniqueList, AlignmentPiece, ReadCollection
from disjointig_resolve.accurate_line import NewLineStorage

sys.path.append("py")

from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.params import Params

def loadAll(handler):
    # type: (TokenReader) -> Tuple[Params, ContigCollection, ReadCollection, DisjointigCollection, NewLineStorage]
    params = Params()
    params.load(handler)
    contigs = ContigCollection()
    contigs.load(handler)
    reads = ReadCollection()
    reads.loadFromFasta(open(params.reads_file, "r"))
    disjointigs = DisjointigCollection()
    disjointigs.load(handler, reads)
    lines = NewLineStorage(disjointigs)
    lines.load(handler, reads, contigs)
    return params, contigs, reads, disjointigs, lines

def saveAll(handler, params, contigs, reads, disjointigs, lines):
    # type: (TokenWriter, Params, ContigCollection, ReadCollection, DisjointigCollection, NewLineStorage) -> None
    params.save(handler)
    contigs.save(handler)
    disjointigs.save(handler)
    lines.save(handler)


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
    saveAll(save_handler.getWriter(), params, contigs, reads, disjointigs, lines)
    print "Resolving"
    extender = LineExtender(disjointigs)
    knotter = LineKnotter(lines)
    while True:
        stop = True
        for line in lines:
            extended = extender.tryExtend(line)
            if extended:
                stop = False
                knotter.tryKnot(line)
        if stop:
            break

    lines.printToFile(sys.stdout)
    lines.printToFasta(open(os.path.join(params.dir, "lines.fasta"), "w"))





#     IMPLEMENT











if __name__ == "__main__":
    main(sys.argv)