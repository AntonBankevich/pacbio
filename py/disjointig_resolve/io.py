from typing import Tuple

from alignment.align_tools import Aligner
from common.save_load import TokenReader, TokenWriter
from common.sequences import ContigCollection, ReadCollection
from disjointig_resolve.accurate_line import NewLineStorage
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.params import Params


def loadAll(handler):
    # type: (TokenReader) -> Tuple[Params, Aligner, ContigCollection, ReadCollection, DisjointigCollection, NewLineStorage, LineDotPlot]
    params = Params()
    params.load(handler)
    aligner = Aligner.load(handler)
    contigs = ContigCollection()
    contigs.load(handler)
    reads = ReadCollection()
    reads.loadFromFasta(open(params.reads_file, "r"))
    disjointigs = DisjointigCollection()
    disjointigs.load(handler, reads)
    lines = NewLineStorage(disjointigs)
    lines.load(handler, reads, contigs)
    dot_plot = LineDotPlot(lines)
    dot_plot.load(handler)
    return params, aligner, contigs, reads, disjointigs, lines, dot_plot


def saveAll(handler, params, aligner, contigs, reads, disjointigs, lines, dot_plot):
    # type: (TokenWriter, Params, Aligner, ContigCollection, ReadCollection, DisjointigCollection, NewLineStorage, LineDotPlot) -> None
    params.save(handler)
    aligner.save(handler)
    contigs.save(handler)
    disjointigs.save(handler)
    lines.save(handler)
    dot_plot.save(handler)