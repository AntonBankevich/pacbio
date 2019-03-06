from common.sequences import ContigCollection, ReadCollection
from disjointig_resolve.accurate_line import NewLineStorage
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