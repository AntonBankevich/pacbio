from typing import Dict, List, Optional, BinaryIO, Callable, Iterator, Generator, Iterable, Union
from common import basic, SeqIO
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, UniqueList, ReadCollection, Contig
from common.alignment_storage import AlignmentPiece, AlignedRead
from disjointig_resolve.smart_storage import AlignmentStorage


class Disjointig(Contig):
    def __init__(self, seq, id, rc=None):
        # type: (str, str, Optional[Disjointig]) -> None
        self.seq = seq
        self.id = id
        if rc is None:
            self.read_alignments = AlignmentStorage() # type: AlignmentStorage
            rc = Disjointig(basic.RC(seq), basic.Reverse(id), self) # type: Disjointig
        else:
            self.read_alignments = self.rc.read_alignments.rc  # type: AlignmentStorage
        Contig.__init__(self, seq, id, rc)
        self.rc = rc # type:Disjointig

    def addAlignments(self, als):
        # type: (Iterable[AlignmentPiece]) -> None
        self.read_alignments.addAll(als)

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.read_alignments.add(al)

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        return self.read_alignments.getAlignmentsTo(seg)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)
        self.read_alignments.save(handler)

    def loadDisjointig(self, handler, disjointigs, reads):
        # type: (TokenReader, DisjointigCollection, ReadCollection) -> None
        self.id = handler.readToken()
        self.rc.id = basic.RC(self.id)
        self.read_alignments.load(handler, reads, disjointigs)
        for al in self.read_alignments:
            read = al.seg_from.contig # type: AlignedRead
            read.addAlignment(al)


class UniqueMarker:
    def __init__(self, disjointigs):
        # type: (DisjointigCollection) -> None
        self.disjointigs = disjointigs

    # Mark unique regions on all disjointigs as correct
    def findUnique(self, contig):
        # type: (NamedSequence) -> Generator[Segment]
        # IMPLEMENT search for unique segments in disjointigs including construction of dot plot
        pass

    def findAllUnique(self, sequences):
        # type: (Iterable[NamedSequence]) -> Generator[Segment]
        for disjointig in UniqueList(sequences):
            for seg in self.findUnique(disjointig):
                yield seg

class DisjointigCollection:
    def __init__(self):
        self.disjointigs = dict() # type: Dict[str, Disjointig]
        self.cnt = 1

    def __iter__(self):
        # type: () -> Iterator[Disjointig]
        return self.disjointigs.values().__iter__()

    def __getitem__(self, item):
        # type: (str) -> Disjointig
        return self.disjointigs[item]

    def __contains__(self, item):
        # type: (Disjointig) -> bool
        return item.id in self.disjointigs

    def containsKey(self, key):
        # type: (str) -> bool
        return key in self.disjointigs

    def add(self, seq, name = None):
        # type: (str, Optional[str]) -> Disjointig
        if name is None:
            name = "D" + str(self.cnt)
            self.cnt += 1
        new_disjointig = Disjointig(seq, name)
        self.disjointigs[name] = new_disjointig
        self.disjointigs[new_disjointig.rc.id] = new_disjointig.rc
        return new_disjointig

    def loadFromFasta(self, handler, save_names = False, int_ids = False, filter = lambda rec: True):
        # type: (BinaryIO, bool, bool, Callable[[NamedSequence], bool]) -> DisjointigCollection
        recs = list(SeqIO.parse_fasta(handler))
        if save_names:
            for rec in recs:
                assert rec.id not in self.disjointigs.keys() and basic.Reverse(rec.id) not in self.disjointigs.keys()
        for rec in recs:
            if not filter(rec):
                continue
            if save_names:
                number = basic.parseNegativeNumber(rec.id)
                if number is not None:
                    self.cnt = max(self.cnt, int(abs(number)) + 1)
                if int_ids:
                    self.add(rec.seq, str(number))
                else:
                    self.add(rec.seq, rec.id)
            else:
                self.add(rec.seq)

    def addAll(self, als):
        # type: (Union[Generator[AlignmentPiece], Iterable[AlignmentPiece]]) -> None
        for al in als:
            dt = al.seg_from.contig # type: Disjointig
            dt.addAlignment(al)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(str(self.cnt))
        handler.writeTokens(str(map(lambda line: line.id, UniqueList(self.disjointigs.values()))))
        for disjointig in UniqueList(self.disjointigs.values()):
            handler.writeTokenLine(disjointig.id)
            handler.writeTokenLine(disjointig.seq)
        for disjointig in UniqueList(self.disjointigs.values()):
            disjointig.save(handler)

    def load(self, handler, reads):
        # type: (TokenReader, ReadCollection) -> None
        self.cnt = int(handler.readToken())
        keys = handler.readTokens()
        for key in keys:
            id = handler.readToken()
            assert id == key
            seq = handler.readToken()
            disjointig = self.add(seq, key)
            assert key == disjointig.id, key + " " + disjointig.id
        for key in keys:
            self.disjointigs[key].loadDisjointig(handler, self, reads)



