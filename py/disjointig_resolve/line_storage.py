from typing import Dict, List, Iterator, Optional, Iterable, BinaryIO

from alignment.align_tools import Aligner
from common import SeqIO
from common.alignment_storage import AlignmentPiece, ReadCollection
from common.save_load import TokenWriter, TokenReader
from common.sequences import ContigStorage, UniqueList, Contig, ContigCollection
from disjointig_resolve.accurate_line import NewLine, ExtensionHandler
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.line_alignments import TwoLineAlignmentStorage
from disjointig_resolve.smart_storage import AlignmentStorage


class NewLineStorage(ContigStorage):
    def __init__(self, disjointigs, aligner):
        # type: (DisjointigCollection, Aligner) -> None
        ContigStorage.__init__(self, [], True)
        self.disjointigs = disjointigs
        self.aligner = aligner
        self.items = dict() # type: Dict[str, NewLine]
        self.cnt = 1
        self.listeners = [] # type: List[LineStorageListener]
        self.name_printer = None

    def __iter__(self):
        # type: () -> Iterator[NewLine]
        return self.items.values().__iter__()

    def __getitem__(self, item):
        # type: (str) -> NewLine
        return self.items[item]

    def addListener(self, listener):
        self.listeners.append(listener)

    def removeListener(self, listener):
        self.listeners.remove(listener)

    def notifyMergedLines(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        for listener in self.listeners:
            listener.FireMergedLines(al1, al2)

    def addNew(self, seq, name = None):
        # type: (str, Optional[str]) -> NewLine
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        new_line = NewLine(seq, name, ExtensionHandler(self.disjointigs, self.aligner))
        self.add(new_line)
        new_line.name_printer = self.name_printer
        new_line.rc.name_printer = self.name_printer
        return new_line

    def fillFromContigs(self, contigs):
        # type: (Iterable[Contig]) -> None
        for contig in UniqueList(contigs):
            line = self.addNew(contig.seq)
            line.initial.add(AlignmentPiece.Identical(line.asSegment(), contig.asSegment()))
            for seg in line.correct_segments:
                line.completely_resolved.add(seg)

    def alignDisjointigs(self):
        for line in self:
            line.disjointig_alignments.clean()
        for al in self.aligner.alignClean(self.disjointigs.unique(), self):
            line = al.seg_to.contig # type: NewLine
            line.disjointig_alignments.add(al)

    # def fillFromDisjointigs(self):
    #     # type: () -> None
    #     for seg in UniqueMarker(self.disjointigs).findAllUnique(self.disjointigs):
    #         line = self.addNew(seg.Seq())
    #         line.initial.add(AlignmentPiece.Identical(line.asSegment(), seg))
    #     #TODO Filter all lines already present in the collection

    def mergeLines(self, alignment, k):
        # type: (AlignmentPiece, int) -> NewLine
        line1 = alignment.seg_from.contig #type: NewLine
        line2 = alignment.seg_to.contig #type: NewLine
        assert line1 != line2

        # Cutting hanging tips of both lines
        al_storage = AlignmentStorage()
        al_storage.add(alignment)
        storage = TwoLineAlignmentStorage(line1, line2)
        line2.addListener(storage)
        line1.addListener(storage.reverse)
        storage.add(alignment)
        if alignment.seg_from.right < len(line1):
            line1.cutRight(alignment.seg_from.right)
        if alignment.seg_to.left > 0:
            line2.rc.cutRight(len(line2) - alignment.seg_to.left)
        alignment = list(storage.content)[0] # type: AlignmentPiece
        line2.removeListener(storage)
        line1.removeListener(storage.reverse)

        # Making sure line sequences match on the overlap
        new_seq = Contig(line1.asSegment().prefix(pos=alignment.seg_from.left).Seq() + line2.seq, "new_seq")
        al2 = AlignmentPiece.Identical(line2.asSegment(), new_seq.asSegment().suffix(length=len(line2)))
        alignment = alignment.compose(al2).reverse()
        assert alignment.seg_to.right == len(line1)
        assert alignment.seg_from.left == al2.seg_to.left
        line1.correctSequence([alignment])

        # Now lines have exact match
        name = "(" + line1.id + "," + line2.id + ")"
        line = self.addNew(new_seq.seq, name)
        assert line.seq.startswith(line1.seq)
        assert line.seq.endswith(line2.seq)
        al1 = AlignmentPiece.Identical(line1.asSegment(), line.asSegment().prefix(length=len(line1)))
        al2 = AlignmentPiece.Identical(line2.asSegment(), line.asSegment().suffix(length=len(line2)))

        line.initial.addAll(line1.initial.targetAsSegment(al1.seg_to).merge(line2.initial.targetAsSegment(al2.seg_to)))
        line.correct_segments.addAll(line1.correct_segments.contigAsSegment(al1.seg_to).
                                     merge(line2.correct_segments.contigAsSegment(al2.seg_to)))
        line.completely_resolved.addAll(line1.completely_resolved.contigAsSegment(al1.seg_to).
                                     merge(line2.completely_resolved.contigAsSegment(al2.seg_to), k))
        line.disjointig_alignments.addAll(line1.disjointig_alignments.targetAsSegment(al1.seg_to).
                                          merge(line2.disjointig_alignments.targetAsSegment(al2.seg_to)))
        for al in line1.read_alignments.targetAsSegment(al1.seg_to).merge(line2.read_alignments.targetAsSegment(al2.seg_to)):
            line.addReadAlignment(al)
        line1.cleanReadAlignments()
        line2.cleanReadAlignments()


        # print line1.read_alignments
        # print line2.read_alignments
        # print line.read_alignments
        self.notifyMergedLines(al1, al2)
        self.remove(line1)
        self.remove(line2)
        return line

    def printToFile(self, handler):
        # type: (BinaryIO) -> None
        for line in self.items:
            handler.write(line.__str__() + "\n")

    def printToFasta(self, handler):
        # type: (BinaryIO) -> None
        for line in UniqueList(self.items.values()):
            SeqIO.write(line, handler, "fasta")

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(str(self.cnt))
        line_ids = map(lambda line: line.id, UniqueList(self.items.values()))
        handler.writeTokens(line_ids)
        for line_id in line_ids:
            line = self.items[line_id]
            handler.writeTokenLine(line.seq)
        for line_id in line_ids:
            line = self.items[line_id]
            line.save(handler)

    def load(self, handler, reads, contigs):
        # type: (TokenReader, ReadCollection, ContigCollection) -> None
        self.cnt = int(handler.readToken())
        keys = handler.readTokens()
        for key in keys:
            self.addNew(handler.readToken(), key)
        for key in keys:
            line = self.items[key]
            line.loadLine(handler, self.disjointigs, reads, contigs)

    def remove(self, line):
        del self.items[line.id]
        del self.items[line.rc.id]


class LineStorageListener:
    def __init__(self):
        pass

    def FireMergedLines(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        pass