import itertools
import os
import sys

from common.save_load import TokenWriter, TokenReader

sys.path.append("py")
from common.seq_records import NamedSequence
from flye_tools.alignment import make_alignment
from common.sequences import ContigCollection, Segment, Contig, \
    ContigStorage
from common.alignment_storage import AlignmentPiece, ReadCollection
from typing import Iterable, Tuple, Generator, BinaryIO
from common import basic, sam_parser, SeqIO, params


class DirDistributor:
    def __init__(self, dir):
        self.dir = dir
        self.cur_dir = 0

    def nextDir(self):
        # type: () -> str
        name = os.path.join(self.dir, str(self.cur_dir))
        self.cur_dir += 1
        assert self.cur_dir <= 10000
        basic.ensure_dir_existance(name)
        return name

    def CheckSequences(self, reads, reads_file):
        # type: (Iterable[NamedSequence], str) -> bool
        if not os.path.exists(reads_file):
            return False
        try:
            for rec, read in itertools.izip_longest(SeqIO.parse_fasta(open(reads_file, "r")), reads):
                if str(rec.id) != str(read.id) or rec.seq != read.seq:
                    return False
            return True
        except:
            return False

    def WriteSequences(self, reads, reads_file):
        # type: (Iterable[NamedSequence], str) -> None
        f = open(reads_file, "w")
        for read in reads:
            SeqIO.write(read, f, "fasta")
        f.close()

    def hash(self, reads):
        # type: (Iterable[NamedSequence]) -> int
        res = 0
        for read in reads:
            res += read.seq.__hash__() + read.id.__hash__()
        return res

    def calculateHash(self, content):
        # type: (list[Tuple[Iterable[NamedSequence], str]]) -> Generator[Tuple[str, str, str]]
        for reads, f_name in content:
            yield f_name, str(self.hash(reads)), str(len(reads))

    def printHash(self, handler, hashs):
        # type: (BinaryIO, list[Tuple[str, str, str]]) -> None
        for rec in hashs:
            handler.write(" ".join(rec) + "\n")

    def compareHash(self, handler, hashs):
        # type: (BinaryIO, list[Tuple[str, str, str]]) -> bool
        lines = handler.readlines()
        if len(lines) != len(hashs):
            return False
        for l, rec in zip(lines, hashs):
            l = l.split()
            if len(l) != len(rec):
                return False
            for s1, s2 in zip(l, rec):
                if s1 != s2:
                    return False
        return True

    def fillNextDir(self, content):
        # type: (list[Tuple[Iterable[NamedSequence], str]]) -> Tuple[str, list[str], bool]
        same = True
        dir = self.nextDir()
        content_files = []
        for reads, f_name in content:
            content_files.append(os.path.join(dir, f_name))
        hash_file = os.path.join(dir, "hashs.txt")
        hashs = list(self.calculateHash(content))
        if os.path.isfile(hash_file) and self.compareHash(open(hash_file, "r"), hashs) and not params.redo_alignments:
            return dir, content_files, True
        self.printHash(open(hash_file, "w"), hashs)
        for reads, f_name in content:
            f_name = os.path.join(dir, f_name)
            self.WriteSequences(reads, f_name)
        return dir, content_files, False

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.dir)
        handler.writeIntLine(self.cur_dir)

    @staticmethod
    def load(handler):
        # type: (TokenReader) -> DirDistributor
        res = DirDistributor(handler.readToken())
        res.cur_dir = handler.readInt()
        return res


class Aligner:
    def __init__(self, dir_distributor, threads = 16):
        # type: (DirDistributor, int) -> Aligner
        self.dir_distributor = dir_distributor
        self.threads = threads

    def alignReadCollection(self, reads_collection, contigs):
        # type: (ReadCollection, Iterable[Contig]) -> None
        contig_collection = ContigCollection(contigs)
        contig_ids = set()
        for contig in contigs:
            if contig.rc.id not in contig_ids:
                contig_ids.add(contig.id)
        read_ids = set()
        for read in reads_collection:
            if read.rc.id not in read_ids:
                read_ids.add(read.id)
        contigs = filter(lambda contig: contig.id in contig_ids, contigs)
        reads = filter(lambda read: read.id in read_ids, reads_collection)
        reads_collection.fillFromSam(self.align(reads, contigs), contig_collection)

    def alignReadsToSegments(self, reads, segments):
        # type: (ReadCollection, Iterable[Segment]) -> None
        segments = list(segments)
        seg_dict = dict()
        for i, seg in enumerate(segments):
            seg_dict[str(i + 1)] = seg
        contigs = map(lambda (i, seg): Contig(seg.Seq(), str(i + 1)), enumerate(segments))
        read_collection = ReadCollection().extendClean(reads)
        self.alignReadCollection(read_collection, ContigCollection(contigs))
        read_collection.contigsAsSegments(seg_dict)
        reads.mergeAlignments(read_collection)

    def realignCollection(self, reads_collection, contigs):
        # type: (ReadCollection, Iterable[Contig]) -> None
        for read in reads_collection:
            read.clean()
        self.alignReadCollection(reads_collection, contigs)

    # def fixExtendedLine(self, line):
    #     # type: (Line) -> None
    #     toFix = []
    #     for read in line.reads:
    #         if not read.noncontradicting(line.asSegment()):
    #             toFix.append(read)
    #     newAlignments = ReadCollection(line.reads.contigs).extend(toFix)
    #     self.alignReadCollection(newAlignments)
    #     for read in newAlignments:
    #         for al in line.reads[read.id].alignments:
    #             if not al.contradicting(line.asSegment()):
    #                 continue
    #             for new_al in read.alignments:
    #                 if new_al.contains(al) and len(new_al) > len(al):
    #                     al.seg_from = new_al.seg_from
    #                     al.seg_to = new_al.seg_to
    #                     al.cigar = new_al.cigar

    def separateAlignments(self, reads, contigs):
        # type: (Iterable[NamedSequence], Iterable[Contig]) -> ReadCollection
        contigs_collection = ContigCollection(list(contigs))
        res = ReadCollection(contigs_collection)
        for read in reads:
            res.addNewRead(NamedSequence(read.seq, read.id)) # remove when all ids are str
        for contig in contigs:
            res.fillFromSam(self.align(res, [contig]), contigs_collection)
        return res

    def align(self, reads, reference):
        # type: (Iterable[NamedSequence], Iterable[Contig]) -> sam_parser.Samfile
        dir, new_files, same = self.dir_distributor.fillNextDir([(reference, "contigs.fasta"), (list(reads), "reads.fasta")])
        contigs_file = new_files[0]
        reads_file = new_files[1]
        alignment_dir = os.path.join(dir, "alignment")
        alignment_file = os.path.join(dir, "alignment.sam")
        basic.ensure_dir_existance(dir)
        basic.ensure_dir_existance(alignment_dir)
        if same and not params.clean and os.path.exists(alignment_file):
            # print "Alignment reused:", alignment_file
            pass
        else:
            # print "Performing alignment:", alignment_file
            make_alignment(contigs_file, [reads_file], self.threads, alignment_dir, "pacbio", alignment_file)
        return sam_parser.Samfile(open(alignment_file, "r"))

    # TODO: make this method accept reference dict of dome sort
    def alignClean(self, reads, ref_storage):
        # type: (Iterable[Contig], ContigStorage) -> Generator[AlignmentPiece]
        reads = list(reads)
        # print "Aligning", reads
        # print list(ref_storage)
        parser = self.align(reads, list(ref_storage.unique()))
        read_storage = ContigStorage(reads, False)
        for rec in parser:
            if rec.is_unmapped:
                continue
            rname = rec.query_name
            seq_from = read_storage[rname]
            seq_to = ref_storage[rec.tname]
            for al in AlignmentPiece.FromSamRecord(seq_from, seq_to, rec).split():
                yield al


    # def matchingAlignment(self, seqs, contig):
    #     # type: (list[str], Contig) -> list[AlignedSequences]
    #     collection = ContigCollection([contig])
    #     res = [] # type: list[AlignedSequences]
    #     for seq in seqs:
    #         res.append(AlignedSequences(seq, contig.seq))
    #     reads = ReadCollection(collection).loadFromSam(
    #         self.align([AlignedRead(SeqIO.SeqRecord(seq, str(i))) for i, seq in enumerate(seqs)], collection))
    #     for read in reads:
    #         tid = int(read.id)
    #         read.sort()
    #         groups = [] #type: list[list[AlignmentPiece]]
    #         group_lens = []
    #         for al in read.alignments:
    #             if al.seg_to.contig != contig:
    #                 continue
    #             found = False
    #             for i, group in enumerate(groups):
    #                 if group[-1].precedes(al, 50):
    #                     group.append(al)
    #                     group_lens[i] += len(al.seg_from)
    #                     found = True
    #                     break
    #             if not found:
    #                 groups.append([al])
    #                 group_lens.append(len(al.seg_from))
    #         best = None
    #         for i in range(len(groups)):
    #             if best == None or group_lens[i] > group_lens[best]:
    #                 best = i
    #         for al in groups[best]:
    #             res[tid].addCigar(al.cigar, al.seg_to.left)
    #     return res

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeIntLine(self.threads)
        self.dir_distributor.save(handler)

    @staticmethod
    def load(handler):
        # type: (TokenReader) -> Aligner
        return Aligner(DirDistributor.load(handler), handler.readInt())


if __name__ == "__main__":
    dir = sys.argv[1]
    query = sys.argv[2]
    target = sys.argv[3]
    aln = Aligner(dir)
    contigs = ContigCollection().loadFromFasta(open(target, "r"))
    for al in aln.alignClean(ReadCollection().loadFromFasta(open(query, "r")), contigs):
        print al



