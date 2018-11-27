import itertools
import os
import sys

sys.path.append("py")
from common.seq_records import NamedSequence
from flye.alignment import make_alignment
from dag_resolve import params
from dag_resolve.repeat_graph import Graph
from common.sequences import AlignedRead, Contig, ContigCollection, ReadCollection, Segment
from typing import Optional, Iterable, Tuple, Generator
from common import basic, sam_parser, SeqIO


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
        # type: (file, list[Tuple[str, str, str]]) -> None
        for rec in hashs:
            handler.write(" ".join(rec) + "\n")

    def compareHash(self, handler, hashs):
        # type: (file, list[Tuple[str, str, str]]) -> bool
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
        if os.path.isfile(hash_file) and self.compareHash(open(hash_file, "r"), hashs):
            return dir, content_files, True
        self.printHash(open(hash_file, "w"), hashs)
        for reads, f_name in content:
            f_name = os.path.join(dir, f_name)
            self.WriteSequences(reads, f_name)
        return dir, content_files, False

class AlignedSequences:
    def __init__(self, seq_from, seq_to):
        # type: (str, str) -> AlignedSequences
        self.seq_from = seq_from
        self.seq_to = seq_to
        self.alignment = [-1] * len(seq_to)
        self.first = None # type: int
        self.last = None # type: int
        self.aligned = False # type: bool
        self.rc = False # type: bool

    def alignedSequence(self):
        return self.seq_from[self.alignment[self.first]:self.alignment[self.last] + 1]

    def AlignedSegments(self):
        positions = []
        res = ""
        for i in range(self.first, self.last + 1):
            if self.alignment[i] == -1:
                continue
            prev = [(pos, self.alignment[pos]) for pos in range(max(i - 20, 0), i + 1) if self.alignment[pos] != -1]
            if len(prev) >= 15 and abs(prev[-1][1] - prev[0][1] - 20) < 5:
                positions.append(i)
        start = self.first
        end = self.first
        for pos in positions:
            if start == -1:
                start = pos
                end = pos
            if pos - end > 5:
                res +=  str((self.alignment[start], self.alignment[pos])) + "->" + str((start, pos)) + "\n"
                start = -1
            else:
                end = pos
        res += str((self.alignment[start], self.alignment[self.last])) + "->" + str((start, self.last)) + "\n"
        return res

    def percentIdentity(self):
        assert self.aligned
        match = 0
        for i in range(self.first, self.last + 1):
            if self.alignment[i] == -1:
                continue
            if self.seq_from[self.alignment[i]] == self.seq_to[i]:
                match += 1
        return float(match) / max(self.last - self.first + 1, self[self.last] - self[self.first] + 1)

    def printToFile(self, handler):
        handler.write("Aligned sequences: " + str(self.first) + " " + str(self.last) + " " +
                      str(self.alignment[self.first]) + " " + str(self.alignment[self.last]) + "\n")
        handler.write(self.seq_to[self.first:self.last + 1] + "\n")
        handler.write(self.seq_from[self.alignment[self.first]:self.alignment[self.last] + 1] + "\n")

    def __getitem__(self, item):
        # type: (int) -> int
        return self.alignment[item]

    def addCigar(self, cigar, pos):
        # type: (str, int) -> None
        for seq_pos, contig_pos in sam_parser.ParseCigar(cigar, len(self.seq_to), pos - 1, True):
            # if contig_pos >= len(self.alignment):
            #     print contig_pos, len(self.alignment), self.rc, pos
            #     print cigar
            #     print list(sam_parser.ParseCigar(cigar, len(self.seq_to), pos))
            self.alignment[contig_pos] = seq_pos
            self.updateBounds(contig_pos)

    def updateBounds(self, pos):
        # type: (int) -> None
        if self.first is None or self.first > pos:
            self.first = pos
        if self.last is None or self.last < pos:
            self.last = pos
            self.aligned = True

    def findNextMatch(self, pos):
        # type: (int) -> Optional[int]
        if not self.aligned:
            return None
        if pos <= self.first:
            return self.first
        while pos <= self.last and self.alignment[pos] == -1:
            pos += 1
        if pos <= self.last:
            return pos
        else:
            return None

    def findPreviousMatch(self, pos):
        # type: (int) -> Optional[int]
        if not self.aligned:
            return None
        if pos >= self.last:
            return self.last
        while pos >= self.first and self.alignment[pos] == -1:
            pos -= 1
        if pos >= self.first:
            return pos
        else:
            return None

    def mapSegment(self, l, r):
        # type: (int, int) -> Optional[Tuple[int, int]]
        assert self.aligned
        if l < self.first:
            return None
        if r > self.last:
            return None
        l_match = self.findNextMatch(l)
        r_match = self.findPreviousMatch(r)
        if l_match > r_match:
            return None
        if l <= l_match and l_match <= r_match and r_match <= r:
            pass
        else:
            print l, l_match, r_match, r, self.first, self.last
            for i in range(min([l, l_match, r_match, r], 1 + max([l, l_match, r_match, r]))):
                print i, self.alignment[i]
        assert l <= l_match and l_match <= r_match and r_match <= r
        return self.alignment[l_match] - (l_match - l), self.alignment[r_match] + (r - r_match)

    def checkHomo(self, pos):
        if pos == 0 or self.seq_to[pos] != self.seq_to[pos -1]:
            return False
        if self.alignment[pos] != -1 and self.alignment[pos - 1] != -1:
            for i in xrange(self.alignment[pos - 1], self.alignment[pos] + 1):
                if self.seq_from[i] != self.seq_to[pos]:
                    return False
        return True

    def divPositions(self, first, last):
        homo_start = first
        homo_active = False
        for i in xrange(first, last):
            if not self.checkHomo(i):
                homo_start = i
                homo_active = False
            elif homo_active:
                yield i
                continue
            if self.alignment[i] == -1 or \
                    (self.alignment[i - 1] != -1 and self.alignment[i] != self.alignment[i - 1] + 1) or \
                            self.seq_to[i] != self.seq_from[self.alignment[i]]:
                if self.checkHomo(i):
                    for j in xrange(homo_start, i + 1):
                        yield j
                else:
                    yield i

    def findSeqPos(self, pos):
        # type: (int) -> Optional[Tuple[Optional[int], Optional[int]]]
        if not self.aligned:
            return None
        if self.alignment[self.last] < pos:
            return (self.last, None)
        if self.alignment[self.first] > pos:
            return (None, self.first)
        l = self.first
        r = self.last
        while l + 1 < r:
            m = (l + r) // 2
            if self.alignment[self.findPreviousMatch(m)] <= pos:
                l = m
            else:
                r = m
        if self.alignment[l] == pos:
            return (l, l)
        if self.alignment[r] == pos:
            return (r, r)
        return (self.findPreviousMatch(l), self.findNextMatch(r))

    def checkPosEquals(self, pos):
        assert pos >= 0 and pos < len(self.seq_to)
        return self.alignment[pos] != -1 and self.seq_from[self.alignment[pos]] == self.seq_to[pos]

    def checkNoIndelAfter(self, pos):
        assert pos >= 0 and pos < len(self.seq_to)
        if pos == len(self.seq_to) - 1:
            return True
        return self.alignment[pos] != -1 and self.alignment[pos + 1] != -1 and self.alignment[pos] + 1 == self.alignment[pos + 1]


def ReadToAlignedSequences(read, contig):
    # type: (AlignedRead, Contig) -> AlignedSequences
    res = AlignedSequences(read.seq, contig.seq)
    for rec in read.alignments:
        if rec.seg_to.contig != contig:
            continue
        res.addCigar(rec.cigar, rec.seg_to.left)
    return res




class Aligner:
    def __init__(self, dir_distributor, threads = 16):
        # type: (DirDistributor, int) -> Aligner
        self.dir_distributor = dir_distributor
        self.cur_alignment = 0
        self.threads = threads

    def alignReadCollection(self, reads_collection, contigs = None):
        # type: (ReadCollection, Iterable[Contig]) -> None
        if contigs is None:
            contigs = reads_collection.contigs
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
        reads_collection.fillFromSam(self.align(reads, contigs))

    def alignReadsToSegments(self, reads, segments):
        # type: (ReadCollection, Iterable[Segment]) -> None
        segments = list(segments)
        seg_dict = dict()
        for i, seg in enumerate(segments):
            seg_dict[str(i + 1)] = seg
        contigs = map(lambda (i, seg): Contig(seg.Seq(), str(i + 1)), enumerate(segments))
        read_collection = ReadCollection(ContigCollection(contigs)).extendClean(reads)
        self.alignReadCollection(read_collection)
        read_collection.contigsAsSegments(seg_dict)
        reads.mergeAlignments(read_collection)

    def realignCollection(self, reads_collection):
        # type: (ReadCollection) -> None
        for read in reads_collection:
            read.clean()
        self.alignReadCollection(reads_collection)

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

    def expandCollection(self, reads_collection, new_reads):
        # type: (ReadCollection, list[AlignedRead]) -> None
        for read in new_reads:
            reads_collection.addNewRead(read)
        reads_collection.fillFromSam(self.align(new_reads, reads_collection.contigs))

    def separateAlignments(self, reads, contigs):
        # type: (Iterable[NamedSequence], Iterable[Contig]) -> ReadCollection
        res = ReadCollection(ContigCollection(list(contigs)))
        for read in reads:
            res.addNewRead(NamedSequence(read.seq, read.id)) # remove when all ids are str
        for contig in contigs:
            res.fillFromSam(self.align(res, [contig]))
        return res

    def align(self, reads, reference):
        # type: (Iterable[NamedSequence], Iterable[Contig]) -> sam_parser.Samfile
        dir, new_files, same = self.dir_distributor.fillNextDir([(reference, "contigs.fasta"), (reads, "reads.fasta")])
        contigs_file = new_files[0]
        reads_file = new_files[1]
        alignment_dir = os.path.join(dir, "alignment")
        alignment_file = os.path.join(dir, "alignment.sam")
        basic.ensure_dir_existance(dir)
        basic.ensure_dir_existance(alignment_dir)
        if same and not params.clean and os.path.exists(alignment_file):
            print "Alignment reused:", alignment_file
        else:
            print "Performing alignment:", alignment_file
            make_alignment(contigs_file, [reads_file], self.threads, alignment_dir, "pacbio", alignment_file)
        return sam_parser.Samfile(open(alignment_file, "r"))

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

    def repairGraphAlignments(self, graph):
        # type: (Graph) -> None
        print "Reparing graph alignments for", len(graph.reads), "reads and the following edges:", map(lambda edge: edge.id, graph.newEdges)
        for rec in self.align(graph.reads, graph.newEdges):
            if not rec.is_unmapped:
                read = graph.reads[rec.query_name]
                graph.E[int(rec.tname)].reads.add(read)
                graph.E[int(rec.tname)].reads.addNewAlignment(rec)
        graph.newEdges = []


if __name__ == "__main__":
    dir = sys.argv[1]
    query = sys.argv[2]
    target = sys.argv[3]
    aln = Aligner(dir)
    contigs = ContigCollection().loadFromFasta(open(target, "r"))
    ReadCollection(contigs).loadFromSam(aln.align(ReadCollection().loadFromFasta(open(query, "r")), contigs)).print_alignments(sys.stdout)



