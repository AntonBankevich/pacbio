import itertools
import os
import sys

from dag_resolve import params
from dag_resolve.repeat_graph import Graph
from dag_resolve.sequences import AlignedRead, Contig, ContigCollection, ReadCollection, AlignmentPiece

sys.path.append("py")
from typing import Optional, Iterable, Tuple

from common import basic, sam_parser, SeqIO
from flye.alignment import make_alignment

class DirDistributor:
    def __init__(self, dir):
        self.dir = dir
        self.cur_dir = 0

    def nextDir(self):
        name = os.path.join(self.dir, str(self.cur_dir))
        self.cur_dir += 1
        assert self.cur_dir <= 2000
        basic.ensure_dir_existance(name)
        return name

    def CheckSequences(self, reads, reads_file):
        # type: (Iterable[AlignedRead], str) -> bool
        if not os.path.exists(reads_file):
            return False
        try:
            for rec, read in itertools.izip_longest(SeqIO.parse_fasta(open(reads_file, "r")), reads):
                if str(rec.id) != str(read.id) or rec.seq != read.seq:
                    return False
            return True
        except:
            return False

    def CheckAndWriteSequences(self, reads, reads_file):
        # type: (Iterable[AlignedRead], str) -> bool
        if self.CheckSequences(reads, reads_file):
            return True
        else:
            f = open(reads_file, "w")
            for read in reads:
                SeqIO.write(read, f, "fasta")
            f.close()
            return False

    def fillNextDir(self, content):
        # type: (list[Tuple[Iterable[AlignedRead], str]]) -> Tuple[str, list[str], bool]
        same = True
        dir = self.nextDir()
        content_files = []
        for reads, f_name in content:
            f_name = os.path.join(dir, f_name)
            if not self.CheckAndWriteSequences(reads, f_name):
                same = False
            content_files.append(f_name)
        return dir, content_files, same


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
        # type: (int) -> Optional[tuple[Optional[int], Optional[int]]]
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
        if rec.seg_to.contig.id != contig.id:
            continue
        if res.aligned and res.rc != rec.rc:
            continue
            # assert False, "Straight and reverse alignments of the same read"
        if not res.aligned and rec.rc:
            res.rc = True
            res.seq_from = basic.RC(res.seq_from)
        res.addCigar(rec.cigar, rec.seg_to.left)
    return res


class Aligner:
    def __init__(self, dir_distributor, threads = 16):
        # type: (DirDistributor, int) -> Aligner
        self.dir_distributor = dir_distributor
        self.cur_alignment = 0
        self.threads = threads

    def align(self, reads, reference):
        # type: (Iterable[AlignedRead], Iterable[Contig]) -> sam_parser.Samfile
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

    def matchingAlignment(self, seqs, contig):
        # type: (list[str], Contig) -> list[AlignedSequences]
        collection = ContigCollection([contig])
        reads = ReadCollection(collection)
        res = [] # type: list[AlignedSequences]
        for i in range(len(seqs)):
            reads.add(AlignedRead(SeqIO.SeqRecord(seqs[i], str(i))))
            res.append(AlignedSequences(seqs[i], contig.seq))
        reads.loadFromSam(self.align(reads, collection))
        for read in reads:
            print read.__str__()
            tid = int(read.id)
            read.sort()
            groups = [] #type: list[list[AlignmentPiece]]
            group_lens = []
            for al in read.alignments:
                if al.rc:
                    continue
                found = False
                for i, group in enumerate(groups):
                    if group[-1].connect(al):
                        group.append(al)
                        group_lens[i] += len(al.seg_from)
                        found = True
                        break
                if not found:
                    groups.append([al])
                    group_lens.append(len(al.seg_from))
            best = None
            for i in range(len(groups)):
                if best == None or group_lens[i] > group_lens[best]:
                    best = i
            for al in groups[best]:
                res[tid].addCigar(al.cigar, al.seg_to.left)
            res[tid].rc = groups[best][0].rc
            if res[tid].rc:
                res[tid].seq_from = basic.RC(res[tid].seq_from)
        return res

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



