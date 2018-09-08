import itertools
import os
import sys

from dag_resolve.sequences import Consensus

sys.path.append("py")
from typing import Optional

from common import basic, sam_parser, SeqIO
from dag_resolve import repeat_graph, sequences
from dag_resolve.params import clean
from flye import polysh_job
from flye.alignment import make_alignment

class DirDistributor:
    def __init__(self, dir):
        self.dir = dir
        self.cur_dir = 0

    def nextDir(self):
        name = os.path.join(self.dir, str(self.cur_dir))
        self.cur_dir += 1
        basic.ensure_dir_existance(name)
        return name


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
        # type: (int, int) -> Optional[tuple[int, int]]
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


class Aligner:
    def __init__(self, dir_distributor, threads = 16):
        # type: (DirDistributor, int) -> Aligner
        self.dir_distributor = dir_distributor
        self.cur_alignment = 0
        self.threads = threads

    def next_dir(self):
        return self.dir_distributor.nextDir()

    def CheckSequences(self, reads, reads_file):
        # type: (sequences.ReadCollection, str) -> bool
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
        # type: (sequences.ReadCollection, str) -> bool
        if self.CheckSequences(reads, reads_file):
            return True
        else:
            f = open(reads_file, "w")
            for read in reads:
                SeqIO.write(read, f, "fasta")
            f.close()
            return False

    def align(self, reads, consensus):
        # type: (sequences.ReadCollection, sequences.ContigCollection) -> sam_parser.Samfile
        dir = self.next_dir()
        contigs_file = os.path.join(dir, "contigs.fasta")
        reads_file = os.path.join(dir, "reads.fasta")
        alignment_dir = os.path.join(dir, "alignment")
        alignment_file = os.path.join(dir, "alignment.sam")
        basic.ensure_dir_existance(dir)
        basic.ensure_dir_existance(alignment_dir)
        same_reads = self.CheckAndWriteSequences(reads, reads_file)
        same_contigs = self.CheckAndWriteSequences(consensus, contigs_file)
        if same_reads and same_contigs and not clean and os.path.exists(alignment_file):
            print "Alignment reused:", alignment_file
        else:
            make_alignment(contigs_file, [reads_file], self.threads, alignment_dir, "pacbio", alignment_file)
        return sam_parser.Samfile(open(alignment_file, "r"))

    def matchingAlignment(self, seqs, contig):
        # type: (list[str], sequences.Contig) -> list[AlignedSequences]
        collection = sequences.ContigCollection([contig])
        reads = sequences.ReadCollection(collection)
        res = [] # type: list[AlignedSequences]
        for i in range(len(seqs)):
            reads.add(sequences.AlignedRead(SeqIO.SeqRecord(seqs[i], str(i))))
            res.append(AlignedSequences(seqs[i], contig.seq))
        for rec in self.align(reads, collection):
            tid = int(rec.query_name)
            if res[tid].aligned and res[tid].rc != rec.rc:
                continue
                # assert False, "Straight and reverse alignments of the same read"
            if not res[tid].aligned and rec.rc:
                res[tid].rc = True
                res[tid].seq_from = basic.RC(res[tid].seq_from)
            res[tid].addCigar(rec.cigar, rec.pos)
        return res

    def ReadToAlignedSequences(self, read, contig):
        # type: (sequences.AlignedRead, sequences.Contig) -> AlignedSequences
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

    def repairGraphAlignments(self, graph):
        # type: (repeat_graph.Graph) -> None
        print "Reparing graph alignments for", len(graph.reads), "reads and the following edges:", map(lambda edge: edge.id, graph.newEdges)
        for rec in self.align(graph.reads, graph.newEdges):
            if not rec.is_unmapped:
                read = graph.reads[rec.query_name]
                graph.E[int(rec.tname)].reads.add(read)
                graph.E[int(rec.tname)].reads.addNewAlignment(rec)
        graph.newEdges = []

class Polisher:
    def __init__(self, aligner, dir_distributor):
        # type: (Aligner, DirDistributor) -> Polisher
        self.aligner = aligner
        self.dir_distributor = dir_distributor

    def polish(self, reads, consensus):
        # type: (sequences.ReadCollection, sequences.Contig) -> str
        dir = self.aligner.next_dir()
        consensus_file_name = os.path.join(dir, "ref.fasta")
        consensus_file = open(consensus_file_name, "w")
        SeqIO.write(consensus, consensus_file, "fasta")
        consensus_file.close()
        reads_file_name = os.path.join(dir, "reads.fasta")
        reads_file = open(reads_file_name, "w")
        reads.print_fasta(reads_file)
        reads_file.close()
        res = polysh_job.polish_from_disk(dir, consensus_file_name, reads_file_name)
        return list(SeqIO.parse_fasta(open(res, "r")))[0].seq

    def polishNoConsensus(self, reads, edge):
        # type: (sequences.ReadCollection, repeat_graph.Edge) -> str
        dir = self.dir_distributor.nextDir()
        consensus_file_name = os.path.join(dir, "ref.fasta")
        consensus_file = open(consensus_file_name, "w")
        best_match = 0
        best = ""
        best_id = None
        for read in reads:
            # print read.id
            for alignment in read.alignments:
                # print alignment.__str__()
                if alignment.seg_to.contig.id == edge.id:
                    if not alignment.rc:
                        if alignment.seg_to.left < 100 and len(read) - alignment.seg_from.left > 2000 and alignment.seg_from.left > best_match:
                            best_match = alignment.seg_from.left
                            best = read.seq[alignment.seg_from.left:]
                            best_id = read.id
                    # else:
                    #     if alignment.seg_from.right > len(best):
                    #         best = basic.RC(read.seq[:alignment.seg_from.right])
                    #         best_id = read.id
        consensus = SeqIO.SeqRecord(best, best_id)
        SeqIO.write(consensus, consensus_file, "fasta")
        consensus_file.close()
        reads_file_name = os.path.join(dir, "reads.fasta")
        reads_file = open(reads_file_name, "w")
        reads.print_fasta(reads_file)
        reads_file.close()
        res = polysh_job.polish_from_disk(dir, consensus_file_name, reads_file_name)
        return list(SeqIO.parse_fasta(open(res, "r")))[0].seq

    def polishAndAnalyse(self, reads, polishing_base):
        # type: (sequences.ReadCollection, sequences.Contig) -> Consensus
        seq = sequences.Contig(self.polish(reads, polishing_base), 0)
        res = [0] * (len(seq) + 1)
        for rec in self.aligner.align(reads, sequences.ContigCollection([seq])):
            if rec.is_unmapped:
                continue
            res[rec.pos - 1] += 1
            res[rec.pos + rec.alen - 1] -= 1
        for i in range(1, len(res)):
            res[i] += res[i - 1]
        return Consensus(seq.seq, res)

    def polishQuiver(self, reads, base_start, pos_start, min_new_len = 1000):
        # type: (sequences.ReadCollection, str, int) -> Optional[Consensus]
        cc = sequences.ContigCollection([sequences.Contig(base_start, 0)])
        reads_to_base = sequences.ReadCollection(cc).loadFromSam(self.aligner.align(reads, cc))
        print "Polishing quiver of", len(reads_to_base), "reads."
        for read in sorted(list(reads_to_base), key = lambda read: len(read))[::-1]:
            print read.__str__()
            for al in read.alignments:
                if al.rc:
                    continue
                if al.seg_to.right > len(base_start) - 50 and len(read) - al.seg_from.right > min_new_len:
                    # print al.__str__()
                    tmp = self.polishAndAnalyse(reads, sequences.Contig(base_start[pos_start:al.seg_to.right] + read.seq[al.seg_from.right:], 0))
                    # print len(tmp.cut()), len(base_start), pos_start, min_new_len
                    if len(tmp.cut()) > len(base_start) - pos_start + min_new_len:
                        return tmp
                    break
        return None


if __name__ == "__main__":
    dir = sys.argv[1]
    query = sys.argv[2]
    target = sys.argv[3]
    aln = Aligner(dir)
    contigs = sequences.ContigCollection().loadFromFasta(open(target, "r"))
    sequences.ReadCollection(contigs).loadFromSam(aln.align(sequences.ReadCollection().loadFromFasta(open(query, "r")), contigs)).print_alignments(sys.stdout)



