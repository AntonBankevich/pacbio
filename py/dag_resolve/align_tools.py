import itertools
import os

from typing import Optional

from common import basic, sam_parser, SeqIO
from dag_resolve import repeat_graph, sequences
from dag_resolve.params import clean
from flye import polysh_job
from flye.alignment import make_alignment


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
        # type: (int, int) -> tuple[int, int]
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
        # type: (int) -> tuple[int, int]
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

class Consensus:
    def __init__(self, seq, cov, full_seq = None):
        # type: (str, list[int]) -> Consensus
        self.seq = seq
        self.cov = cov
        if full_seq is None:
            self.full_seq = seq
        else:
            self.full_seq = full_seq

    def printQuality(self, handler, cov_threshold = 15):
        # type: (file, int) -> None
        for c, a in zip(self.seq, self.cov):
            if a < cov_threshold:
                handler.write(c.lower())
            else:
                handler.write(c.upper())
        handler.write("\n")

    def cut(self, cov_threshold = 15, length = None):
        l = 0
        while l < len(self.seq) and self.cov[l] >= cov_threshold:
            l += 1
        if length is not None:
            l = min(l, length)
        return Consensus(self.seq[:l], self.cov[:l], self.seq)

    def __len__(self):
        return len(self.seq)

class Aligner:
    def __init__(self, dir, threads = 16):
        self.dir = dir
        self.cur_alignment = 0
        self.threads = threads

    def next_dir(self):
        self.cur_alignment += 1
        name = os.path.join(self.dir, str(self.cur_alignment - 1))
        basic.ensure_dir_existance(name)
        return name

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
            print "Alignment reused"
        else:
            make_alignment(contigs_file, [reads_file], self.threads, alignment_dir, "pacbio", alignment_file)
        return sam_parser.Samfile(open(alignment_file, "r"))

    def matchingAlignment(self, seqs, contig):
        # type: (list[str], sequences.Contig) -> list[AlignedSequences]
        collection = sequences.ContigCollection([contig])
        reads = sequences.ReadCollection(collection)
        res = [] # type: list[AlignedSequences]
        for i in range(len(seqs)):
            reads.add(sequences.Read(SeqIO.SeqRecord(seqs[i], str(i))))
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
        # type: (sequences.Read, sequences.Contig) -> AlignedSequences
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

    def polish(self, reads, consensus):
        # type: (sequences.ReadCollection, sequences.Contig) -> str
        dir = self.next_dir()
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
        dir = self.next_dir()
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
                        if alignment.seg_to.left < 100 and len(read) - alignment.seg_from.left > 3500 and alignment.seg_from.left > best_match:
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

    def polishAndAnalyse(self, reads, consensus, polishing_base = None):
        # type: (sequences.ReadCollection, Edge, sequences.Contig) -> Consensus
        if polishing_base is not None:
            seq = sequences.Contig(self.polish(reads, polishing_base), "Noname")
        else:
            seq = sequences.Contig(self.polishNoConsensus(reads, consensus), "Noname")
        res = [0] * (len(seq) + 1)
        for rec in self.align(reads, sequences.ContigCollection([seq])):
            if rec.is_unmapped:
                continue
            res[rec.pos - 1] += 1
            res[rec.pos + rec.alen - 1] -= 1
        for i in range(1, len(res)):
            res[i] += res[i - 1]
        return Consensus(seq.seq, res)

class SquareRecord:
    def __init__(self, val = 0):
        self.sub_score = val
        self.ins_score = val
        self.del_score = val
        self.sub_info = None
        self.del_info = None
        self.ins_info = None

    def updateIns(self, val, info):
        if val < self.ins_score:
            self.ins_score = val
            # self.ins_info = info

    def updateDel(self, val, info):
        if val < self.del_score:
            self.del_score = val
            # self.del_info = info

    def updateSub(self, val, info):
        if val < self.sub_score:
            self.sub_score = val
            # self.sub_info = info

    def best(self):
        return min(self.sub_score, self.ins_score, self.del_score)

    def bestAction(self):
        if self.sub_score < self.ins_score and self.sub_score < self.del_score:
            return "sub"
        if self.ins_score < self.del_score:
            return "ins"
        return "del"

    def getInfo(self, action):
        if action == "sub":
            return self.sub_info
        if action == "ins":
            return self.ins_info
        if action == "del":
            return self.del_info

class AccurateAligner:
    def __init__(self):
        self.ins_score = 10
        self.del_score = 6
        self.sub_score = 10
        self.homo_score = 4
        self.switch_core = 1
        self.center_score = 20
        self.inf = 10000000

    def align(self, query, pattern, match_position = None):
        # type: (str, str, int) -> int
        # assert len(query) <= len(pattern)
        query = query.upper()
        pattern = pattern.upper()
        res = []# type: list[list[SquareRecord]]
        for i in xrange(len(query) + 1):
            res.append([])
            for j in xrange(len(pattern) + 1):
                res[-1].append(SquareRecord(self.inf))
            res[-1][0].ins_score = self.ins_score * i
            res[-1][0].ins_info = (i - 1, 0, "ins")
            res[-1][0].del_score = self.inf
            res[-1][0].sub_score = self.inf
        for j in xrange(len(pattern) + 1):
            res[0][j].ins_score = self.inf
            res[0][j].del_score = self.inf
            res[0][j].sub_score = 0
        for i in xrange(len(query)):
            for j in xrange(len(pattern)):
                #making sure that matching positions match
                if match_position is not None and j == match_position:
                    res[i + 1][j + 1].del_score = self.inf
                    res[i + 1][j + 1].ins_score = self.inf
                    if query[i] == pattern[j]:
                        if i > 0 and j > 0 and query[i - 1] == pattern[j - 1]:
                            res[i + 1][j + 1].updateSub(res[i][j].sub_score, (i, j, "sub"))
                        else:
                            res[i + 1][j + 1].updateSub(res[i][j].best() + self.switch_core, (i, j, res[i][j].bestAction()))
                    continue
                #calculating scores for substitution case
                if query[i] == pattern[j]:
                    if i > 0 and j > 0 and query[i - 1] == pattern[j - 1]:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score, (i, j, "sub"))
                    else:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score + self.switch_core, (i, j, "sub"))
                    res[i + 1][j + 1].updateSub(res[i][j].del_score + self.switch_core, (i, j, "del"))
                    res[i + 1][j + 1].updateSub(res[i][j].ins_score + self.switch_core, (i, j, "ins"))
                else:
                    if i > 0 and j > 0 and query[i - 1] == pattern[j - 1]:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score + self.sub_score + self.switch_core, (i, j, "sub"))
                    else:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score + self.sub_score, (i, j, "sub"))
                    res[i + 1][j + 1].updateSub(res[i][j].del_score + self.sub_score, (i, j, "del"))
                    res[i + 1][j + 1].updateSub(res[i][j].ins_score + self.sub_score, (i, j, "ins"))
                #calculating scores for homopolymer case
                if i > 0 and query[i] == query[i - 1] and pattern[j] == query[i]:
                    res[i + 1][j + 1].updateIns(res[i][j + 1].best() + self.homo_score, (i, j + 1, res[i][j + 1].bestAction()))
                if j > 0 and pattern[j] == pattern[j - 1] and query[i] == pattern[j]:
                    res[i + 1][j + 1].updateDel(res[i + 1][j].best() + self.homo_score, (i + 1, j, res[i + 1][j].bestAction()))
                #calculating scores for insertion case
                if i > 0 and query[i - 1] == pattern[j]:
                    res[i + 1][j + 1].updateIns(res[i][j + 1].sub_score + self.ins_score + self.switch_core, (i, j + 1, "sub"))
                res[i + 1][j + 1].updateIns(res[i][j + 1].best() + self.ins_score, (i, j + 1, res[i][j + 1].bestAction()))
                #calculating scores for deletion case
                if j > 0 and query[i] == pattern[j - 1]:
                    res[i + 1][j + 1].updateDel(res[i + 1][j].sub_score + self.del_score + self.switch_core, (i + 1, j, "sub"))
                res[i + 1][j + 1].updateDel(res[i + 1][j].best() + self.del_score, (i + 1, j, res[i + 1][j].bestAction()))
                #making sure that matching positions match
        best = self.inf
        cur_query = len(query)
        # cur_pattern = None
        # cur_action = None
        for i in xrange(len(pattern) + 1):
            if best > res[cur_query][i].best():
                cur_pattern = i
                cur_action = res[cur_query][i].bestAction()
                best = res[cur_query][i].best()
        # alignment = [[],[]]
        # info_list = []
        # while cur_query > 0:
        #     info_list.append((cur_query, cur_pattern, cur_action))
        #     print cur_query, cur_pattern, cur_action
        #     cur_query, cur_pattern, cur_action = res[cur_query][cur_pattern].getInfo(cur_action)
        # info_list = info_list[::-1]
        # for cur_query, cur_pattern, cur_action in info_list:
        #     if cur_action == "sub":
        #         alignment[0].append(query[cur_query - 1])
        #         alignment[1].append(pattern[cur_pattern - 1])
        #     elif cur_action == "ins":
        #         alignment[0].append(query[cur_query - 1])
        #         alignment[1].append("-")
        #     else:
        #         alignment[0].append("-")
        #         alignment[1].append(pattern[cur_pattern - 1])
        # print "".join(alignment[0])
        # print "".join(alignment[1])

        # print "Accurate alignment:", query, pattern, best
        return min(best, self.inf)








