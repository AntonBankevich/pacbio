from common import sam_parser, SeqIO


class ContigCollection():
    def __init__(self):
        self.contigs = dict()
        self.main = None

    def add(self, contig):
        if "main" in contig.info:
            self.main = contig
        self.contigs[contig.id] = contig

    def filter(self, condition):
        res = ContigCollection()
        for contig in self.contigs.values():
            if condition(contig):
                res.add(contig)
        return res

    def incoming(self):
        return self.filter(lambda contig: "in" in contig.info)

    def outgoing(self):
        return self.filter(lambda contig: "out" in contig.info)

    def main(self):
        return self.main

    def print_names(self, handler):
        for contig in self:
            handler.write(contig.id + " " + " ".join(contig.info) + "\n")

    def __iter__(self):
        return self.contigs.values().__iter__()

    def __getitem__(self, contig_id):
        return self.contigs[contig_id]

class Contig:
    def __init__(self, seq, id, info):
        self.seq = seq
        self.id = id
        self.info = info

    def end(self, len):
        return Segment(self, self.__len__() - len, self.__len__())

    def start(self, len):
        return Segment(self, 0, len)

    def as_segment(self):
        return Segment(self, 0, len(self))

    def print_fasta(self, handler):
        SeqIO.write(self, handler, "fasta")

    def __len__(self):
        return len(self.seq)

class Segment:
    def __init__(self, contig, left, right):
        self.contig = contig
        self.left = left
        self.right = right

    def inter(self, other):
        return self.contig.id == other.contig.id and not (self.right < other.left or self.left > other.right)

    def contains(self, other):
        return self.contig.id == other.contig.id and self.left <= other.left and self.right >= other.right

    def subcontig(self):
        return Contig(self.contig.seq[self.left:self.right], self.contig.id + "[" + str(self.left) + "," + str(self.right) + "]", self.contig.info)

    def __str__(self):
        return self.contig.id + "[" + str(self.left) + ":" + str(self.right) + "]"


class AlignmentPiece:
    def __init__(self, seg_from, seg_to, cigar = None):
        self.seg_from = seg_from
        self.seg_to = seg_to
        self.cigar = cigar

    def __str__(self):
        return "(" + str(self.seg_from) + "->" + str(self.seg_to) + ")"

class Read:
    def __init__(self, rec):
        self.id = rec.id
        self.seq = rec.seq
        self.alignments = []

    def __len__(self):
        return len(self.seq)

    def AddSamAlignment(self, rec, contig):
        cigar_list = list(sam_parser.CigarToList(rec.cigar))
        ls = 0
        rs = 0
        if cigar_list[0][0] in "HS":
            ls = cigar_list[0][1]
        if cigar_list[-1][0] in "HS":
            rs = cigar_list[-1][1]
        self.alignments.append(AlignmentPiece(Segment(self, ls, len(self.seq) - rs), Segment(contig, rec.pos, rec.pos + rec.alen), rec.cigar))

    def inter(self, other):
        for ap in self.alignments:
            if ap.seg_to.inter(other):
                return True
        return False

    def contains(self, other):
        for ap in self.alignments:
            if ap.seg_to.contains(other):
                return True
        return False

    def contains_start(self, contig):
        return self.inter(Segment(contig, 0, 200))

    def contains_end(self, contig):
        return self.inter(Segment(contig, len(contig) - 200, len(contig)))

class ReadCollection:
    def __init__(self, contigs):
        self.reads = dict()
        self.contigs = contigs

    def addNewRead(self, rec):
        rec.id = rec.id.split()[0]
        self.reads[rec.id] = Read(rec)

    def addNewAlignment(self, rec):
        if rec.is_unmapped:
            return None
        rname = rec.query_name.split()[0]
        if rname not in self.reads:
            return None
        self.reads[rname].AddSamAlignment(rec, self.contigs[rec.tname])

    def add(self, read):
        self.reads[read.id] = read

    def filter(self, condition):
        res = ReadCollection(self.contigs)
        for read in self.reads.values():
            if condition(read):
                res.add(read)
        return res

    def inter(self, segment):
        return self.filter(lambda read: read.inter(segment))

    def contain(self, segment):
        return self.filter(lambda read: read.contain(segment))

    def __iter__(self):
        return self.reads.values().__iter__()

    def __getitem__(self, read_id):
        return self.reads[read_id]

    def print_fasta(self, hander):
        for read in self:
            SeqIO.write(read, hander, "fasta")

    def print_alignments(self, handler):
        for read in self:
            handler.write(read.id + "\n")
            for rec in read.alignments:
                handler.write(str(rec) + "\n")
