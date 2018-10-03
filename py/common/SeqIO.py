import itertools
import sys
import gzip
from typing import Generator, Union, Callable, Optional

from common import basic


class Reader:
    def __init__(self, handler):
        # type: (file) -> Reader
        self.handler = handler
        self.cash = None # type: str

    def Push(self, value):
        # type: (str) -> None
        assert(self.cash is None)
        self.cash = value

    def FillCash(self):
        if self.cash is None:
            while True:
                self.cash = self.handler.readline()
                if self.cash != "\n":
                    self.cash = self.cash.strip()
                    return

    def TrashCash(self):
        self.cash = None

    def Top(self):
        # type: () -> str
        self.FillCash()
        return self.cash

    def Readline(self):
        # type: () -> str
        self.FillCash()
        result = self.cash
        self.cash = None
        return result

    def ReadUntill(self, f):
        # type: (Callable[[str], bool]) -> str
        result = []
        while True:
            next_line = self.Readline()
            if next_line == "" or f(next_line):
                self.Push(next_line)
                return "".join(result)
            result.append(next_line)

    def ReadUntillFill(self, buf_size):
        # type: (int) -> str
        cnt = 0
        result = []
        while True:
            next_line = self.Readline()
            cnt += len(next_line)
            result.append(next_line)
            if next_line == "" or cnt >= buf_size:
                assert(cnt == buf_size)
                return "".join(result)

    def EOF(self):
        # type: () -> bool
        return self.Top() == ""


class NamedSequence:
    def __init__(self, seq, id):
        # type: (str, Union[str, int]) -> NamedSequence
        self.id = id
        self.seq = seq.upper()

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def RC(self):
        # type: () -> NamedSequence
        return NamedSequence(basic.RC(self.seq), basic.Reverse(self.id))

    def __getitem__(self, item):
        return self.seq[item]


class SeqRecord(NamedSequence):
    def __init__(self, seq, id, qual = None, info = None):
        # type: (str, str, Optional[str], Optional[str]) -> SeqRecord
        assert qual is None or len(qual) == len(seq)
        NamedSequence.__init__(self, seq, id)
        self.qual = qual
        self.info = info

    def RC(self):
        # type: () -> SeqRecord
        qual = self.qual
        if qual is not None:
            qual = qual[::-1]
        return SeqRecord(basic.RC(self.seq), basic.Reverse(self.id), qual)

    def __getitem__(self, key):
        # type: (int) -> str
        return self.seq[key]

    def QualSubseq(self, l, r):
        # type: (int, int) -> Optional[str]
        if self.qual is not None:
            return self.qual[l: r]
        return None

    def subseq(self, l, r):
        # type: (int, int) -> SeqRecord
        if l != 0 or r != len(self.seq):
            return SeqRecord(self.seq[l:r], self.id + "(" + str(l + 1) +"-" + str(r) + ")", self.QualSubseq(l, r))
        else:
            return self


def parse(handler, file_name):
    file_type = file_name.split(".")[-1]
    assert file_type in ["fasta", "fastq", "fa", "fq"]
    if file_type in ["fasta", "fa"]:
        return parse_fasta(handler)
    if file_type in ["fastq", "fq"]:
        return parse_fastq(handler)


def parse_fasta(handler):
    # type: (file) -> Generator[SeqRecord]
    reader = Reader(handler)
    while not reader.EOF():
        rec_id = reader.Readline().strip()
        info = None
        pos = rec_id.find(" ")
        if pos != -1:
            info = rec_id[pos + 1:]
            rec_id = rec_id[:pos]
        assert(rec_id[0] == '>')
        rec_seq = reader.ReadUntill(lambda s: s.startswith(">"))
        yield SeqRecord(rec_seq, rec_id[1:], None, info)


def parse_fastq(handler):
    reader = Reader(handler)
    while not reader.EOF():
        rec_id = reader.Readline()
        assert(rec_id[0] == '@')
        rec_seq = reader.ReadUntill(lambda s: s.startswith("+"))
        tmp = reader.Readline()
        assert(tmp[0] == '+')
        rec_qual = reader.ReadUntillFill(len(rec_seq))
        yield SeqRecord(rec_seq, rec_id[1:], rec_qual)


def write(rec, handler, file_type):
    if file_type == "fasta":
        handler.write(">" + str(rec.id) + "\n")
        handler.write(rec.seq + "\n")
    elif file_type == "fastq":
        handler.write("@" + rec.id + "\n")
        handler.write(rec.seq + "\n")
        handler.write("+" + "\n")
        handler.write(rec.qual + "\n")


def FilterContigs(input_handler, output_handler, f, file_type):
    for contig in parse(input_handler, file_type):
        if f(contig):
            write(contig, output_handler, file_type)


def RemoveNs(input_handler, output_handler):
    for contig in parse(input_handler, "fasta"):
        l = 0
        while l < len(contig) and contig[l] == 'N':
            l += 1
        r = len(contig)
        while r > l and contig[r - 1] == 'N':
            r -= 1
        if r > l:
            output_handler.write(SeqRecord(contig.seq[l:r], contig.id))
