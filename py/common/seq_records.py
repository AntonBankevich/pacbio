from typing import Optional

from common import basic


class NamedSequence:
    def __init__(self, seq, id):
        # type: (str, str) -> NamedSequence
        self.id = id # type: str
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