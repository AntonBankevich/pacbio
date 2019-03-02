from typing import List

from common.save_load import TokenWriter, TokenReader


class Params:
    def __init__(self):
        self.reads_file = None
        self.contigs_file = None
        self.disjointigs_file = None
        self.load_from = None
        self.dir = None
        self.save_dir = None

    def parse(self, args):
        # type: (List[str]) -> Params
        return self

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.reads_file)
        handler.writeTokenLine(self.contigs_file)
        handler.writeTokenLine(self.disjointigs_file)
        handler.writeTokenLine(self.load_from)
        handler.writeTokenLine(self.dir)
        handler.writeTokenLine(self.save_dir)

    def load(self, handler):
        # type: (TokenReader) -> None
        self.reads_file = handler.readToken()
        self.contigs_file = handler.readToken()
        self.disjointigs_file = handler.readToken()
        self.load_from = handler.readToken()
        self.dir = handler.readToken()
        self.save_dir = handler.readToken()
