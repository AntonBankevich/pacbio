import datetime
import os

from typing import BinaryIO, Iterable, Generator, Optional

from common import basic


class TokenWriter:
    def __init__(self, handler):
        # type: (BinaryIO) -> None
        self.handler = handler
        self.new_line = True

    def writeToken(self, token):
        # type: (str) -> None
        assert token.find(" ") == -1 and token.find("\n") == -1 and token != "", "Token:" + token
        if not self.new_line:
            self.handler.write(" ")
        self.new_line = False
        self.handler.write(token)

    def writeInt(self, val):
        # type: (str) -> None
        self.handler.write(str(val))

    def writeTokenLine(self, token):
        # type: (str) -> None
        if not self.new_line:
            self.handler.write(" ")
        self.new_line = False
        self.handler.write(token)
        self.newLine()

    def writeIntLine(self, val):
        # type: (int) -> None
        self.writeToken(str(val))

    def newLine(self):
        self.handler.write("\n")
        self.new_line = True

    def writeTokens(self, tokens):
        # type: (Iterable[str]) -> None
        self.writeToken("List")
        for token in tokens:
            self.writeToken(token)
        self.newLine()

class TokenReader:
    def __init__(self, handler):
        # type: (BinaryIO) -> None
        self.handler = handler
        self.line = None
        self.pos = None

    def readToken(self):
        # type: () -> Optional[str]
        if self.line is None or self.pos == len(self.line):
            line = self.handler.readline().split()
            pos = 0
        self.pos += 1
        if self.line[self.pos - 1] == "None":
            return None
        else:
            return self.line[self.pos - 1]

    def readInt(self):
        # type: () -> int
        return int(self.readToken())

    def readTokens(self):
        # type: () -> Generator[str]
        check = self.readToken()
        assert check == "List"
        for token in self.line[self.pos:]:
            yield token
        self.pos = len(self.line)

class SaveHandler:
    def __init__(self, dir, clean = False):
        # type: (str, bool) -> None
        self.dir = dir
        if clean:
            basic.recreate(self.dir)
        else:
            basic.ensure_dir_existance(self.dir)
        self.cnt = 0
        for name in os.listdir(self.dir):
            num = basic.parseNumber(name)
            if num is not None and num >= self.cnt:
                self.cnt = num + 1

    def getWriter(self, suffix = None):
        name = str(self.cnt) + "_" + datetime.datetime.now().strftime('%Y:%m:%d:%H:%M')
        if suffix is not None:
            name += "_" + suffix
        return TokenWriter(open(os.path.join(self.dir, name), "w"))