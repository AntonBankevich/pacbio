import heapq
import shutil
import sys

import os

from typing import Callable, Union

rc = dict()
rc['A'] = 'T'
rc['C'] = 'G'
rc['G'] = 'C'
rc['T'] = 'A'
rc['a'] = 'T'
rc['c'] = 'G'
rc['g'] = 'C'
rc['t'] = 'A'


def RC(s):
    # type: (str) -> str
    res = []
    for a in s[::-1]:
        # sys.stderr.write(a)
        res.append(rc[a])
    return "".join(res)

def ensure_dir_existance(dir):
    # type: (str) -> None
    if not os.path.exists(dir):
        os.makedirs(dir)

def recreate(dir):
    if os.path.exists(dir):
        shutil.rmtree(dir)
    ensure_dir_existance(dir)

def parseNumber(s, pos=0):
    # type: (str, int) -> Union[None, float, int]
    while pos < len(s) and s[pos] not in "0123456789":
        pos += 1
    if pos == len(s):
        return None
    pos1 = pos
    while pos1 < len(s) and s[pos1] in "0123456789.":
        pos1 += 1
    res = s[pos:pos1]
    if "." in res:
        return float(res)
    else:
        return int(res)

def parseNegativeNumber(s, pos=0):
    # type: (str, int) -> Union[None, float, int]
    while pos < len(s) and s[pos] not in "0123456789":
        pos += 1
    minus = False
    if pos > 0 and s[pos - 1] == '-':
        minus = True
    if pos == len(s):
        return None
    pos1 = pos
    while pos1 < len(s) and s[pos1] in "0123456789.":
        pos1 += 1
    if minus:
        pos -= 1
    res = s[pos:pos1]
    if "." in res:
        return float(res)
    else:
        return int(res)

def best2(arr, better = lambda x, y: x < y):
    # type: (list[int], Callable[[int, int], bool]) -> tuple[int,int]
    assert len(arr) >= 2
    res = [0, 1]
    if better(arr[1], arr[0]):
        res = [1, 0]
    for i, val in enumerate(arr[2:], 2):
        if better(val, arr[res[1]]):
            res[1] = i
            if better(val, arr[res[0]]):
                res = res[::-1]
    return (res[0], res[1])

def merge(*iterables):
    prev = None
    for cur in heapq.merge(*iterables):
        if cur != prev:
            yield cur
            prev = cur

class OStreamWrapper:
    def __init__(self, *streams):
        self.streams = list(streams)

    def write(self, string):
        for stream in self.streams:
            stream.write(string)

    def writelines(self, lines):
        for line in lines:
            for stream in self.streams:
                stream.write(line)

    def flush(self):
        for stream in self.streams:
            stream.flush()

def Link(arr, dist):
    if len(arr) == 0:
        return
    arr = sorted([(pos, i) for i, pos in enumerate(arr)])
    left = arr[0]
    prev = arr[0]
    for val in arr[1:]:
        if val - prev > dist:
            yield left, prev
            left = val
        prev = val
    yield prev, arr[-1]

def Reverse(val):
    if isinstance(val, str):
        assert not val.startswith("--")
        if val.startswith("-"):
            return val[1:]
        else:
            return "-" + val
    elif isinstance(val, int):
        return -val
    assert False, "Tried to reverse an object that is neither number nor a string"