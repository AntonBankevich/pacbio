import sys

import os

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

def parseNumber(s, pos=0):
    # type: (str, int) -> int
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
    # type: (str, int) -> int
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

def diff(s1, s2):
    # type: (str, str) -> int
    assert False, "Diff not implemented"