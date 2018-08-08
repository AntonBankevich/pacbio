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
    res = []
    for a in s[::-1]:
        # sys.stderr.write(a)
        res.append(rc[a])
    return "".join(res)

def ensure_dir_existance(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)