import os

from common import basic


def parseUPaths(flye_dir):
    id_map = dict()
    for s in open(os.path.join(flye_dir, "flye.log"), "r").readlines():
        if s.find("UPath") != -1:
            s = s.split()
            id_map[s[4][:-1]] = s[5::2]
            id_map["-" + s[4][:-1]] = map(lambda val: basic.Reverse(val), s[5::2])
    return id_map

