import sys

import os

from typing import Dict

from dag_resolve.repeat_graph import Graph, DotParser
from dag_resolve.sequences import ContigCollection


def analyse(name, myprog, trestile):
    # type: (str, str) -> Dict[str, int]
    log = open(myprog, "r").readlines()[-30:]
    res = dict()
    for s in log:
        s = s.strip()
        if s.startswith("Stat: "):
            s = s[6:].split(":")
            res[s[0]] = int(s[1])
    tres = open(trestile, "r").readlines()[1:]
    tr = 0
    for s in tres:
        if s.find("True") != -1:
            tr += 1
    res["Trestile"] = tr
    res["Name"] = name
    return res


if __name__ == "__main__":
    indir = sys.argv[1]
    outdir = sys.argv[2]
    col = ["Name", "Unique edges", "Repeat edges", "Repeat components", "2in2out", "Resolved edges", "Unique edge connections", "Trestile"]
    for s in col:
        sys.stdout.write(s + "\t")
    for dir in os.listdir(indir):
        tres = os.path.join(indir, dir, "4-trestle", "trestle_summary.txt")
        out = os.path.join(outdir, dir, "log.info")
        tmp = analyse(dir, out, tres)
        for s in col:
            sys.stdout.write(tmp[s] + "\t")
        sys.stdout.write("\n")
