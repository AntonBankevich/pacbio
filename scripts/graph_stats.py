import os
import sys
sys.path.append("py")
from common.dot_parser import DotParser


def graphStats(dir):
    for fn in os.listdir(dir):
        if fn.endswith(".dot") or fn.endswith(".gv"):
            ucovs = []
            sys.stdout.write(fn + " ")
            fn = os.path.join(dir, fn)
            rep_len = 0
            utmp = 0
            ulen = 0
            for val in DotParser(open(fn, "r")).parse():
                if val[4].unique:
                    ulen += val[4].len
                    utmp += val[4].len * val[4].cov
                    ucovs.append(val[4].cov)
                else:
                    rep_len += val[4].len * val[4].cov
            tmp = utmp / ulen
            for val in DotParser(open(fn, "r")).parse():
                if val[4].unique:
                    pass
                else:
                    print val[4].len, val[4].cov, (tmp), val[4].cov / (tmp), val[4].len * val[4].cov / (tmp)
            sys.stdout.write(" ".join(map(str, [tmp, rep_len / (tmp), ulen])) + "\n")

if __name__ == "__main__":
    dir = sys.argv[1]
    graphStats(dir)
