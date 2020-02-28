import sys

import os

sys.path.append("py")
from common.SimpleGraph import SimpleGraph
from common import basic


g = SimpleGraph()
g.ReadDot(sys.argv[1])
basic.ensure_dir_existance(sys.argv[2])
cnt = 0
oppa = []
for comp in g.Split(1000000):
    if len(comp) < 3:
        if len(g.v[ comp[0]].inc) + len(g.v[comp[0]].out) + len(g.v[comp[-1]].inc) + len(g.v[comp[-1]].out) <= 2:
            pass
        else:
            oppa.extend(comp)
        if len(oppa) > 30:
            comp = list(oppa)
            oppa = []
        else:
            continue
    print cnt, len(comp)
    f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w")
    g.Draw(comp, f)
    f.close()
    cnt += 1
f = open(os.path.join(sys.argv[2], "small.dot"), "w")
g.Draw(oppa, f)
f.close()





