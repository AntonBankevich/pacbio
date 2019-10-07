import sys

import os

sys.path.append("py")
from common.SimpleGraph import SimpleGraph
from common import basic


g = SimpleGraph()
g.ReadDot(sys.argv[1])
basic.ensure_dir_existance(sys.argv[2])
cnt = 0
for comp in g.Split(10):
    if len(comp) < 3:
        continue
    print cnt, len(comp)
    f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w   ")
    g.Draw(comp, f)
    f.close()
    cnt += 1





