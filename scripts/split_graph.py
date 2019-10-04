import sys

import os

from common.SimpleGraph import SimpleGraph

sys.path.append("py")
from common import basic


g = SimpleGraph()
g.ReadDot(sys.argv[1])
basic.ensure_dir_existance(sys.argv[2])
for cnt, comp in enumerate(g.Split(50)):
    if len(comp) < 3:
        continue
    print len(comp)
    f = open(os.path.join(sys.argv[2], str(cnt) + ".dot"), "w   ")
    g.Draw(comp, f)
    f.close()
    cnt += 1





