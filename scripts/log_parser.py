import sys
sys.path.append("py")
from common import basic

ll = []
start = False
for s in open(sys.argv[1], "r").readlines():
    if s.startswith("New depth 24218"):
        start = True
    if not start:
        continue
    if s.startswith("Tail"):
        continue
    if s.startswith("New tails"):
        break
    if s.find("Fail") != -1:
        ll = []
    elif s.find("success") != -1:
        call = basic.parseNumber(s)
        for tmp in map(eval, ll):
            if tmp[call] < tmp[1 - call]:
                print min(tmp), max(tmp), 0
            else:
                print min(tmp), max(tmp), 1
        ll = []
    elif s.find("[") != -1:
        ll.append(s.strip())

