import sys
sys.path.append("py")
import basic

f = open(sys.argv[1], "r")
s = map(lambda x: x.strip(), f.readlines())
res = [[]]
cur = 0
while cur < len(s):
    

