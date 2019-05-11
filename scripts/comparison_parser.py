import sys
sys.path.append("py")

f = open(sys.argv[1], "r")
s = map(lambda x: x.strip(), f.readlines())
res = [[]]
cur = 0
while cur < len(s):
    if s[cur].find("Comparison results") != -1:
        res[-1].append(s[cur])
    elif s[cur].find("hampion") != -1:
        res[-1].append(s[cur])
        res.append([])
    cur += 1
res = res[:-1]
fights = []
for tmp in res:
    if tmp[-1].find("False") != -1:
        continue
    if tmp[-1].find("Champion") == -1:
        champ = None
    else:
        champ = tmp[-1][tmp[-1].find("NCTC"):tmp[-1].find("[")]
    for f in tmp[:-1]:
        x = (champ is not None and f.find(champ) != -1)
        f = f[f.find("Comparison results"):].split()[2:]
        if "None" in f:
            continue
        fights.append([x, f[0], f[1], f[2], f[3]])
for f in fights:
    print " ".join(map(str, f))

