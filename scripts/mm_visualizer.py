import sys

f = sys.argv[1]
s = list(open(f, "r").readlines())
s = [tmp.strip() for tmp in s]
tmp = "".join(s)
s = [tmp[:len(tmp) / 2], tmp[len(tmp) / 2:]]
print len(s)
print map(len, s)
print s[0].strip()
for a, b in zip(s[0], s[1]):
    if a == b:
        sys.stdout.write("|")
    else:
        sys.stdout.write("-")
print ""
print s[1].strip()
