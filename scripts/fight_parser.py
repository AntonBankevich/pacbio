import sys

def main(args):
    res = dict()
    cnt = 0
    win = dict()
    for s in open(args[1], "r"):
        if not s.startswith("Fight"):
            continue
        cnt += 1
        cur = 0
        l = s.find("((", cur) + 2
        r = s.find("[", cur)
        rn = s[l:r]
        cur = r
        l = s.find("->", cur) + 2
        cur = l
        r = s.find("[", cur)
        c1 = s[l:r]
        cur = r
        l = s.find("->", cur) + 2
        cur = l
        r = s.find("[", cur)
        c2 = s[l:r]
        cur = r
        cur = s.find("Comparison results: ", cur) + len("Comparison results: ")
        vals = s[cur:].split(" ")[:4]
        if c1 < c2:
            key = (rn,c1,c2)
        else:
            key = (rn, c2, c1)
        if s.find("Winner") != -1:
            cur = s.find("Winner")
            l = s.find("->", cur) + 2
            cur = l
            r = s.find("[", cur)
            w = s[l:r]
            win[key] = w
        if vals[0] == "None" or vals[0] == "1000":
            continue
        if key not in res:
            res[key] = []
        res[key].append(vals)
        # if cnt % 1000 == 0:
        #     for a, b in res.items():
        #         print a
        #         print b
    for a, b in res.items():
        tmp = filter(lambda s: str(sorted(s)) != str(sorted(b[-1])), b)
        if len(tmp) == 0 or a not in win:
            continue
        if a in win:
            print a, win[a]
        else:
            continue
            print a
        if len(tmp) > 0:
            print tmp
            print b

if __name__ == "__main__":
    main(sys.argv)