import sys


def pathStats(fn, k):
    f = open(fn, "r")
    arr = map(int, f.readlines())
    res = [0] * 1000
    for i in range(len(arr)):
        s = 0
        j = i
        while s <= k and j < len(arr):
            s += arr[j]
            j += 1
        res[min(len(res) - 1, j - i)] += 1
    for x in res:
        print x




if __name__ == "__main__":
    pathStats(sys.argv[1], int(sys.argv[2]))