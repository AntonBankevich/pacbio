class AccurateAligner:
    def __init__(self):
        self.ins_score = 10
        self.del_score = 6
        self.sub_score = 10
        self.homo_score = 4
        self.switch_core = 1
        self.center_score = 20
        self.inf = 10000000

    def align(self, query, pattern, match_position = None):
        # type: (str, str, int) -> int
        # assert len(query) <= len(pattern)
        query = query.upper()
        pattern = pattern.upper()
        res = []# type: list[list[SquareRecord]]
        for i in xrange(len(query) + 1):
            res.append([])
            for j in xrange(len(pattern) + 1):
                res[-1].append(SquareRecord(self.inf))
            res[-1][0].ins_score = self.ins_score * i
            res[-1][0].ins_info = (i - 1, 0, "ins")
            res[-1][0].del_score = self.inf
            res[-1][0].sub_score = self.inf
        for j in xrange(len(pattern) + 1):
            res[0][j].ins_score = self.inf
            res[0][j].del_score = self.inf
            res[0][j].sub_score = 0
        for i in xrange(len(query)):
            for j in xrange(len(pattern)):
                #making sure that matching positions match
                if match_position is not None and j == match_position:
                    res[i + 1][j + 1].del_score = self.inf
                    res[i + 1][j + 1].ins_score = self.inf
                    if query[i] == pattern[j]:
                        if i > 0 and j > 0 and query[i - 1] == pattern[j - 1]:
                            res[i + 1][j + 1].updateSub(res[i][j].sub_score, (i, j, "sub"))
                        else:
                            res[i + 1][j + 1].updateSub(res[i][j].best() + self.switch_core, (i, j, res[i][j].bestAction()))
                    continue
                #calculating scores for substitution case
                if query[i] == pattern[j]:
                    if i > 0 and j > 0 and query[i - 1] == pattern[j - 1]:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score, (i, j, "sub"))
                    else:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score + self.switch_core, (i, j, "sub"))
                    res[i + 1][j + 1].updateSub(res[i][j].del_score + self.switch_core, (i, j, "del"))
                    res[i + 1][j + 1].updateSub(res[i][j].ins_score + self.switch_core, (i, j, "ins"))
                else:
                    if i > 0 and j > 0 and query[i - 1] == pattern[j - 1]:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score + self.sub_score + self.switch_core, (i, j, "sub"))
                    else:
                        res[i + 1][j + 1].updateSub(res[i][j].sub_score + self.sub_score, (i, j, "sub"))
                    res[i + 1][j + 1].updateSub(res[i][j].del_score + self.sub_score, (i, j, "del"))
                    res[i + 1][j + 1].updateSub(res[i][j].ins_score + self.sub_score, (i, j, "ins"))
                #calculating scores for homopolymer case
                if i > 0 and query[i] == query[i - 1] and pattern[j] == query[i]:
                    res[i + 1][j + 1].updateIns(res[i][j + 1].best() + self.homo_score, (i, j + 1, res[i][j + 1].bestAction()))
                if j > 0 and pattern[j] == pattern[j - 1] and query[i] == pattern[j]:
                    res[i + 1][j + 1].updateDel(res[i + 1][j].best() + self.homo_score, (i + 1, j, res[i + 1][j].bestAction()))
                #calculating scores for insertion case
                if i > 0 and query[i - 1] == pattern[j]:
                    res[i + 1][j + 1].updateIns(res[i][j + 1].sub_score + self.ins_score + self.switch_core, (i, j + 1, "sub"))
                res[i + 1][j + 1].updateIns(res[i][j + 1].best() + self.ins_score, (i, j + 1, res[i][j + 1].bestAction()))
                #calculating scores for deletion case
                if j > 0 and query[i] == pattern[j - 1]:
                    res[i + 1][j + 1].updateDel(res[i + 1][j].sub_score + self.del_score + self.switch_core, (i + 1, j, "sub"))
                res[i + 1][j + 1].updateDel(res[i + 1][j].best() + self.del_score, (i + 1, j, res[i + 1][j].bestAction()))
                #making sure that matching positions match
        best = self.inf
        cur_query = len(query)
        # cur_pattern = None
        # cur_action = None
        for i in xrange(len(pattern) + 1):
            if best > res[cur_query][i].best():
                cur_pattern = i
                cur_action = res[cur_query][i].bestAction()
                best = res[cur_query][i].best()
        # alignment = [[],[]]
        # info_list = []
        # while cur_query > 0:
        #     info_list.append((cur_query, cur_pattern, cur_action))
        #     print cur_query, cur_pattern, cur_action
        #     cur_query, cur_pattern, cur_action = res[cur_query][cur_pattern].getInfo(cur_action)
        # info_list = info_list[::-1]
        # for cur_query, cur_pattern, cur_action in info_list:
        #     if cur_action == "sub":
        #         alignment[0].append(query[cur_query - 1])
        #         alignment[1].append(pattern[cur_pattern - 1])
        #     elif cur_action == "ins":
        #         alignment[0].append(query[cur_query - 1])
        #         alignment[1].append("-")
        #     else:
        #         alignment[0].append("-")
        #         alignment[1].append(pattern[cur_pattern - 1])
        # print "".join(alignment[0])
        # print "".join(alignment[1])

        # print "Accurate alignment:", query, pattern, best
        return min(best, self.inf)


class SquareRecord:
    def __init__(self, val = 0):
        self.sub_score = val
        self.ins_score = val
        self.del_score = val
        self.sub_info = None
        self.del_info = None
        self.ins_info = None

    def updateIns(self, val, info):
        if val < self.ins_score:
            self.ins_score = val
            # self.ins_info = info

    def updateDel(self, val, info):
        if val < self.del_score:
            self.del_score = val
            # self.del_info = info

    def updateSub(self, val, info):
        if val < self.sub_score:
            self.sub_score = val
            # self.sub_info = info

    def best(self):
        return min(self.sub_score, self.ins_score, self.del_score)

    def bestAction(self):
        if self.sub_score < self.ins_score and self.sub_score < self.del_score:
            return "sub"
        if self.ins_score < self.del_score:
            return "ins"
        return "del"

    def getInfo(self, action):
        if action == "sub":
            return self.sub_info
        if action == "ins":
            return self.ins_info
        if action == "del":
            return self.del_info