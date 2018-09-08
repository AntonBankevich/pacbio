class Phasing:
    def __init__(self, states = None):
        # type: (list[DivergenceState]) -> Phasing
        if states is None:
            states = []
        self.states = states

    def add(self, state):
        # type: (DivergenceState) -> DivergenceState
        self.states.append(state)
        return self.states[-1]

    def printToFile(self, handler):
        # type: (file) -> None
        handler.write(str(self))
        handler.write("\n")

    def __str__(self):
        return "".join(map(lambda state:state.char(), self.states))

    def ambibuousRate(self):
        res = 0
        for state in self.states:
            if state.isAmbiguous():
                res += 1
        return float(state) / len(self.states)

    def called(self):
        res = 0
        for state in self.states:
            if not state.isAmbiguous():
                res += 1
        return res

    def matching(self, other):
        # type: (Phasing) -> int
        res = 0
        for ph1, ph2 in zip(self.__iter__(), other.__iter__()):
            assert ph1.divergence == ph2.divergence
            if not ph1.isAmbiguous() and ph1 == ph2:
                res += 1
        return res

    def diff(self, other):
        # type: (Phasing) -> int
        res = 0
        for ph1, ph2 in zip(self.__iter__(), other.__iter__()):
            assert ph1.divergence == ph2.divergence
            if not ph1.isAmbiguous() and not ph2.isAmbiguous() and ph1 != ph2:
                res += 1
        return res

    def diff2(self, other1, other2):
        # type: (Phasing, Phasing) -> tuple[int,int, int]
        d1, d2, diff = 0, 0, 0
        for ph, ph1, ph2 in zip(self.__iter__(), other1.__iter__(), other2.__iter__()):
            if not ph.isAmbiguous() and ph1 != ph2:
                diff += 1
                if ph != ph1:
                    d1 += 1
                if ph != ph2:
                    d2 += 1
        return (d1, d2, diff)

    def __getitem__(self, item):
        # type: (int) -> DivergenceState
        return self.states[item]

    def __iter__(self):
        # type: () -> Iterator[DivergenceState]
        return self.states.__iter__()

    def __len__(self):
        # type: () -> int
        return self.states.__len__()