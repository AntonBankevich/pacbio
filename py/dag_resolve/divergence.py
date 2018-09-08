class DivergenceState:
    def __init__(self, div, seq, id):
        # type: (Divergence, str) -> DivergenceState
        self.divergence = div
        self.seq = seq
        self.id = id

    def isAmbiguous(self):
        # type: () -> bool
        return self.seq is None

    def char(self):
        if self.isAmbiguous():
            return "-"
        return str(self.id)

    def __eq__(self, other):
        # type: (DivergenceState) -> bool
        return self.seq == other.seq and self.divergence == other.divergence

    def __ne__(self, other):
        return not self.__eq__(other)


class Divergence:
    def __init__(self, edge, pos, states = []):
        # type: (repeat_graph.Edge, tuple[int, int], list[str]) -> Divergence
        self.edge = edge
        self.pos = pos
        self.states = [] # type: list[DivergenceState]
        self.ambiguous = DivergenceState(self, None, -1)
        for state in states:
            self.addState(state)
        self.statistics = Statistics()

    def addState(self, state):
        # type: (str) -> DivergenceState
        self.states.append(DivergenceState(self, state, len(self.states)))
        return self.states[-1]

    def __eq__(self, other):
        # type: (Divergence) -> bool
        return self.edge.id == other.edge.id and self.pos == other.pos

    def __ne__(self, other):
        return not self.__eq__(other)

    def printToFile(self, handler, delim = " "):
        # type: (file) -> None
        handler.write(str(self.edge.id) + ": " + str(self.pos) + "\n")
        for state in self.states:
            handler.write(state.seq + " ")
        handler.write("\n")


class Statistics:
    def __init__(self):
        self.correct = 0 # type: int
        self.wrong = 0 # type: int
        self.ambig = 0 # type: int
        self.called = 0 # type: int

    def __str__(self):
        return str((self.correct, self.wrong, self.ambig, self.called))