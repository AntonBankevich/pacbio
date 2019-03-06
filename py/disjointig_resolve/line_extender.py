from typing import Optional

from common.sequences import Segment
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.disjointigs import DisjointigCollection

k = 1000

class LineExtender:
    def __init__(self, disjointigs):
        # type: (DisjointigCollection) -> None
        self.disjointigs = disjointigs

    def tryExtend(self, line):
        # type: (NewLine) -> bool
        line.completely_resolved.mergeSegments(k)
        result = []
        for i in range(len(line.completely_resolved) - 1):
            seg = line.completely_resolved[i]
            new_segment = self.attemptCleanResolution(seg, line.completely_resolved[i + 1].left + k)
            if new_segment != seg:
                result.append(new_segment)
        seg = line.completely_resolved[-1]
        new_segment = self.attemptCleanResolution(seg, len(line))
        if new_segment != seg:
            result.append(new_segment)
        if len(result) > 0:
            line.completely_resolved.addAll(result)
            line.completely_resolved.mergeSegments(k)
            return True
        else:
            return False

    def attemptCleanResolution(self, resolved, bound):
        # type: (Segment, int) -> Optional[Segment]
        # IMPLEMENT!
        pass