import sys
import traceback

from typing import Dict, List, Any

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common.alignment_storage import AlignmentPiece
from common.line_align import Scorer
from common.save_load import TokenReader
from common.sequences import Contig
from disjointig_resolve.smart_storage import SegmentStorage


class Tester:

    def __init__(self, aligner):
        # type: (Aligner) -> None
        self.aligner = aligner
        self.polisher = Polisher(aligner, aligner.dir_distributor)
        testList = [SegmentStorageTest, AlignmentPieceTest, AlignmentPolishingTest]
        self.tests = dict([(c.__name__, c) for c in testList])

    def testAll(self, fname):
        params = self.readParams(fname)
        for name, cases in params.items():
            self.tests[name]().testAll(cases)
        # self.testAlignmentStorage()
        # self.testAlignmentComposition()
        # self.testCorrectionMapping()
        # self.testDotPlotConstruction()
        # self.testDotPlotModification()
        # self.testUniqueRegionMarking()
        # self.testReadRecruitment()
        # self.testLineCorrection()
        # self.testLineExtension()
        # self.testLineMerging()
        # self.testLineKnotting()
        # self.testPotentialReadsConstruction()
        # self.testExtensionConsensus()
        # self.testSaves()

    def readParams(self, fname):
        # type: (str) -> Dict[str, List[List[str]]]
        params = dict()  # type: Dict[str, List[List[str]]]
        handler = TokenReader(open(fname, "r"))
        for i in range(handler.readInt()):
            name = handler.readToken()
            assert name is not None
            params[name] = []
            for j in range(handler.readInt()):
                params[name].append(list(handler.readTokens()))
        return params

class SimpleTest:
    def __init__(self):
        self.scorer = Scorer()

    def testAll(self, params):
        # type: (List[List[str]]) -> bool
        print "Starting test", self.__class__.__name__
        fails = []
        try:
            self.testManual()
        except AssertionError as e:
            _, _, tb = sys.exc_info()
            fails.append(("Manual", tb, e.message))
        for tn, instance in enumerate(params):
            try:
                self.testCase(instance)
            except AssertionError:
                _, _, tb = sys.exc_info()
                fails.append([str(tn), tb])
        if len(fails) == 0:
            print "Finished test", self.__class__.__name__ + ": Passed"
            return True
        else:
            print "Finished test", self.__class__.__name__ + ": Failed"
            for tn, tb, message in fails:
                print "Failed test " + tn + ":"
                traceback.print_tb(tb)
                print "Message:", message
            return False

    def testCase(self, instance):
        pass

    def testManual(self):
        pass


class SegmentStorageTest(SimpleTest):
    def testManual(self):
        contig = Contig("ACGT", "test")
        storage = SegmentStorage()
        storage.add(contig.segment(0,1))
        storage.add(contig.segment(1,2))
        storage.add(contig.segment(2,3))
        storage.add(contig.segment(3,4))
        assert str(storage) == "Storage+:[test[0:1], test[1:2], test[2:4-1], test[3:4-0]]", str(storage)
        assert str(storage.rc) == "Storage-:[-test[0:1], -test[1:2], -test[2:4-1], -test[3:4-0]]", str(storage.rc)
        storage.mergeSegments(1)
        assert str(storage) == "Storage+:[test[0:1], test[1:2], test[2:4-1], test[3:4-0]]", str(storage)
        storage.mergeSegments()
        assert str(storage) == "Storage+:[test[0:4-0]]", str(storage)
        assert str(storage.rc) == "Storage-:[-test[0:4-0]]", str(storage.rc)
        contig = Contig("ACGTACGTACGTACGT", "test")
        storage = SegmentStorage()
        storage.add(contig.segment(0,5))
        storage.add(contig.segment(10,15))
        assert storage.find(contig.segment(5, 10)) == contig.segment(0,5), str(storage.find(contig.segment(5, 10)))
        assert storage.find(contig.segment(6, 10)) == contig.segment(10,15), str(storage.find(contig.segment(6, 10)))
        assert storage.find(contig.segment(5, 9)) == contig.segment(0,5), str(storage.find(contig.segment(5, 9)))
        assert storage.find(contig.segment(0, 16)) == contig.segment(0,5), str(storage.find(contig.segment(0, 16)))

class AlignmentPieceTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTACGTA", "from")
        contig2 = Contig("ACTACGTACGTACAT", "to")
        al1 = AlignmentPiece(contig1.asSegment(), contig2.segment(0, 8), "2M1I6M")
        al2 = AlignmentPiece(contig1.segment(0, 8), contig2.segment(7, 15), "8M")
        glued = AlignmentPiece.GlueOverlappingAlignments([al1, al2])
        assert glued.cigar == "2M1I5M8M", str(glued) + " " + glued.cigar
        assert  glued.seg_from.Seq() == "ACGTACGTACGTACGT", str(glued) + " " + glued.cigar
        assert al1.reduce(query=contig1.segment(0, 2)).cigar == "2M"
        assert al1.reduce(query=contig1.segment(0, 3)).cigar == "2M"
        assert al1.reduce(query=contig1.segment(0, 4)).cigar == "2M1I1M"

class AlignmentPolishingTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTTAAACGT", "from")
        contig2 = Contig("ACGTTTAACGT", "to")
        al = AlignmentPiece.Identical(contig1.asSegment(), contig2.asSegment())
        al1 = self.scorer.polyshAlignment(al)
        assert al1.cigar == "5M1D2M1I3M", str(al1.asMatchingStrings())
