import random
import sys
import inspect
import traceback
from StringIO import StringIO
from string import ascii_lowercase, ascii_uppercase

from typing import Dict, List, Any

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common.alignment_storage import AlignmentPiece
from common.line_align import Scorer
from common.save_load import TokenReader, TokenWriter
from common.seq_records import NamedSequence
from common.sequences import Contig
from disjointig_resolve.accurate_line import NewLine, NewLineStorage
from disjointig_resolve.correction import Correction
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import DotPlot, LineDotPlot
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage


class Tester:

    def __init__(self, aligner):
        # type: (Aligner) -> None
        self.aligner = aligner
        self.polisher = Polisher(aligner, aligner.dir_distributor)
        testList = []
        for name, obj in inspect.getmembers(sys.modules[__name__]):
            if inspect.isclass(obj) and name.endswith("Test"):
                testList.append(obj)
        self.tests = dict([(c.__name__, c) for c in testList])

    def testAll(self, fname):
        params = self.readParams(fname)
        for name, cases in params.items():
            self.tests[name]().testAll(cases, self.aligner)
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

class TestDataset:
    def __init__(self, genome = "", letter_size = 500, error_rate = 0.05, mutation_rate = 0.005, seed = 0):
        random.seed(seed)
        self.reads = [] # type: List[NamedSequence]
        self.disjointigs = [] # type: List[NamedSequence]
        self.contigs = [] # type: List[NamedSequence]
        self.letter_size = 500
        self.error_rate = error_rate
        self.mutation_rate = mutation_rate
        self.alphabet = dict()
        for c1, c2 in zip(ascii_lowercase, ascii_uppercase):
            seq = self.generate(letter_size)
            self.alphabet[c1] = seq
            self.alphabet[c2] = self.mutate(seq, self.mutation_rate)
        self.genome = Contig(self.translate(genome), genome)

    def translate(self, seq):
        return "".join(map(lambda c: self.alphabet[c], seq))

    def addRead(self, read):
        self.reads.append(NamedSequence(self.mutate(self.translate(read), self.error_rate), "R" + str(len(self.reads)) + "_" + read))

    def addDisjointig(self, disjointig):
        self.disjointigs.append(NamedSequence(self.mutate(self.translate(disjointig), self.mutation_rate), "D" + str(len(self.disjointigs)) + "_" + disjointig))

    def addContig(self, contig):
        self.contigs.append(NamedSequence(self.translate(contig), "C" + str(len(self.contigs)) + "_" + contig))

    def generate(self, letter_size):
        # type: (int) -> str
        return "".join([random.choice(["A", "C", "G", "T"]) for i in range(letter_size)])

    def mutate(self, seq, rate):
        # type: (str, float) -> str
        res = []
        for c in seq:
            if random.random() < rate:
                vars = ["A", "C", "G", "T"]
                vars.remove(c)
                res.append(random.choice([random.choice(vars), "", c + c]))
            else:
                res.append(c)
        return "".join(res)

    def saveStructure(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.genome.id)
        handler.writeInt(len(self.reads))
        for read in self.reads:
            handler.writeToken(read.id)
        handler.writeInt(len(self.disjointigs))
        for disjointig in self.disjointigs:
            handler.writeToken(disjointig.id)
        handler.writeInt(len(self.contigs))
        for contig in self.contigs:
            handler.writeToken(contig.id)

    @staticmethod
    def loadStructure(handler):
        # type: (TokenReader) -> TestDataset
        random.seed(0)
        res = TestDataset(handler.readToken())
        for i in range(handler.readInt()):
            res.addRead(handler.readToken())
        for i in range(handler.readInt()):
            res.addDisjointig(handler.readToken())
        for i in range(handler.readInt()):
            res.addContig(handler.readToken())
        return res



class SimpleTest:
    def __init__(self):
        self.aligner = None # type: Aligner
        self.scorer = Scorer()

    def testAll(self, params, aligner):
        # type: (List[List[str]], Aligner) -> bool
        self.aligner = aligner
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
            except AssertionError as e:
                _, _, tb = sys.exc_info()
                fails.append([str(tn), tb, e.message])
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
        assert al1.cigar == "4M1D2M1I4M", str(al1.asMatchingStrings())
        contig1 = Contig("ACATGATCACT", "from")
        contig2 = Contig("ACGTGAAACGT", "to")
        al = AlignmentPiece.Identical(contig1.asSegment(), contig2.asSegment())
        al1 = self.scorer.polyshAlignment(al)
        assert al1.cigar == "6M1I3M1D1M", str(al1.asMatchingStrings())


class AlignmentStorageTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTACGTACGT", "from")
        contig2 = Contig("ACGTACGTACGT", "to")
        al1 = AlignmentPiece.Identical(contig1.segment(0, 4), contig2.segment(0, 4))
        al2 = AlignmentPiece.Identical(contig1.segment(0, 4), contig2.segment(4, 8))
        al3 = AlignmentPiece.Identical(contig1.segment(4, 8), contig2.segment(8, 12))
        storage = AlignmentStorage()
        storage.addAll([al1, al2, al3])
        assert str(storage.items) == "[(from[0:4]->to[0:4]:1.000), (from[0:4]->to[4:12-4]:1.000), (from[4:12-4]->to[8:12-0]:1.000)]"
        assert str(list(storage.rc)) == "[(-from[4:12-4]->-to[0:4]:1.000), (-from[8:12-0]->-to[4:12-4]:1.000), (-from[8:12-0]->-to[8:12-0]:1.000)]"
        assert str(list(storage.calculateCoverage())) == "[(to[0:12-0], 1)]"
        assert str(list(storage.filterByCoverage(0, 1))) == "[]"
        assert str(list(storage.filterByCoverage(1, 2))) == "[to[0:12-0]]"
        assert str(list(storage.filterByCoverage(2))) == "[]"
        storage.addAndMergeRight(al3)
        assert str(list(storage)) == "[(from[0:4]->to[0:4]:1.000), (from[0:4]->to[4:12-4]:1.000), (from[4:12-4]->to[8:12-0]:1.000)]"
        al4 = AlignmentPiece.Identical(contig1.segment(2, 8), contig2.segment(2, 8))
        al5 = AlignmentPiece.Identical(contig1.segment(4, 10), contig2.segment(4, 10))
        storage.addAll([al4, al5])
        assert str(list(storage.calculateCoverage())) == "[(to[0:2], 1), (to[2:4], 2), (to[4:12-4], 3), (to[8:12-2], 2), (to[10:12-0], 1)]"
        assert str(list(storage.filterByCoverage(2,3))) == "[to[2:4], to[8:12-2]]"
        assert str(list(storage.filterByCoverage(2))) == "[to[2:12-2]]"
        assert str(list(storage.getAlignmentsTo(contig2.segment(2, 3)))) == "[(from[0:4]->to[0:4]:1.000), (from[2:12-4]->to[2:12-4]:1.000)]"
        assert str(list(storage.getAlignmentsTo(contig2.segment(2, 6)))) == "[(from[2:12-4]->to[2:12-4]:1.000)]"


class AlignmentCompositionTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTACGTACGT", "c1")
        contig2 = Contig("ACGTAGGTACGT", "c2")
        contig3 = Contig("ACTTACGTACGT", "c3")
        al1 = AlignmentPiece.Identical(contig1.asSegment(), contig2.asSegment())
        al2 = AlignmentPiece.Identical(contig2.asSegment(), contig3.asSegment())
        al3 = al1.compose(al2)
        assert str(al3) == "(c1[0:12-0]->c3[0:12-0]:0.92)"
        assert al3.cigar == "12M"
        al4 = al1.reverse()
        al5 = al4.composeTargetDifference(al2)
        assert str(al5) == "(c1[0:12-0]->c3[0:12-0]:0.92)"
        assert al5.cigar == "12M"


class CorrectionMappingTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTAAAAGGGTACGT", "c1")
        contig2 = Contig("ACGTAAGGGGGTACGT", "c2")
        al = self.scorer.polyshAlignment(AlignmentPiece.Identical(contig1.segment(5, 12), contig2.segment(5, 12)))
        corr = Correction(contig1, contig2, [al])
        assert corr.mapPositionsUp(range(len(contig2))) == [0, 1, 2, 3, 4, 5, 8, 9, 9, 9, 10, 11, 12, 13, 14, 15]
        assert corr.mapPositionsDown(range(len(contig1))) == [0, 1, 2, 3, 4, 5, 6, 6, 6, 9, 10, 11, 12, 13, 14, 15]
        al2 = AlignmentPiece.Identical(contig2.segment(0, 4))
        al3 = AlignmentPiece.Identical(contig2.segment(6, 8))
        al4 = AlignmentPiece.Identical(contig2.segment(6, 16))
        al5 = AlignmentPiece.Identical(contig2.segment(7, 16))
        assert str(corr.composeQueryDifferences([al2, al3, al4, al5])) == "[(c2[0:4]->c1[0:4]:1.000), (c2[6:7]->c1[8:9]:1.000), (c2[6:16-0]->c1[8:16-0]:0.80), (c2[9:16-0]->c1[9:16-0]:1.000)]"


class DotPlotConstructionTest(SimpleTest):
    def testCase(self, instance):
        data = TokenReader(StringIO(" ".join(instance)))
        dataset = TestDataset.loadStructure(data)
        disjointigs = DisjointigCollection()
        for dis in dataset.disjointigs:
            disjointigs.addNew(dis.seq, dis.id)
        dp = DotPlot(disjointigs)
        dp.construct(self.aligner)
        save = StringIO()
        save_handler = TokenWriter(save)
        dp.save(save_handler)
        tmp = save.getvalue()
        save_str = tmp.replace(" ", "").replace("\n", "")
        result = data.readToken()
        if save_str != result:
            for dis in disjointigs:
                print list(dp.allInter(dis.asSegment()))
        assert save_str == result, "\n" + save_str + "\n" + result

class DotPlotModificationTest(SimpleTest):
    def testManual(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line1 = lines.addNew("ACGTAAAAGGGTACGT", "c1")
        line2 = lines.addNew("ACGTAAGGGGGTACGT", "c2")
        al = self.scorer.polyshAlignment(AlignmentPiece.Identical(line1.asSegment(), line2.asSegment()))
        dp = LineDotPlot(lines, self.aligner)
        dp.addAlignment(al)
        print "before"
        dp.printAll(sys.stdout)
        alignment = AlignmentPiece.Identical(Contig("GC", "tmp").asSegment(), line2.segment(0, 2))
        line2.correctSequence([alignment])
        print "after"
        dp.printAll(sys.stdout)
        Enforce start and end matches for all alignments



# LAUNCH TANYA!!!