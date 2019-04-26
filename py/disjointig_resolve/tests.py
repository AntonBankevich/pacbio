import random
import sys
import inspect
import traceback
from StringIO import StringIO
from string import ascii_lowercase, ascii_uppercase

from typing import Dict, List, Any, Tuple

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import params
from common.alignment_storage import AlignmentPiece, ReadCollection
from common.line_align import Scorer
from common.save_load import TokenReader, TokenWriter
from common.seq_records import NamedSequence
from common.sequences import Contig
from disjointig_resolve.accurate_line import NewLine, NewLineStorage
from disjointig_resolve.correction import Correction
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import DotPlot, LineDotPlot
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage
from disjointig_resolve.unique_marker import UniqueMarker


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
        fail = 0
        for name, cases in params.items():
            result = self.tests[name]().testAll(cases, self.aligner)
            if not result:
                fail += 1
        print fail, "tests failed"
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
    def __init__(self, genome = "", letter_size = 550, error_rate = 0.05, mutation_rate = 0.005, seed = 0):
        random.seed(seed)
        self.reads = [] # type: List[NamedSequence]
        self.disjointigs = [] # type: List[NamedSequence]
        self.contigs = [] # type: List[NamedSequence]
        self.letter_size = letter_size
        self.error_rate = error_rate
        self.mutation_rate = mutation_rate
        self.alphabet = dict()
        self.matches = dict()
        for c1, c2 in zip(ascii_lowercase, ascii_uppercase):
            seq = self.generate(self.letter_size)
            self.alphabet[c1] = seq
            seq, matches = self.mutate(seq, self.mutation_rate)
            self.alphabet[c2] = seq
            self.matches[c1] = matches
            self.matches[c2] = [(b, a) for a, b in matches]
        self.genome = Contig(self.translate(genome), genome)

    def translate(self, seq):
        return "".join(map(lambda c: self.alphabet[c], seq))

    def addRead(self, read_seq):
        self.reads.append(NamedSequence(self.mutate(self.translate(read_seq), self.error_rate)[0], "R" + str(len(self.reads)) + "_" + read_seq))

    def addDisjointig(self, disjointig_seq):
        # type: (str) -> str
        self.disjointigs.append(NamedSequence(self.mutate(self.translate(disjointig_seq), self.mutation_rate)[0], "D" + str(len(self.disjointigs)) + "_" + disjointig_seq))
        return self.disjointigs[-1].id

    def addContig(self, contig_seq):
        # type: (str) -> str
        name = "C" + str(len(self.contigs)) + "_" + contig_seq
        self.contigs.append(NamedSequence(self.translate(contig_seq), name))
        return name

    def generateReads(self, length = 5, cov = 15, circular = False):
        genome = self.genome.id
        if circular:
            genome = genome + genome[0:length - 1]
        for i in range(0, len(genome) - length + 1):
            for j in range((cov + length - 1) / length):
                self.addRead(genome[i:i + length])


    def generate(self, letter_size):
        # type: (int) -> str
        return "".join([random.choice(["A", "C", "G", "T"]) for i in range(letter_size)])

    def genAll(self, aligner):
        # type: (Aligner) -> Tuple[NewLineStorage, LineDotPlot, ReadCollection]
        disjointigs = DisjointigCollection()
        for dis in self.disjointigs:
            disjointigs.addNew(dis.seq, dis.id)
        lines = NewLineStorage(disjointigs, aligner)
        for line in self.contigs:
            lines.addNew(line.seq, line.id)
        dp = LineDotPlot(lines, aligner)
        dp.construct(aligner)
        lines.alignDisjointigs()
        reads = ReadCollection()
        for read in self.reads:
            reads.addNewRead(read)
        disjointigs.addAlignments(aligner.alignClean(reads, disjointigs))
        return lines, dp, reads

    def mutate(self, seq, rate):
        # type: (str, float) -> Tuple[str, List[Tuple[int, int]]]
        res = [seq[0]]
        matches = []
        matches.append((0,0))
        cur = 1
        for i, c in enumerate(seq):
            if i == 0 or i == len(seq) - 1:
                continue
            if random.random() < rate:
                vars = ["A", "C", "G", "T"]
                vars.remove(c)
                res.append(random.choice([random.choice(vars), "", c + c]))
                cur += len(res[-1])
            else:
                res.append(c)
                matches.append((cur, i))
                cur += 1
        res.append(seq[-1])
        matches.append((len(seq) - 1, cur))
        return "".join(res), matches

    def saveStructure(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.genome.id)
        handler.writeInt(len(self.reads))
        for read in self.reads:
            handler.writeToken(read.id.split("_")[-1])
        handler.writeInt(len(self.disjointigs))
        for disjointig in self.disjointigs:
            handler.writeToken(disjointig.id.split("_")[-1])
        handler.writeInt(len(self.contigs))
        for contig in self.contigs:
            handler.writeToken(contig.id.split("_")[-1])

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

    def assertResult(self, res, ethalon):
        # type: (str, str) -> None
        assert res.replace(" ", "") == ethalon, res.replace(" ", "")
        pass

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
        assert str(storage) == "ReadStorage+:[test[0:1], test[1:2], test[2:4-1], test[3:4-0]]", str(storage)
        assert str(storage.rc) == "ReadStorage-:[-test[0:1], -test[1:2], -test[2:4-1], -test[3:4-0]]", str(storage.rc)
        storage.mergeSegments(1)
        assert str(storage) == "ReadStorage+:[test[0:1], test[1:2], test[2:4-1], test[3:4-0]]", str(storage)
        storage.mergeSegments()
        assert str(storage) == "ReadStorage+:[test[0:4-0]]", str(storage)
        assert str(storage.rc) == "ReadStorage-:[-test[0:4-0]]", str(storage.rc)
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
        test_result = tmp.replace(" ", "").replace("\n", "")
        ethalon = data.readToken()
        if test_result != ethalon:
            for dis in disjointigs:
                print list(dp.allInter(dis.asSegment()))
        assert test_result == ethalon, "\n" + test_result + "\n" + ethalon

class DotPlotModificationTest(SimpleTest):
    def testManual(self):
        self.test1()
        self.test2()
        self.test3()
        self.test4()

    def test1(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line1 = lines.addNew("ACGTAAAAGGGTACGT", "c1")
        line2 = lines.addNew("ACGTAAGGGGGTACGT", "c2")
        al = self.scorer.polyshAlignment(AlignmentPiece.Identical(line1.asSegment(), line2.asSegment()))
        dp = LineDotPlot(lines, self.aligner)
        dp.addAlignment(al)
        alignment = AlignmentPiece.Identical(Contig("AGG", "tmp").asSegment(), line2.segment(0, 3))
        line2.correctSequence([alignment])
        assert str(list(dp.alignmentsToFrom[line2.id][line1.id])) == "[(c1[0:16-0]->c2[0:16-0]:0.81)]"

    def test2(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line = lines.addNew("ACGTACGTACGT", "c")
        dp = LineDotPlot(lines, self.aligner)
        al1 = AlignmentPiece.Identical(line.segment(0, 8), line.segment(4, 12))
        al2 = AlignmentPiece.Identical(line.segment(0, 4), line.segment(8, 12))
        dp.addAlignment(al1)
        dp.addAlignment(al2)
        alignment = AlignmentPiece.Identical(Contig("AGG", "tmp").asSegment(), line.segment(4, 7))
        line.correctSequence([alignment])
        assert str(list(dp.auto_alignments[
                            "c"])) == "[(c[0:12-4]->c[4:12-0]:0.75), (c[0:4]->c[8:12-0]:1.000), (c[4:12-0]->c[0:12-4]:0.75), (c[8:12-0]->c[0:4]:1.000), (c[0:12-0]->c[0:12-0]:1.000)]"

    def test3(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line = lines.addNew("ACGTACGTACGT", "c")
        dp = LineDotPlot(lines, self.aligner)
        al1 = AlignmentPiece.Identical(line.segment(0, 8), line.segment(4, 12))
        al2 = AlignmentPiece.Identical(line.segment(0, 4), line.segment(8, 12))
        dp.addAlignment(al1)
        dp.addAlignment(al2)
        alignment = AlignmentPiece.Identical(Contig("TCC", "tmp").asSegment(), line.segment(3, 6))
        line.correctSequence([alignment])
        assert str(list(dp.auto_alignments["c"])) == "[(c[1:12-4]->c[5:12-0]:0.86), (c[0:4]->c[8:12-0]:1.000), (c[5:12-0]->c[1:12-4]:0.86), (c[8:12-0]->c[0:4]:1.000), (c[0:12-0]->c[0:12-0]:1.000)]"

    def test4(self):
        dataset = TestDataset("abcABC")
        name = dataset.addContig("abcAB")
        lines, dp, reads = dataset.genAll(self.aligner)
        line = lines[name]
        line.extendRight(dataset.alphabet["C"])
        assert str(list(dp.auto_alignments[line.id])) == "[(C0_abcAB[1650:3302-0]->C0_abcAB[0:1650]:0.995), (C0_abcAB[0:1650]->C0_abcAB[1650:3302-0]:0.995), (C0_abcAB[0:3302-0]->C0_abcAB[0:3302-0]:1.000)]"


class PotentialReadsTest(SimpleTest):
    def testManual(self):
        dataset = TestDataset("abcdefgh")
        dataset.addDisjointig("abcdefgh")
        name = dataset.addContig("abcdefgh")
        dataset.generateReads(4, 2, True)
        lines, dp, reads = dataset.genAll(self.aligner)
        line = lines[name]
        assert str(list(line.getRelevantAlignmentsFor(line.asSegment()))) == "[(R0_abcd[0:2196-4]->C0_abcdefgh[0:2195]:0.96), (R1_bcde[0:2202-0]->C0_abcdefgh[550:2750]:0.96), (R2_cdef[0:2187-0]->C0_abcdefgh[1100:3300]:0.96), (R3_defg[1:2189-0]->C0_abcdefgh[1651:3850]:0.97), (R4_efgh[0:2204-0]->C0_abcdefgh[2200:4400-0]:0.97), (R5_fgha[0:1644]->C0_abcdefgh[2750:4400-0]:0.96), (R5_fgha[1645:2200-0]->C0_abcdefgh[1:550]:0.95), (R6_ghab[1099:2198-0]->C0_abcdefgh[0:1100]:0.97), (R6_ghab[0:1099]->C0_abcdefgh[3300:4400-0]:0.97), (R7_habc[550:2192-8]->C0_abcdefgh[0:1643]:0.96), (R7_habc[0:550]->C0_abcdefgh[3850:4400-0]:0.97)]", str(list(line.getRelevantAlignmentsFor(line.asSegment())))

class UniqueRegionMarkingTest(SimpleTest):
    # def testManual(self):
        # self.test1()
        # self.test2()
        # self.test3()

    def testCase(self, instance):
        data = TokenReader(StringIO(" ".join(instance)))
        dataset = TestDataset.loadStructure(data)
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker().markAllUnique(lines, dp)
        ethalon1 = data.readToken()
        ethalon2 = data.readToken()
        line = lines[dataset.contigs[0].id]
        self.assertResult(str(line.correct_segments), ethalon1)
        self.assertResult(str(line.completely_resolved), ethalon2)

    def test1(self):
        dataset = TestDataset("abcdefgh")
        dataset.addDisjointig("abcdefgh")
        name = dataset.addContig("abcdefgh")
        dataset.generateReads(4, 15, True)
        dataset.saveStructure(TokenWriter(sys.stdout))
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker().markAllUnique(lines, dp)
        assert str(lines[name].correct_segments) == "ReadStorage+:[C0_abcdefgh[0:4400-0]]", str(lines[name].correct_segments)
        assert str(lines[name].completely_resolved) == "ReadStorage+:[C0_abcdefgh[0:4400-0]]", str(lines[name].completely_resolved)

    def test2(self):
        dataset = TestDataset("abcdefghcdeijk")
        dataset.addDisjointig("abcdefghcdeijk")
        name = dataset.addContig("abcdefgh")
        dataset.generateReads(4, 15, True)
        dataset.saveStructure(TokenWriter(sys.stdout))
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker().markAllUnique(lines, dp)
        assert str(lines[name].correct_segments) == "ReadStorage+:[C0_abcdefgh[0:4400-0]]", str(lines[name].correct_segments)
        assert str(lines[name].completely_resolved) == "ReadStorage+:[C0_abcdefgh[0:1341], C0_abcdefgh[2500:4400-0]]", str(lines[name].completely_resolved)

    def test3(self):
        dataset = TestDataset("abcdefghcdeijkeflmn")
        dataset.addDisjointig("abcdefghcdeijkeflmn")
        name = dataset.addContig("abcdefgh")
        dataset.generateReads(4, 15, True)
        dataset.saveStructure(TokenWriter(sys.stdout))
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker().markAllUnique(lines, dp)
        assert str(lines[name].correct_segments) == "ReadStorage+:[C0_abcdefgh[0:1650], C0_abcdefgh[1651:4400-0]]"
        assert str(lines[name].completely_resolved) == "ReadStorage+:[C0_abcdefgh[0:1350], C0_abcdefgh[3050:4400-0]]"


# LAUNCH TANYA!!!