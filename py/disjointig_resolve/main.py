import itertools
import os
import shutil
import subprocess
import sys
import time
import traceback

from typing import Iterable

sys.path.append("py")


from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.smart_storage import AlignmentStorage
from common.dot_parser import DotParser
from disjointig_resolve.unique_marker import UniqueMarker
from alignment.align_tools import Aligner, DirDistributor
from alignment.polishing import Polisher
from common import basic, SeqIO, sam_parser, params
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.tests import Tester
from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.sequences import ContigCollection, Contig, ContigStorage
from common.alignment_storage import ReadCollection, AlignmentPiece, AlignedRead
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.saves_io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.cl_params import Params


def CreateLog(dir):
    old_logs_dir = os.path.join(dir, "old")
    basic.ensure_dir_existance(old_logs_dir)
    log_file = os.path.join(dir, "log.info")
    if os.path.isfile(log_file):
        num = len(os.listdir(old_logs_dir))
        shutil.copy(log_file, os.path.join(old_logs_dir, str(num) + ".log"))
    log = open(log_file, "w")
    sys.stdout = basic.OStreamWrapper(sys.stdout, log)
    sys.stdout.prefix = lambda s: time.strftime("%I:%M:%S") + "  "
    sys.stderr = sys.stdout


def prepare_disjointigs_file(disjointigs_file, disjointigs_file_list):
    recs = []
    for fn in disjointigs_file_list:
        for rec in SeqIO.parse_fasta(open(fn, "r")):
            rec.append(rec)
    h = open(disjointigs_file, "w")
    for rec in recs:
        SeqIO.write(rec, h, "fasta")
    h.close()


def printToFile(als, dir, name):
    # type: (Iterable[AlignmentPiece], str, str) -> None
    f = open(os.path.join(dir, name + ".txt"), "w")
    for al in als:
        f.write(" ".join([str(al.seg_from.contig.id), str(al.seg_to.contig.id), str(al.seg_from.left), str(al.seg_from.right), str(al.seg_to.left), str(al.seg_to.right), str(len(list(al.split())))]) + "\n")
    f.close()

def printVals(vals, dir, name):
    # type: (Iterable[int], str, str) -> None
    f = open(os.path.join(dir, name + ".txt"), "w")
    for val in vals:
        f.write(str(val))
    f.close()


def countStats(reads, lines, disjointigs, aligner, dir):
    # type: (ReadCollection, NewLineStorage, DisjointigCollection, Aligner, str) -> None
    tmp = AlignmentStorage()
    cnt = 0
    for al in aligner.alignClean(reads, lines):
        line = al.seg_to.contig # type: NewLine
        line.addReadAlignment(al)
        tmp.addAll(al.split())
        cnt += 1
        if cnt % 10000 == 0:
            print cnt
    line = lines["1"]
    als = line.read_alignments
    print list(als.calculateWindowedCoverage(500))
    als = line.read_alignments.filter(lambda al: al.contradictingRTC(line.asSegment(), 500))
    print list(als.calculateWindowedCoverage(500))
    print list(tmp.calculateWindowedCoverage(500))
    als = tmp.filter(lambda al: al.contradictingRTC(line.asSegment(), 500))
    print list(als.calculateWindowedCoverage(500))








def main(args):
    cl_params = Params().parse(args)
    cl_params.check()
    CreateLog(cl_params.dir)
    print " ".join(cl_params.args)
    sys.stdout.info("Started")
    sys.stdout.info("Params:", " ".join(args))
    if cl_params.test:
        aligner = Aligner(DirDistributor(cl_params.alignmentDir()))
        Tester(aligner).testAll("tests/cases.txt")
        sys.stdout.write("Finished\n")
        print (time.strftime("%d.%m.%Y  %I:%M:%S"))
        return
    sys.stdout.info("Preparing initial state")
    save_handler = SaveHandler(cl_params.save_dir)
    if cl_params.load_from is not None:
        # if cl_params.new_disjointigs:
        #     new_disjointigs = cl_params.disjointigs_file_list
        # else:
        #     new_disjointigs = []
        # print new_disjointigs
        sys.stdout.info("Loading initial state from saves")
        cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot = loadAll(TokenReader(open(cl_params.load_from, "r")))
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
        # for line in lines:
        #     if line.knot is None and line.id.endswith("l") and not line.id.startswith("-"):
        #         other = lines[line.id[:-1] + "r"]
        #         line.tie(other, 0, "")
        # writer = save_handler.getWriter()
        # print "Save details:", writer.info
        # saveAll(writer, cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot)
        # dot_plot = LineDotPlot(lines, aligner)
        # dot_plot.construct(aligner)
        dot_plot.printAll(sys.stdout)
        # dp1 = LineDotPlot(lines, aligner)
        # dp1.construct(aligner)
        # dp1.printAll(sys.stdout)
        # dot_plot = dp1
        # if len(new_disjointigs) > 0:
        #     sys.stdout.info("Creating disjointig collection")
        #     for line in lines:
        #         line.disjointig_alignments.clean()
        #     for f in new_disjointigs:
        #         disjointigs.loadFromFasta(open(f, "r"))
        #     print "Disjointigs:"
        #     print list(disjointigs)
        #     lines.alignDisjointigs()
        #
        #     sys.stdout.info("Aligning reads to disjointigs")
        #     disjointigs.addAlignments(aligner.alignClean(reads, disjointigs))
        #     good_reads = set()
        #     for dis in disjointigs:
        #         for al in dis.read_alignments:
        #             good_reads.add(al.seg_from.contig.id)
        #             good_reads.add(al.seg_from.contig.rc.id)
        #     bad_reads = []
        #     brf = open(os.path.join(cl_params.dir, "br.fasta"), "w")
        #     bad = 0
        #     for read in reads:
        #         if read.id not in good_reads:
        #             bad += 1
        #             bad_reads.append(read)
        #             SeqIO.write(read, brf, "fasta")
        #     brf.close()
        #     sys.stdout.info("Fraction of reads without full alignment to disjointigs:", float(bad) / len(reads))


    else:
        aligner = Aligner(DirDistributor(cl_params.alignmentDir()))
        # polisher = Polisher(aligner, aligner.dir_distributor)
        # consensus = SeqIO.parse_fasta(open("results/NCTC9002/alignment/148/ref.fasta", "r")).next()
        # consensus = Contig(consensus.seq, consensus.id)
        # for al in sam_parser.Samfile(open("results/NCTC9002/alignment/148/work/polish/minimap_1.sam", "r")):
        #     al = AlignmentPiece.FromSamRecord(Contig(al.seq, "1"), consensus, al)
        #     print al
        #     print "\n".join(al.asMatchingStrings())
        # polisher.polish(ReadCollection().loadFromFasta(open("results/NCTC9002/alignment/148/reads.fasta", "r")), consensus)


        sys.stdout.info("Creating read collection")
        reads = ReadCollection()
        reads.loadFromFasta(open(cl_params.reads_file, "r"), downsample=params.downsample)
        # if cl_params.stats:
        #     countStats(reads, lines, disjointigs, aligner, cl_params.dir)
        #     return

        sys.stdout.info("Creating disjointig collection")
        # prepare_disjointigs_file(cl_params.disjointigs_file, cl_params.disjointigs_file_list)
        disjointigs = DisjointigCollection()
        for f in cl_params.disjointigs_file_list:
            disjointigs.loadFromFasta(open(f, "r"))

        sys.stdout.info("Extending disjointig collection")
        clen = 5000000
        bad_reads = reads.cleanCopy()
        tlen0 = sum(map(len, bad_reads))
        good_reads = set()
        for al in aligner.alignAndSplit(reads, disjointigs):
            if not al.contradictingRTC(al.seg_to.contig.asSegment(), 500):
                good_reads.add(al.seg_from.contig.id)
        sys.stdout.info("Fraction of reads without full alignment to disjointigs:", 1 - float(len(good_reads)) / len(reads))
        rf = os.path.join(cl_params.dir, "badreads.fasta")
        bad_reads = bad_reads.filter(lambda read: read.id not in good_reads)
        tlen = sum(map(len, bad_reads))
        bad_reads.print_fasta(open(rf, "w"))
        l = tlen * clen / tlen0
        assembly_dir = os.path.join(cl_params.dir, "assembly0")
        subprocess.check_call(["./bin/flye", "-o", assembly_dir, "-t", "8", "--pacbio-raw", rf, "--genome-size", str(l),
             "--no-trestle"])
        df = os.path.join(assembly_dir, "10-consensus", "consensus.fasta")
        disjointigs.loadFromFasta(open(df, "r"))

        sys.stdout.info("Aligning reads to disjointigs")
        disjointigs.addAlignments(aligner.alignClean(reads, disjointigs))

        sys.stdout.info("Creating contig collection")
        unique = [str(val[0]) for val in DotParser(open(cl_params.graph_file, "r")).parse() if val[4].unique]
        contigs = ContigCollection()
        contigs.loadFromFasta(open(cl_params.contigs_file, "r"), num_names=True)
        contigs = contigs.filter(lambda contig: contig.id in unique)
        sys.stdout.info("Created", len(contigs), "initial contigs")

        sys.stdout.info("Creating line collection")
        lines = NewLineStorage(disjointigs, aligner)
        # lines.fillFromContigs(contigs)
        lines.splitFromContigs(contigs)
        lines.alignDisjointigs()
        # lines.fillFromDisjointigs()

        sys.stdout.info("Constructing line dot plot")
        dot_plot = LineDotPlot(lines, aligner)
        dot_plot.construct(aligner)
        dot_plot.printAll(sys.stdout)

        sys.stdout.info("Marking unique regions")
        UniqueMarker(aligner).markAllUnique(lines, dot_plot, reads)
        for line in lines.unique():
            print line, line.completely_resolved

        sys.stdout.info("Updating sequences and resolved segments.")
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
        extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))

        print "Saving initial state"
        try:
            writer = save_handler.getWriter()
            print "Save details:", writer.info
            saveAll(writer, cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot)
        except Exception as e:
            _, _, tb = sys.exc_info()
            sys.stdout.warn("Could not write save")
            traceback.print_tb(tb)
            print "Message:", e.message

    print "Disjointig alignments"
    for line in lines:
        print line.disjointig_alignments
    sys.stdout.info("Resolving")
    # for line in lines:
    #     cov = 0
    #     for al in line.read_alignments:
    #         cov += len(al.seg_to)
    #     print line.id, float(cov) / len(line)
    # while "52" in lines.items and "-11l" in lines.items:
    #     extender.tryExtend(lines["52"])
    #     if "-11l" in lines.items:
    #         extender.tryExtend(lines["-11l"])
    # return
    # extender.tryExtend(lines["-24"])
    # print "initial", lines["-16l"].initial
    # print "resolved", lines["-16l"].completely_resolved
    # print "initial", lines["-15l"].initial
    # print "resolved", lines["-15l"].completely_resolved
    # print "dotplot:", list(dot_plot.alignmentsToFrom["-16l"]["-15l"])
    # print "correct", lines["-16l"].correct_segments
    # print "disjointig alignments", lines["-16l"].disjointig_alignments
    # # print "Reads to disjointigs", list(disjointigs["-D1"].read_alignments.allInter(disjointigs["-D1"].segment(4953235-22220,4953235-17220)))
    # print "Selected read alignments", list(lines["-16l"].read_alignments.allInter(lines["-16l"].asSegment().suffix(length=4000)))
    # print "Detected read alignments", list(lines["-16l"].getRelevantAlignmentsFor(lines["-16l"].completely_resolved[-1].suffix(length=1000)))
    # print "Bruteforce:", list(AlignmentStorage().addAll(aligner.alignClean(reads, ContigStorage().addAll([lines["-16l"]]))).allInter(lines["-16l"].completely_resolved[-1].suffix(length=1000)))
    # print "Bruteforce:", list(AlignmentStorage().addAll(aligner.alignAndSplit(reads, ContigStorage().addAll([lines["-16l"]]))).allInter(lines["-16l"].completely_resolved[-1].suffix(length=1000)))
    # extender.updateAllStructures(lines["-16l"].completely_resolved)
    cnt = 0
    stop = False
    while not stop:
        stop = True
        if cl_params.focus is None:
            keys = list(lines.items.keys())
        else:
            keys = cl_params.focus
        for line_id in keys:
            if line_id not in lines.items:
                if cl_params.focus is None:
                    continue
                else:
                    for key in lines.items.keys:
                        if basic.parseLineName(key)[-1].startswith(line_id):
                            line_id = key
                            break
            line = lines[line_id]
            sys.stdout.info("Investigating", line)
            extended = extender.tryExtend(line)
            if extended:
                cnt += 1
                stop = False
            if cnt > 10:
                cnt = 0
                print "Saving current state"
                writer = save_handler.getWriter()
                print "Save details:", writer.info
                saveAll(writer, cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot)
    lines.printToFasta(open(os.path.join(cl_params.dir, "lines.fasta"), "w"))
    # print "Disjointig alignments"
    # for line in lines
    sys.stdout.info("Finished")


if __name__ == "__main__":
    main(sys.argv)