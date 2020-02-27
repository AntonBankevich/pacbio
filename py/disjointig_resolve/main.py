import itertools
import os
import subprocess
import sys
import time
import traceback

from typing import Iterable


sys.path.append("py")
sys.path.append(".")

from common.basic import CreateLog
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.smart_storage import AlignmentStorage
from alignment.align_tools import Aligner, DirDistributor
from alignment.polishing import Polisher
from common import basic, SeqIO, params
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.tests import Tester
from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.alignment_storage import ReadCollection, AlignmentPiece
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.saves_io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.cl_params import Params, parseFlyeDir
from disjointig_resolve.initialization import CreateLineCollection, CreateDisjointigCollection, CreateContigCollection, \
    CreateReadCollection, constructDisjointigs, LoadLineCollection


def prepare_disjointigs_file(disjointigs_file, disjointigs_file_list):
    recs = []
    for fn in disjointigs_file_list:
        for rec in SeqIO.parse_fasta(open(fn, "r")):
            recs.append(rec)
    h = open(disjointigs_file, "w")
    for rec in recs:
        SeqIO.write(rec, h, "fasta")
    h.close()


def printToFile(als, dir, name):
    # type: (Iterable[AlignmentPiece], str, str) -> None
    f = open(os.path.join(dir, name + ".txt"), "w")
    for al in als:
        f.write(" ".join([str(al.seg_from.contig.id), str(al.seg_to.contig.id), str(al.seg_from.left), str(al.seg_from.right), str(al.seg_to.left), str(al.seg_to.right), str(len(list(al.splitRead())))]) + "\n")
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
    for al in aligner.overlapAlign(reads, lines):
        line = al.seg_to.contig # type: NewLine
        line.addReadAlignment(al)
        tmp.addAll(al.splitRead())
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
    start = time.time()
    cl_params = Params().parse(args)
    cl_params.check()
    CreateLog(cl_params.dir)
    print " ".join(cl_params.args)
    sys.stdout.info("Started")
    sys.stdout.info("Params:", " ".join(args))
    sys.stdout.info("Version:", subprocess.check_output(["git", "rev-parse", "HEAD"]))
    sys.stdout.info("Modifications:")
    print subprocess.check_output(["git", "diff"])
    if cl_params.test:
        aligner = Aligner(DirDistributor(cl_params.alignmentDir()))
        Tester(aligner).testAll("tests/cases.txt")
        sys.stdout.write("Finished\n")
        print (time.strftime("%d.%m.%Y  %I:%M:%S"))
        return
    if cl_params is not None:
        print "Focus:", str(cl_params.focus)
    sys.stdout.info("Preparing initial state")
    save_handler = SaveHandler(cl_params.save_dir)
    if cl_params.load_from is not None:
        # tmp = cl_params.focus
        sys.stdout.info("Loading initial state from saves")
        cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot = loadAll(TokenReader(open(cl_params.load_from, "r")))
        cl_params.parse(args)
        # cl_params.focus = tmp
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
        dot_plot.printAll(sys.stdout)
        printState(lines)
    else:
        aligner = Aligner(DirDistributor(cl_params.alignmentDir()))
        polisher = Polisher(aligner, aligner.dir_distributor)

        reads = CreateReadCollection(cl_params.reads_file, cl_params.cut_reads, cl_params.downsample)


        if cl_params.contigs_file is None:
            assembly_dir = os.path.join(cl_params.dir, "assembly_initial")
            reads_file = os.path.join(cl_params.dir, "actual_reads.fasta")
            reads.print_fasta(open(reads_file, "w"))
            subprocess.check_call(["./bin/flye", "--meta", "-o", assembly_dir, "-t", str(cl_params.threads), "--" + params.technology + "-raw", reads_file, "--genome-size", str(params.expected_size), "--min-overlap", str(params.k)])
            cl_params.set_flye_dir(assembly_dir, cl_params.mode)
        elif len(cl_params.disjointigs_file_list) == 0:
            assembly_dir = os.path.join(cl_params.dir, "assembly_initial")
            reads_file = os.path.join(cl_params.dir, "actual_reads.fasta")
            reads.print_fasta(open(reads_file, "w"))
            disjointigs_file = constructDisjointigs(reads, params.expected_size, assembly_dir)
            # graph_file, contigs_file, disjointigs_file, rep_dir, graph_file_after, contigs_file_after = parseFlyeDir(assembly_dir)
            cl_params.disjointigs_file_list.append(disjointigs_file)
            params.min_contra_for_break = 8

        contigs = CreateContigCollection(cl_params.graph_file, cl_params.contigs_file, cl_params.min_cov, aligner, polisher, reads)

        disjointigs = CreateDisjointigCollection(cl_params.disjointigs_file_list, cl_params.dir, aligner, reads)

        if cl_params.init_file is None:
            dot_plot, lines = CreateLineCollection(cl_params.dir, aligner, contigs, disjointigs, reads, cl_params.split)
        else:
            dot_plot, lines = LoadLineCollection(cl_params.dir, cl_params.init_file, aligner, contigs, disjointigs, reads)

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

    EACL(aligner, cl_params, contigs, disjointigs, dot_plot, extender, lines, reads, save_handler)

    dot_plot.printAll(sys.stdout)

    print "Final result:"
    lines.printToFasta(open(os.path.join(cl_params.dir, "lines.fasta"), "w"))
    lines.printKnottedToFasta(open(os.path.join(cl_params.dir, "lines_knotted.fasta"), "w"))
    printState(lines)

    # print "Disjointig alignments"
    # for line in lines
    sys.stdout.info("Finished")
    secs = int(time.time() - start)
    days = secs / 60 / 60 / 24
    hours = secs / 60 / 60 % 24
    mins = secs / 60 % 60
    sys.stdout.info("Finished in %d days, %d hours, %d minutes" % (days, hours, mins))


def printState(lines):
    print "Lines:"
    for line in lines:
        print line
    print "Chains:"
    for chain in lines.chains():
        if chain[-1].knot is not None:
            print "->" + "->".join([line.id for line in chain]) + "->"
        else:
            print "->".join([line.id for line in chain])


def EACL(aligner, cl_params, contigs, disjointigs, dot_plot, extender, lines, reads, save_handler):
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
                    for key in lines.items.keys():
                        if str(basic.parseNegativeNumber(basic.parseLineName(key)[-1])) == line_id and lines[key].knot is None:
                            line_id = key
                            break
                    if line_id not in lines.items.keys():
                        continue
            line = lines[line_id]
            sys.stdout.info("Investigating", line)
            extended = extender.processLine(line)
            if extended > 0:
                cnt += 1
                stop = False
            if cnt > 10:
                cnt = 0
                printState(lines)
                print "Saving current state"
                writer = save_handler.getWriter()
                print "Save details:", writer.info
                saveAll(writer, cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot)


if __name__ == "__main__":
    main(sys.argv)
# 56 73 77 167 76 12 84 1 78
