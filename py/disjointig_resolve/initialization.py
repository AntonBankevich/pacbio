import itertools
import os
import subprocess
import sys

from typing import Dict, List

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import params
from common.alignment_storage import ReadCollection, AlignmentPiece
from common.dot_parser import DotParser
from common.sequences import ContigCollection, ContigStorage, Contig
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.analysis import CoverageAnalyser
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.line_extender import LineExtender
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.unique_marker import UniqueMarker


def CreateLineCollection(dir, aligner, contigs, disjointigs, reads, split):
    sys.stdout.info("Creating line collection")
    lines = NewLineStorage(disjointigs, aligner)
    if split:
        lines.splitFromContigs(contigs)
    else:
        lines.fillFromContigs(contigs)
    lines.writeToFasta(open(os.path.join(dir, "initial_lines.fasta"), "w"))
    lines.alignDisjointigs()
    # lines.fillFromDisjointigs()
    sys.stdout.info("Marking unique regions")
    UniqueMarker(aligner).markAllUnique(lines, reads)
    for line in lines.unique():
        print line, line.completely_resolved
    if not split:
        print "Splitting lines into parts"
        line_list = list(lines.unique()) # type: List[NewLine]
        while len(line_list) > 0:
            line = line_list.pop()
            print "Splitting line", line
            if len(line.completely_resolved) > 1:
                print "Splitted line because of two resolved segments:", str(line.completely_resolved)
                left = line.completely_resolved[1].left
                right = line.completely_resolved[0].right
                if left >= right:
                    right = left + 1
                line1, line2 = lines.splitLine(line.segment(left, right))
                line_list.extend([line1, line2])
            else:
                if line.initial[-1].seg_to.right + 5000 < len(line) - 5000:
                    print "Cut line of the right because too long unresolved segment:", line, str(line.completely_resolved)
                    line.cutRight(line.initial[-1].seg_to.right + 3000)
                if line.initial[0].seg_to.left > 5000:
                    print "Cut line of the left because too long unresolved segment:", line, str(line.completely_resolved)
                    line.rc.cutRight(len(line) - line.initial[0].seg_to.left + 3000)
                if len(line.completely_resolved[0]) > 40000:
                    print "Splitted line because it is too long", str(line.completely_resolved)
                    line12, line3 = lines.splitLine(line.completely_resolved[0].suffix(length=10000).prefix(length=1000))
                    line1, line2 = lines.splitLine(line12.completely_resolved[0].prefix(length=10000).suffix(length=1000))
                    line1.tie(line2, -1000, "")
                    line2.tie(line3, -1000, "")
        line_list = sorted(lines.unique(), key = lambda line: line.id)
        print "Final list of lines:"
        for line in line_list:
            print line, line.completely_resolved
    analyser = CoverageAnalyser(aligner, reads)
    sys.stdout.info("Analysing k-mer coverage by reads.")
    analyser.printAnalysis(analyser.analyseLines(lines))
    sys.stdout.info("Constructing line dot plot")
    dot_plot = LineDotPlot(lines, aligner)
    dot_plot.construct(aligner)
    dot_plot.printAll(sys.stdout)
    sys.stdout.info("Updating sequences and resolved segments.")
    knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
    extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
    extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))
    lines.writeToFasta(open(os.path.join(dir, "initial_prolonged_lines.fasta"), "w"))
    return dot_plot, extender, lines


def CreateDisjointigCollection(d_files, dir, aligner, reads):
    sys.stdout.info("Creating disjointig collection")
    disjointigs = DisjointigCollection()
    for f in d_files:
        disjointigs.loadFromFasta(open(f, "r"))
    sys.stdout.info("Extending disjointig collection")
    clen = 5000000
    bad_reads = reads.cleanCopy()
    tlen0 = sum(map(len, bad_reads))
    good_reads = set()
    for al in aligner.localAlign(reads, disjointigs):
        if not al.contradictingRTC(al.seg_to.contig.asSegment(), params.bad_end_length) and len(al.seg_from.contig) > len(al) - 2 * params.bad_end_length:
            good_reads.add(al.seg_from.contig.id)
    sys.stdout.info("Fraction of reads without full alignment to disjointigs:", 1 - float(len(good_reads)) / len(reads))
    rf = os.path.join(dir, "badreads.fasta")
    bad_reads = bad_reads.filter(lambda read: read.id not in good_reads)
    tlen = sum(map(len, bad_reads))
    bad_reads.print_fasta(open(rf, "w"))
    l = tlen * clen / tlen0
    assembly_dir = os.path.join(dir, "assembly0")
    code = subprocess.call(["./bin/flye", "--meta", "-o", assembly_dir, "-t", "8", "--" + params.technology + "-raw", rf, "--genome-size", str(l),
                           "--no-trestle", "--min-overlap", str(params.k)])
    if code == 0:
        df = os.path.join(assembly_dir, "10-consensus", "consensus.fasta")
        disjointigs.loadFromFasta(open(df, "r"))
        print "Disjointigs:"
        for dis in disjointigs:
            print dis.id, len(dis)
        disjointigs.writeToFasta(open(os.path.join(dir, "disjointigs.fasta"), "w"))
    else:
        print "Could not assemble new disjointigs"
    sys.stdout.info("Aligning reads to disjointigs")
    disjointigs.addAlignments(aligner.localAlign(reads, disjointigs))
    return disjointigs


def CreateContigCollection(graph_file, contigs_file, min_cov, aligner, polisher, reads):
    sys.stdout.info("Creating contig collection")
    if graph_file is not None:
        nonunique = [str(val[0]) for val in DotParser(open(graph_file, "r")).parse() if not (val[4].unique and val[4].cov >= min_cov)]
    else:
        nonunique = []
    contigs = ContigCollection()
    contigs.loadFromFasta(open(contigs_file, "r"), num_names=True)
    contigs = contigs.filter(lambda contig: contig.id not in nonunique and len(contig) > params.k + 20)
    sys.stdout.info("Created", len(contigs), "initial contigs")
    sys.stdout.info("Polishing contigs")
    polished_contigs = polisher.polishMany(reads, list(contigs.unique()))
    contigs = ContigCollection().addAll(polished_contigs)
    sys.stdout.info("Extending short lines")
    ExtendShortLines(contigs, reads, aligner, polisher)
    return contigs


def CreateReadCollection(reads_file, downsample):
    sys.stdout.info("Creating read collection")
    num = params.downsample
    if downsample < 1:
        sys.stdout.info("Downsampling:", downsample)
        reads = ReadCollection()
        reads.loadFromFasta(open(reads_file, "r"))
        num = int(reads.__len__() * downsample)
    reads = ReadCollection()
    reads.loadFromFasta(open(reads_file, "r"), num)
    return reads


# TODO: many short edges are not really unique. Need to address it properly
def ExtendShortLines(contigs, reads, aligner, polisher):
    # type: (ContigStorage, ReadCollection, Aligner, Polisher) -> None
    short_contigs = ContigStorage()
    als = dict() # type: Dict[str, List[AlignmentPiece]]
    for contig in contigs.unique():
        if len(contig) < params.k + 500:
            short_contigs.add(contig)
            als[contig.id] = []
            als[contig.rc.id] = []
    for al in aligner.overlapAlign(reads, short_contigs):
        if al.seg_to.left <= 20 and al.rc.seg_to.left <= 20:
            added = False
            for i, al1 in enumerate(als[al.seg_to.contig.id]):
                if al1.seg_from.contig.id == al.seg_from.contig.id:
                    added = True
                    if al.percentIdentity() > al1.percentIdentity():
                        als[al.seg_to.contig.id][i] = al
                        als[al.seg_to.contig.rc.id][i] = al.rc
                    break
            if not added:
                als[al.seg_to.contig.id].append(al)
                als[al.seg_to.contig.rc.id].append(al.rc)
    for contig in short_contigs.unique():
        if len(als[contig.id]) > 0:
            tmp_contig, new_als = polisher.polishEnd(als[contig.id], params.reliable_coverage)
            r = len(tmp_contig) - len(contig)
            tmp_contig, new_als = polisher.polishEnd([al.rc for al in new_als], params.reliable_coverage)
            l = len(tmp_contig) - len(contig) - r
        else:
            tmp_contig, new_als = contig, als[contig.id]
            l = 0
            r = 0
        if l > params.k / 2 and r > params.k / 2:
            tmp_contig.seq = tmp_contig.seq[l - params.k / 2:-r + params.k / 2]
        else:
            tmp_contig.seq = tmp_contig.seq[max(0, l - params.k):-max(1, r - params.k)]
        if len(tmp_contig) > params.k + 500:
            sys.stdout.info("Prolonged contig", contig.id, "for", l, "and", r, "nucleotides from left and right")
            contigs.add(Contig(tmp_contig.seq, contig.id))
        else:
            sys.stdout.warn("Could not prolong contig", contig.id, "enough. Removing it.")
            contigs.remove(contig)


def constructDisjointigs(reads, total_length, dir):
    # type: (ReadCollection, int, str) -> str
    reads_file = os.path.join(dir, "reads.fasta")
    disjointigs_file = os.path.join(dir, "disjointigs.fasta")
    log_file = os.path.join(dir, "log.txt")
    reads.print_fasta(open(reads_file, "w"))
    subprocess.check_call(["./bin/flye-assemble", reads_file, disjointigs_file, str(total_length), "flye/config/bin_cfg/asm_raw_reads.cfg", "-v", "1500", "-t", "16", "-u", "-l", log_file])
    return disjointigs_file
