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


def CreateLineCollection(aligner, contigs, disjointigs, reads, split):
    sys.stdout.info("Creating line collection")
    lines = NewLineStorage(disjointigs, aligner)
    if split:
        lines.splitFromContigs(contigs)
    else:
        lines.fillFromContigs(contigs)
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
    if not split:
        line_list = list(lines.unique()) # type: List[NewLine]
        while len(line_list) > 0:
            line = line_list.pop()
            if len(line.completely_resolved) > 1:
                left = line.completely_resolved[1].left
                right = line.completely_resolved[0].right
                if left >= right:
                    right = left + 1
                line1, line2 = lines.splitLine(line.segment(left, right))
                line_list.extend([line1, line2])
            elif len(line.completely_resolved[0]) > 40000:
                line1, line2 = lines.splitLine(line.completely_resolved[0].prefix(length=1000))
                line2, line3 = lines.splitLine(line2.completely_resolved[0].suffix(length=1000))
                line1.tie(line2, -1000, "")
                line2.tie(line3, -1000, "")
        line_list = sorted(lines.unique(), key = lambda line: line.id)
        print "Final list of lines:"
        for line in line_list:
            print line, line.completely_resolved
    analyser = CoverageAnalyser(aligner, reads)
    sys.stdout.info("Analysing k-mer coverage by reads.")
    analyser.printAnalysis(analyser.analyseLines(lines))
    sys.stdout.info("Updating sequences and resolved segments.")
    knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
    extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
    extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))
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
        if not al.contradictingRTC(al.seg_to.contig.asSegment(), 500):
            good_reads.add(al.seg_from.contig.id)
    sys.stdout.info("Fraction of reads without full alignment to disjointigs:", 1 - float(len(good_reads)) / len(reads))
    rf = os.path.join(dir, "badreads.fasta")
    bad_reads = bad_reads.filter(lambda read: read.id not in good_reads)
    tlen = sum(map(len, bad_reads))
    bad_reads.print_fasta(open(rf, "w"))
    l = tlen * clen / tlen0
    assembly_dir = os.path.join(dir, "assembly0")
    subprocess.check_call(["./bin/flye", "-o", assembly_dir, "-t", "8", "--pacbio-raw", rf, "--genome-size", str(l),
                           "--no-trestle"])
    df = os.path.join(assembly_dir, "10-consensus", "consensus.fasta")
    disjointigs.loadFromFasta(open(df, "r"))
    print "Disjointigs:"
    for dis in disjointigs:
        print dis.id, len(dis)
    sys.stdout.info("Aligning reads to disjointigs")
    disjointigs.addAlignments(aligner.localAlign(reads, disjointigs))
    return disjointigs


def CreateContigCollection(graph_file, contigs_file, min_cov, aligner, polisher, reads):
    sys.stdout.info("Creating contig collection")
    unique = [str(val[0]) for val in DotParser(open(graph_file, "r")).parse() if val[4].unique and val[4].cov >= min_cov]
    contigs = ContigCollection()
    contigs.loadFromFasta(open(contigs_file, "r"), num_names=True)
    contigs = contigs.filter(lambda contig: contig.id in unique)
    sys.stdout.info("Created", len(contigs), "initial contigs")
    sys.stdout.info("Polishing contigs")
    polished_contigs = polisher.polishMany(reads, list(contigs.unique()))
    contigs = ContigCollection().addAll(polished_contigs)
    sys.stdout.info("Extending short lines")
    ExtendShortLines(contigs, reads, aligner, polisher)
    return contigs


def CreateReadColection(reads_file):
    sys.stdout.info("Creating read collection")
    reads = ReadCollection()
    reads.loadFromFasta(open(reads_file, "r"), downsample=params.downsample)
    return reads


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
        als[al.seg_to.contig.id].append(al)
        als[al.seg_to.contig.rc.id].append(al.rc)
    for contig in short_contigs.unique():
        tmp_contig, new_als = polisher.polishEnd(als[contig.id], params.reliable_coverage)
        r = len(tmp_contig) - len(contig)
        tmp_contig, new_als = polisher.polishEnd([al.rc for al in new_als], params.reliable_coverage)
        l = len(tmp_contig) - len(contig) - r
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