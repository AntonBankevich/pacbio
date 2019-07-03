import os
import sys

from typing import Optional, List, Iterable, Tuple

from alignment.align_tools import Aligner, DirDistributor
from common import basic, SeqIO, params
from common.line_align import Scorer
from common.seq_records import NamedSequence
from common.sequences import Consensus, Contig, ContigCollection, Segment, Contig, ContigStorage
from common.alignment_storage import AlignmentPiece, AlignedRead, ReadCollection
from flye_tools.polish import PolishException
from flye_tools.polysh_job import JobPolishing


class Polisher:
    def __init__(self, aligner, dir_distributor):
        # type: (Aligner, DirDistributor) -> Polisher
        self.aligner = aligner
        self.dir_distributor = dir_distributor

    def polish(self, reads, consensus):
        # type: (Iterable[NamedSequence], Contig) -> str
        dir, new_files, same = self.dir_distributor.fillNextDir([([consensus], "ref.fasta"), (reads, "reads.fasta")])
        consensus_file_name = new_files[0]
        reads_file_name = new_files[1]
        args = FakePolishingArgs()
        basic.ensure_dir_existance(os.path.join(dir, "work"))
        job = JobPolishing(args, os.path.join(dir, "work"), os.path.join(dir, "log.info"), [reads_file_name], consensus_file_name, "polish")
        polished_file = job.out_files["contigs"]
        if same and not params.clean and os.path.exists(polished_file):
            print "Polishing reused:", polished_file
        else:
            job.run()
        return list(SeqIO.parse_fasta(open(polished_file, "r")))[0].seq

    def polishAndAnalyse(self, reads, polishing_base, reliable_start = None):
        # type: (ReadCollection, Contig, Optional[int]) -> Consensus
        if reliable_start is None:
            reliable_start = len(polishing_base)
        seq = Contig(self.polish(reads, polishing_base), "contig")
        res = [0] * (len(seq) + 1)
        alignment = ReadCollection().extendClean(reads)
        self.aligner.alignReadCollection(alignment, [seq])
        contra = 0
        ok = 0
        late = 0
        for read in alignment:
            for al in read.alignmentsTo(seq.asSegment()):# type: AlignmentPiece
                if al.contradicting(seq.asSegment()):
                    contra += 1
                elif al.seg_to.left > reliable_start:
                    late += 1
                else:
                    res[al.seg_to.left] += 1
                    res[al.seg_to.right] -= 1
                    ok += 1
        for i in range(1, len(res)):
            res[i] += res[i - 1]
        print "Polyshed and analysed using", len(alignment), "reads. Ok:", ok, "late:", late, "contra:", contra
        if contra > 10 or contra > ok / 2:
            for read in alignment:
                print read
                for al in read.alignmentsTo(seq.asSegment()):
                    if al.contradictingRTC(seq.asSegment()):
                        print "contra_al:", al
                    elif al.seg_to.left > reliable_start:
                        print "late_al:", al
                    else:
                        print "ok_al:", al
        return Consensus(seq.seq, res)

    def polishQuiver(self, reads, base_start, pos_start, min_new_len = 3000):
        # type: (ReadCollection, str, int, int) -> Optional[Consensus]
        cc = ContigCollection([Contig(base_start, "base_start")])
        reads_to_base = ReadCollection().extendClean(reads).fillFromSam(self.aligner.align(reads, cc, "overlap"), cc)
        # for read in reads:
        #     print read
        #     if read.id in reads_to_base.reads:
        #         print reads_to_base.reads[read.id]
        #     else:
        #         print "No alignment"
        print "Polishing quiver of", len(reads_to_base), "reads."
        best = None
        for read in sorted(list(reads_to_base.__iter__()), key = lambda read: len(read))[::-1]:
            for al in read.alignments:
                if al.seg_to.contig.id == "base_start" and al.seg_to.right > len(base_start) - 50 and len(read) - al.seg_from.right > 1000:
                    base_conig = Contig(base_start[pos_start:al.seg_to.right] + read.seq[al.seg_from.right:], "base")
                    if best is not None and len(base_conig) < len(best.cut()):
                        continue
                    candidate = self.polishAndAnalyse(reads, base_conig, al.seg_to.right - pos_start)
                    if len(candidate.cut()) > len(base_start) - pos_start + min_new_len: #len(candidate.cut()) is too slow
                        print "Final polishing base alignment:", al
                        print base_conig.seq
                        print candidate.seq
                        return candidate
                    if best is None and len(candidate.cut()) > len(base_start) - pos_start or best is not None and len(candidate.cut()) > len(best.cut()):
                        print "Updated best polyshing base alignment:", al
                        print base_conig.seq
                        print candidate.seq
                        best = candidate
                    break
        if best is None:
            if len(base_start) - pos_start > 500:
                best = self.polishAndAnalyse(reads, Contig(base_start[pos_start:], "base"), len(base_start) - pos_start + 100)
            else:
                best = self.polishAndAnalyse(reads, Contig(base_start, "base"), len(base_start) + 100).suffix(pos_start)
        return best

    # Returns alignment of polished sequence to old sequence
    def polishSegment(self, seg, als):
        # type: (Segment, List[AlignmentPiece]) -> AlignmentPiece
        print "Polishing segment", seg
        w = 900
        r = 50
        first = seg.left / w
        last = min(seg.right + w - 1, len(seg.contig)) / w
        segs = []
        for i in range(first, last + 1):
            segs.append(Segment(seg.contig, max(0, i * w - r), min(len(seg.contig), (i + 1) * w + r)))
        als_by_segment = [[] for i in range(last - first + 1)]
        for al in als:
            l = al.seg_to.left / w
            r = al.seg_to.right / w + 1
            for i in range(max(0, l - first - 1), min(last - first + 1, r - first + 2)):
                if al.seg_to.inter(segs[i]):
                    als_by_segment[i].append(al)
        res_als = []
        for seg1, seg_als in zip(segs, als_by_segment):
            res_als.append(self.polishSmallSegment(seg1, seg_als))
        # print "Res als", res_als
        res = AlignmentPiece.GlueOverlappingAlignments(res_als)
        return res.reduce(target=seg)

    # Returns alignment of polished sequence to old sequence
    def polishSmallSegment(self, seg, als):
        # type: (Segment, List[AlignmentPiece]) -> AlignmentPiece
        ok = False
        for al in als:
            if al.seg_to.contains(seg):
                ok = True
        if not ok:
            sys.stdout.log(params.LogPriority.warning, "Warning", seg, "has no covering reads")
            return AlignmentPiece.Identical(seg.asContig().asSegment(), seg)
        reads = []
        start = basic.randomSequence(200)
        end = basic.randomSequence(200)
        for al in als:
            new_seq = ""
            al = al.reduce(target=seg)
            if al.seg_to.left < seg.left + 3:
                new_seq += start
            new_seq += al.seg_from.Seq()
            if al.seg_to.right > seg.right - 3:
                new_seq += end
            reads.append(NamedSequence(new_seq, al.seg_from.contig.id))
        base = Contig(start + seg.Seq() + end, "base")
        polished = None
        try:
            polished = Contig(self.polish(reads, base), "polished")
        except PolishException:
            sys.stdout.log(params.LogPriority.warning, "Warning", seg, "has a sequence very different from reads. Using reads to correct.")
            for al, read in zip(als, reads):
                if al.seg_to.contains(seg):
                    try:
                        polished = Contig(self.polish(reads, Contig(read.seq, read.id)), "polished")
                        break
                    except PolishException:
                        pass
        if polished is None:
            sys.stdout.log(params.LogPriority.warning, "Warning", seg, "could not be corrected even though some reads cover it.")
            polished = seg.asContig()
        al = self.aligner.overlapAlign([polished], ContigStorage([base])).next()
        mapping = AlignmentPiece.Identical(base.segment(len(start), len(base) - len(end)), seg)
        return al.compose(mapping)

    def polishEnd(self, als, min_cov = 4):
        # type: (List[AlignmentPiece], int) -> Tuple[Contig, List[AlignmentPiece]]
        scorer = Scorer()
        contig = als[0].seg_to.contig
        new_contig = contig.asSegment().asContig()
        relevant_als = [al.changeTargetContig(new_contig) for al in als if al.rc.seg_to.left < 100]
        finished_als = []
        while True:
            tmp = []
            for al in relevant_als:
                if al.seg_to.inter(new_contig.asSegment().suffix(length=20)) and al.rc.seg_from.left > 100:
                    tmp.append(al)
                else:
                    finished_als.append(al)
            relevant_als = tmp
            # TODO replace with position search in cigar
            if len(relevant_als) < min_cov:
                break
            start = basic.randomSequence(params.flanking_size) + new_contig.asSegment().suffix(length=params.flanking_size).Seq()
            reduced_read_list = [
                AlignedRead.new(start + al.seg_from.contig.asSegment().suffix(pos=al.seg_from.right).Seq(), str(i) + "_" + al.seg_from.contig.id)
                for i, al in enumerate(relevant_als)]
            # print "RRL:", reduced_read_list
            reduced_reads = ReadCollection(reduced_read_list)
            found = False
            for base_al in relevant_als:
                if base_al.rc.seg_from.left < params.flanking_size:
                    continue
                # Base consists of copy of the previous 200 nucleotides and a segment of read of length at most 500
                base_segment = base_al.seg_from.contig.segment(base_al.seg_from.right,
                                                     min(len(base_al.seg_from.contig), base_al.seg_from.right + params.window_size))
                base = Contig(start + base_segment.Seq(), "base")
                for read in reduced_read_list:
                    read.clean()
                polished_base = Contig(self.polish(reduced_reads, base), "polished_base")
                for al in self.aligner.localAlign(reduced_reads, ContigStorage().addAll([polished_base])):
                    reduced_reads.reads[al.seg_from.contig.id].addAlignment(al)
                # self.aligner.alignReadCollection(reduced_reads, [polished_base])
                candidate_alignments = []
                # print "RRL1:", reduced_read_list
                for read in reduced_read_list:
                    candidate_alignments.append(None)
                    for al in read.alignmentsTo(polished_base.asSegment()):
                        if al.seg_to.left == 0 and ((candidate_alignments[-1] is None or candidate_alignments[-1].seg_to.right < al.seg_to.right)):
                            candidate_alignments[-1] = al
                positions = []
                # print "CA:", candidate_alignments
                for i, al in enumerate(candidate_alignments):
                    if al is None:
                        print al, reduced_read_list[i]
                        print reduced_read_list[i].alignments
                        print reduced_read_list[i].seq
                        print polished_base.seq
                    assert al is not None, reduced_read_list[i]
                    positions.append(al.seg_to.right)
                positions = sorted(positions)[::-1]
                num = max(min_cov, len(relevant_als)  * 8 / 10)
                if num >= len(positions):
                    continue
                print "Chose to use read", base_al, "Alignments:"
                print map(str, reduced_read_list)
                cutoff_pos = max(positions[num - 1], len(start))
                if cutoff_pos > len(start) + 100:
                    found = True
                    new_contig_candidate = Contig(new_contig.seq + polished_base[len(start):cutoff_pos], "candidate")
                    embedding = AlignmentPiece.Identical(polished_base.segment(len(start), cutoff_pos), new_contig_candidate.asSegment().suffix(pos=len(new_contig)))
                    read_mappings = []
                    for al1, al2 in zip(candidate_alignments, relevant_als):
                        seg_from = al2.seg_from.contig.asSegment().suffix(length = len(al1.seg_from.contig) - len(start))
                        seg_to = al1.seg_from.contig.asSegment().suffix(length = len(al1.seg_from.contig) - len(start))
                        read_mappings.append(AlignmentPiece.Identical(seg_from, seg_to))
                    # print [(al2, al1, embedding) for al1, al2 in zip(candidate_alignments, read_mappings)]
                    embedded_alignments = []
                    for al1, al2 in zip(candidate_alignments, read_mappings):
                        # print al2, al1, embedding
                        # print al2.compose(al1)
                        # print al2.compose(al1).compose(embedding)
                        if al1.seg_from.right <= len(start) + 5:
                            embedded_alignments.append(None)
                        else:
                            embedded_alignments.append(al2.compose(al1).compose(embedding))
                    # candidate_alignments = [al2.compose(al1).compose(embedding) for al1, al2 in zip(candidate_alignments, read_mappings)]
                    corrected_relevant_alignments = [al.targetAsSegment(new_contig_candidate.asSegment().prefix(len(new_contig))) for al in relevant_als]
                    relevant_als = []
                    for al1, al2 in zip(corrected_relevant_alignments, embedded_alignments):
                        if al2 is None:
                            al = al1
                        else:
                            al = al1.mergeDistant(al2)
                            if al is None:
                                al = al1
                            elif al1.seg_from.dist(al2.seg_from) >= 10 or al1.seg_to.dist(al2.seg_to) >= 10:
                                al = scorer.polyshAlignment(al)
                        relevant_als.append(al)
                    finished_als = [al.targetAsSegment(new_contig_candidate.asSegment().prefix(len(new_contig))) for al in finished_als]
                    new_contig = new_contig_candidate
                    break
            if not found:
                break
        return new_contig, relevant_als + finished_als


class FakePolishingArgs:
    def __init__(self):
        self.num_iters = params.num_iters
        self.platform = "pacbio"
        self.threads = 8
