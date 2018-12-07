import os

from typing import Optional

from alignment.align_tools import Aligner, DirDistributor
from common import basic, SeqIO
from dag_resolve import params
from common.sequences import Consensus, ReadCollection, Contig, ContigCollection, AlignmentPiece
from flye_tools.polysh_job import JobPolishing


class Polisher:
    def __init__(self, aligner, dir_distributor):
        # type: (Aligner, DirDistributor) -> Polisher
        self.aligner = aligner
        self.dir_distributor = dir_distributor

    def polish(self, reads, consensus):
        # type: (ReadCollection, Contig) -> str
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
        alignment = ReadCollection(ContigCollection([seq])).extendClean(reads)
        self.aligner.alignReadCollection(alignment)
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
        reads_to_base = ReadCollection(cc).extendClean(reads).fillFromSam(self.aligner.align(reads, cc))
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


class FakePolishingArgs:
    def __init__(self):
        self.num_iters = params.num_iters
        self.platform = "pacbio"
        self.threads = 8