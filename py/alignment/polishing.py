import os

from typing import Optional

from alignment.align_tools import Aligner, DirDistributor
from common import basic, SeqIO
from dag_resolve import params
from dag_resolve.sequences import Consensus, ReadCollection, Contig, ContigCollection
from flye.polysh_job import JobPolishing


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

    def polishAndAnalyse(self, reads, polishing_base):
        # type: (ReadCollection, Contig) -> Consensus
        seq = Contig(self.polish(reads, polishing_base), 0)
        res = [0] * (len(seq) + 1)
        for rec in self.aligner.align(reads, ContigCollection([seq])):
            if rec.is_unmapped:
                continue
            res[rec.pos - 1] += 1
            res[rec.pos + rec.alen - 1] -= 1
        for i in range(1, len(res)):
            res[i] += res[i - 1]
        return Consensus(seq.seq, res)

    def polishQuiver(self, reads, base_start, pos_start, min_new_len = 1000):
        # type: (ReadCollection, str, int) -> Optional[Consensus]
        cc = ContigCollection([Contig(base_start, 0)])
        reads_to_base = ReadCollection(cc).loadFromSam(self.aligner.align(reads, cc))
        print "Polishing quiver of", len(reads_to_base), "reads."
        best = None
        for read in sorted(list(reads_to_base.__iter__()), key = lambda read: len(read))[::-1]:
            for al in read.alignments:
                if al.rc:
                    continue
                if al.seg_to.right > len(base_start) - 50 and len(read) - al.seg_from.right > 1000:
                    tmp = self.polishAndAnalyse(reads, Contig(base_start[pos_start:al.seg_to.right] + read.seq[al.seg_from.right:], 0))
                    if len(tmp.cut()) > len(base_start) - pos_start + min_new_len:
                        # print read.__str__(), al.seg_to.right, len(base_start), pos_start, al.seg_from.right, len(tmp)
                        return tmp
                    if best is None and len(tmp.cut()) > len(base_start) - pos_start or best is not None and len(tmp.cut()) > len(best.cut()):
                        # print read.__str__(), al.seg_to.right, len(base_start), pos_start, al.seg_from.right, len(tmp)
                        best = tmp
                    break
        if best is None:
            best = self.polishAndAnalyse(reads, Contig(base_start[pos_start:], 0))
        return best


class FakePolishingArgs:
    def __init__(self):
        self.num_iters = params.num_iters
        self.platform = "pacbio"
        self.threads = 8