import os
import random
import sys

sys.path.append("py")

from common import basic, SeqIO
from common.seq_records import NamedSequence
from disjointig_resolve import dataset_simulation


def simulate1(dir, mutation, genome):
    print "Simulating", genome
    ds = dataset_simulation.TestDataset(genome, 5000, mutation_rate=mutation)
    genome_seq = ds.mutate(ds.genome, mutation / 2)[0]
    f = open(os.path.join(dir, genome + str(mutation * 100) + ".fasta"), "w")
    SeqIO.write(NamedSequence(genome_seq, genome), f, "fasta")
    f.close()

def simulate2(dir, mutation, error_rate, genome):
    print "Simulating", genome
    ds = dataset_simulation.TestDataset(genome, 4000, mutation_rate=mutation)
    genome_seq = ds.mutate(ds.genome, mutation / 2)[0]
    total = 0
    f = open(os.path.join(dir, genome + ".fasta"), "w")
    SeqIO.write(NamedSequence(genome_seq, genome), f, "fasta")
    f.close()
    f = open(os.path.join(dir, "reads.fasta"), "w")
    cnt = 0
    while total < len(genome_seq) * 30:
        l = random.randint(3000, 3500)
        pos = random.randint(0, len(genome_seq) - l)
        seq = ds.mutate(genome_seq[pos:pos + l], error_rate)[0]
        SeqIO.write(NamedSequence(seq, str(cnt)), f, "fasta")
        cnt += 1
        total += len(seq)
    f.close()


if __name__ == "__main__":
    dir = sys.argv[1]
    mutation = float(sys.argv[2])
    error_rate = float(sys.argv[3])
    basic.ensure_dir_existance(dir)
    genomes = sys.argv[4:]
    for genome in genomes:
        simulate2(dir, mutation, error_rate, genome)


#AbcdefghIbcdefghJ AbcdefIcdefJdefKefHfL LabcdefHabmdefIamcdefJabcdmfK LabcdefHabdefIacdefJabcdfK LabcdefHabgdefIagcdefJabcdgfK
