import os
import sys

sys.path.append("py")

from common import basic, SeqIO
from common.seq_records import NamedSequence
from disjointig_resolve import dataset_simulation


def simulate(dir, mutation, genome):
    print "Simulating", genome
    ds = dataset_simulation.TestDataset(genome, 5000, mutation_rate=mutation)
    genome_seq = ds.mutate(ds.genome, mutation / 2)[0]
    f = open(os.path.join(dir, genome + str(mutation * 100) + ".fasta"), "w")
    SeqIO.write(NamedSequence(genome_seq, genome), f, "fasta")
    f.close()


if __name__ == "__main__":
    dir = sys.argv[1]
    mutation = float(sys.argv[2])
    basic.ensure_dir_existance(dir)
    genomes = sys.argv[3:]
    for genome in genomes:
        simulate(dir, mutation, genome)


#AbcdefghIbcdefghJ AbcdefIcdefJdefKefHfL LabcdefHabmdefIamcdefJabcdmfK LabcdefHabdefIacdefJabcdfK LabcdefHabgdefIagcdefJabcdgfK
