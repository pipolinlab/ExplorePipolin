import os
import sys
from collections import namedtuple
import subprocess

os.getcwd()
os.chdir('./tfm/NCBI_strains_analysis')

with open('pipolins.csv') as inf:
    fields = inf.readline().strip().split(',')
    PolB_Genomes = namedtuple('PolB_Genomes', fields, defaults=(None,) * len(fields))
    polB_genomes = []
    for line in inf:
        genome = line.strip().split(',')
        polB_genomes.append(PolB_Genomes(id=genome[0], polB_s=genome[1], polB_e=genome[2],
                                         attL_s=genome[3], attL_e=genome[4],
                                         attR_s=genome[5], attR_e=genome[6],
                                         attM_s=genome[7], attM_e=genome[8]))

# abricate resfinder
os.mkdir('./abricate_resfinder')
print(sys.path)

for i_g, genome in enumerate(polB_genomes):
    print(i_g)
    subprocess.run(['/home/liubov/repos/abricate/bin/abricate', '--db', 'resfinder',
                    f'pipolin_regions/{genome.id}_selected.fa',
                    '>', f'./abricate_resfinder/{genome.id}.tab'])
# ERROR: Could not find 'any2fasta'. Please install it and ensure it is in the PATH.
