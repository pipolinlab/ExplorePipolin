import os
import shutil
from collections import namedtuple
import subprocess

# print(os.getcwd())
# # go to the working directory
# os.chdir('tfm/NCBI_strains_analysis')
#
#
# def parse_blast_hittable(ht_file, id=True):
#     # :return: list of tuples with subject id, start and end positions
#     with open(ht_file) as inf:
#         entries = []
#         for line in inf:
#             if line[0] != '#':
#                 entry = line.strip().split(sep='\t')
#                 if id:
#                     entries.append((entry[1], entry[8], entry[9]))
#                 else:
#                     entries.append((entry[8], entry[9]))
#     return entries


# def get_seq_from_fasta(f_file):
#     with open(f_file) as inf:
#         seq = ''
#         for line in inf:
#             if line[0] != '>':
#                 seq += line.strip()
#     return seq


# # create a namedtuple object to store info
# fields = ('id', 'polB_s', 'polB_e', 'attL_s', 'attL_e', 'attR_s', 'attR_e', 'attM_s', 'attM_e')
# PolB_Genomes = namedtuple('PolB_Genomes', fields, defaults=(None, ) * len(fields))

# # Parse blast output to get genome ids and start/end positions of piPolBs
# polB_genomes = []
# for hit in parse_blast_hittable('./blastN_hittable_oct16.txt'):
#     polB_genomes.append(PolB_Genomes(id=hit[0], polB_s=hit[1], polB_e=hit[2]))
# # print(len(polB_genomes))   # 72
# # print(polB_genomes[0].id)

# # Download all the required genomes
# os.mkdir('genomes')
# os.chdir('./genomes')
# for genome in polB_genomes:
#     subprocess.run(['/home/liubov/.local/bin/ncbi-acc-download', '-F', 'fasta', genome.id], )
# os.chdir('..')

# # Reference pipolin
# ref_pipolin = PolB_Genomes(id='NZ_JNMI01000006.1', polB_s=80192, polB_e=82792)
# # from gtfm/from_Modesto/annotation_pipolins.xlsx
# ref_pipolin = ref_pipolin._replace(attL_s=64242)
# ref_pipolin = ref_pipolin._replace(attL_e=64374)
# ref_pipolin = ref_pipolin._replace(attR_s=89933)
# ref_pipolin = ref_pipolin._replace(attR_e=90063)

# # Since attL and attR are very similar, I will use only one of them to blast for atts
# ref_seq = get_seq_from_fasta('./NZ_JNMI01000006.1.fa')
# # print(len(ref_seq))   # 165336
# # extract "reference" att sequence
# with open('ref_att.fa', 'w') as ouf:
#     print('>att-site', file=ouf)
#     print(ref_seq[ref_pipolin.attL_s - 1:ref_pipolin.attL_e - 1], file=ouf)

# # BLAST each genome against att
# os.mkdir('att_blast')
# for genome in polB_genomes:
#     with open(f'./att_blast/{genome.id}_hits.txt', 'w') as ouf:
#         subprocess.run(['blastn', '-query', f'{ref_pipolin.id}_att.fa',
#                         '-subject', f'./genomes/{genome.id}.fa', '-outfmt', '6'],
#                        stdout=ouf)
# # Parse blast output
# for i_g, genome in enumerate(polB_genomes):
#     att_regions = parse_blast_hittable(f'./att_blast/{genome.id}_hits.txt', id=False)
#     att_regions.sort()   # we can since the att regions are not overlapping in our case
#     # check whether polB within the att regions:
#     if polB_genomes[i_g].polB_s >= sorted(att_regions[0])[0]:
#         if polB_genomes[i_g].polB_e <= sorted(att_regions[-1])[1]:
#             polB_genomes[i_g] = polB_genomes[i_g]._replace(attL_s=att_regions[0][0],
#                                                            attL_e=att_regions[0][1],
#                                                            attR_s=att_regions[-1][0],
#                                                            attR_e=att_regions[-1][1])
#             if len(att_regions) == 3:
#                 polB_genomes[i_g] = polB_genomes[i_g]._replace(attM_s=att_regions[1][0],
#                                                                attM_e=att_regions[1][1])
#         else:
#             continue
#     else:
#         continue


# copy the contig with reference pipolin as a genome
# shutil.copyfile('./NZ_JNMI01000006.1.fa', './genomes/NZ_JNMI01000006.1.fa')
# # add reference pipolin as well
# polB_genomes.append(ref_pipolin)
# # store the data as a csv file
# with open('pipolins.csv', 'w') as ouf:
#     print(','.join(PolB_Genomes._fields), file=ouf)
#     for i_g, _ in enumerate(polB_genomes):
#         pip = ['None' if i is None else str(i) for i in polB_genomes[i_g]]
#         print(','.join(pip), file=ouf, )


# # Create pipolin regions for annotation
# os.mkdir('pipolin_regions')
# for genome in polB_genomes:
#     genome_seq = get_seq_from_fasta(f'./genomes/{genome.id}.fa')
#
#     with open(f'./pipolin_regions/{genome.id}_selected.fa', 'w') as ouf:
#         print(f'>{genome.id} selected region {int(genome.attR_e) - int(genome.attL_s) + 100}', file=ouf)
#         print(f'{genome_seq[int(genome.attL_s) - 51:int(genome.attR_e) + 49]}', file=ouf)
#
# annotate pipolins with prokka
# os.mkdir('./prokka_annotations')
# for i_g, genome in enumerate(polB_genomes):
#     print(i_g)
#     subprocess.run(['/home/liubov/repos/prokka/bin/prokka', '--outdir', './prokka_annotations',
#                     '--prefix', f'{genome.id}', '--rfam', '--rawproduct', '--cdsrnaolap', '--cpus', '4',
#                     '--locustag', f'{genome.id}', f'./pipolin_regions/{genome.id}_selected.fa'])
#
