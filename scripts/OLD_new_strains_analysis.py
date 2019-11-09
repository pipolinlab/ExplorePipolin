import os
# import shutil
from collections import namedtuple
# import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

print(os.getcwd())
# go to the working directory
os.chdir('tfm/new_strains_analysis')


def parse_blast_hittable(ht_file, id=True):
    # :return: list of tuples with subject id, start and end positions
    with open(ht_file) as inf:
        entries = []
        for line in inf:
            if line[0] != '#':
                entry = line.strip().split(sep='\t')
                if id:
                    entries.append((entry[1], entry[8], entry[9]))
                else:
                    entries.append((entry[8], entry[9]))
    return entries


# read contigs from fasta files:
files_w_contigs = []
strain_names = []
for file in os.listdir('./contigs'):
    strain_names.append(file.split(sep='.')[0])
    files_w_contigs.append(os.path.join('./contigs', file))

contigs_by_strain = []
for file in files_w_contigs:
    contigs_by_strain.append(SeqIO.to_dict(SeqIO.parse(file, 'fasta')))

# create a namedtuple object to store info
fields = ('strain_id',
          'polB_1_s', 'polB_1_e', 'polB_1_node',
          'polB_2_s', 'polB_2_e', 'polB_2_node',
          'att1_s', 'att1_e', 'att1_node',
          'att2_s', 'att2_e', 'att2_node',
          'att3_s', 'att3_e', 'att3_node')
NewStrains = namedtuple('NewStrains', fields, defaults=(None, ) * len(fields))

# # run BLAST against ref_att.fa
# os.mkdir('att_blast')
# for strain in strain_names:
#     with open(f'./att_blast/{strain}_hits.txt', 'w') as ouf:
#         subprocess.run(['blastn', '-query', f'../NCBI_strains_analysis/ref_att.fa',
#                         '-subject', f'./contigs/{strain}.fa', '-outfmt', '6'],
#                        stdout=ouf)
#
# # run BLAST against ref pi-polB
# os.mkdir('polB_blast')
# for strain in strain_names:
#     with open(f'./polB_blast/{strain}_hits.txt', 'w') as ouf:
#         subprocess.run(['blastn', '-query', f'../piPolBnt.fa',
#                         '-subject', f'./contigs/{strain}.fa', '-outfmt', '6'],
#                        stdout=ouf)

# parse att BLAST output
new_strains = []
for i_s, strain in enumerate(strain_names):
    att_regions = parse_blast_hittable(f'./att_blast/{strain}_hits.txt')
    att_regions.sort(key=lambda i: i[1:])

    polB_regions = parse_blast_hittable(f'./polB_blast/{strain}_hits.txt')

    new_strains.append(NewStrains(strain_id=strain,
                                  polB_1_s=polB_regions[0][1], polB_1_e=polB_regions[0][2],
                                  polB_1_node=polB_regions[0][0],
                                  att1_s=att_regions[0][1], att1_e=att_regions[0][2],
                                  att1_node=att_regions[0][0]))
    if len(polB_regions) == 2:
        new_strains[i_s] = new_strains[i_s]._replace(polB_2_s=polB_regions[1][1],
                                                     polB_2_e=polB_regions[1][2],
                                                     polB_2_node=polB_regions[1][0])
    if len(att_regions) > 1:
        new_strains[i_s] = new_strains[i_s]._replace(att2_s=att_regions[1][1],
                                                     att2_e=att_regions[1][2],
                                                     att2_node=att_regions[1][0])
        if len(att_regions) == 3:
            new_strains[i_s] = new_strains[i_s]._replace(att3_s=att_regions[2][1],
                                                         att3_e=att_regions[2][2],
                                                         att3_node=att_regions[2][0])

# # store the data as a csv file
# with open('new_pipolins.csv', 'w') as ouf:
#     print(','.join(NewStrains._fields), file=ouf)
#     for i_s, _ in enumerate(new_strains):
#         pip = ['None' if i is None else str(i) for i in new_strains[i_s]]
#         print(','.join(pip), file=ouf, )


os.mkdir('./pipolin_regions')

for i_s, strain in enumerate(strain_names):
    if new_strains[i_s].att1_node == new_strains[i_s].att2_node:
        node_seq = contigs_by_strain[i_s][new_strains[i_s].att1_node]
        pipoline_ends = (sorted([new_strains[i_s].att1_s, new_strains[i_s].att1_e])[0],
                       sorted([new_strains[i_s].att2_s, new_strains[i_s].att2_e])[0])
        with open(f'./pipolin_regions/{strain}_selected.fa', 'w') as ouf:
            record = SeqRecord(seq=node_seq.seq[int(pipoline_ends[0]) - 51:int(pipoline_ends[1]) + 49],
                               id=node_seq.id,
                               name=node_seq.name,
                               description=node_seq.description)
            SeqIO.write(record, ouf, 'fasta')

    elif new_strains[i_s].att2_node is None:
        node_seq = contigs_by_strain[i_s][new_strains[i_s].att1_node]
        end1 = 0 if int(new_strains[i_s].polB_1_s) < int(new_strains[i_s].att1_s) else len(node_seq) - 1
        single_att = sorted([int(new_strains[i_s].att1_s), int(new_strains[i_s].att1_e)])
        end2 = single_att[0] if single_att[0] < end1 else single_att[1]
        pipoline_ends = sorted([end1, end2])
        with open(f'./pipolin_regions/{strain}_selected.fa', 'w') as ouf:
            record = SeqRecord(seq=node_seq.seq[pipoline_ends[0] - 51:pipoline_ends[1] + 49],
                               id=node_seq.id, name=node_seq.name,
                               description=node_seq.description)
            SeqIO.write(record, ouf, 'fasta')
    # TODO rewrite this with better structure!

