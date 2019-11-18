# TODO: rewrite the script and rerun the analysis!!!

import os
from Bio import SeqIO

os.getcwd()
# go to the working directory
os.chdir('/home/liubov/Documents/tfm/new_strains_analysis')

roary_groups = {}
with open('./roary_output/gene_presence_absence.csv') as inf:
    # in this file:
    #   - the first column is a common gene name or group id
    #   - the last 67 columns are presence/absence of this gene for each genome
    #     (if present -- then unique gene id, if absent -- empty)
    _ = inf.readline()   # skip a header line
    for line in inf:
        group = line.strip().replace('\t', ',').split(sep=',')
        group = [i.replace('"', '') for i in group]   # remove quotation marks from items
        roary_groups[group[0]] = [group[i] for i in range(-24, 0) if len(group[i]) > 0]
# print(len(roary_groups))   # 233

faa_files = []
for file in os.listdir('./prokka_annotations'):
    if file.endswith('.faa'):
        faa_files.append(os.path.join('./prokka_annotations', file))
# print(len(faa_files))   # 67

faa_records = {}
for file in faa_files:
    faa_records.update(SeqIO.to_dict(SeqIO.parse(file, 'fasta')))

os.mkdir('./roary_groups_aa')
for group_id, genes in roary_groups.items():
    seqs = [faa_records[gene] for gene in genes if gene in faa_records]
    with open(f'./roary_groups_aa/{group_id}.faa', 'w') as ouf:
        SeqIO.write(seqs, ouf, 'fasta')


ffn_files = []
for file in os.listdir('./prokka_annotations'):
    if file.endswith('.ffn'):
        ffn_files.append(os.path.join('./prokka_annotations', file))
# print(len(ffn_files))   # 67

ffn_records = {}
for file in ffn_files:
    ffn_records.update(SeqIO.to_dict(SeqIO.parse(file, 'fasta')))

os.mkdir('./roary_groups_nucl')
for group_id, genes in roary_groups.items():
    seqs = [ffn_records[gene] for gene in genes if gene in ffn_records]
    with open(f'./roary_groups_nucl/{group_id}.ffn', 'w') as ouf:
        SeqIO.write(seqs, ouf, 'fasta')
