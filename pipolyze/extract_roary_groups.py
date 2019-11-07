import os
from Bio import SeqIO

os.getcwd()
# go to the working directory
os.chdir('./tfm/NCBI_strains_analysis')

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
        roary_groups[group[0]] = [group[i] for i in range(-67, 0) if len(group[i]) > 0]
# print(len(roary_groups))   # 233

faa_files = []
for file in os.listdir('./prokka_annotations'):
    if file.endswith('.faa'):
        faa_files.append(os.path.join('./prokka_annotations', file))
# print(len(faa_files))   # 67

faa_records = {}
for file in faa_files:
    faa_records.update(SeqIO.to_dict(SeqIO.parse(file, 'fasta')))

os.mkdir('./roary_groups')
for group_id, genes in roary_groups.items():
    seqs = [faa_records[gene] for gene in genes if gene in faa_records]
    with open(f'./roary_groups/{group_id}.faa', 'w') as ouf:
        SeqIO.write(seqs, ouf, 'fasta')
