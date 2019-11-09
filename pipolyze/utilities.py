import os
import subprocess
from Bio import SearchIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def ncbi_acc_download(ids):
    for genome_id in ids:
        # TODO fix local path!
        subprocess.run(['/home/liubov/.local/bin/ncbi-acc-download', '-F', 'fasta', genome_id])


def blast_genomes_against_seq(genomes_dir, seq):
    os.chdir(os.path.dirname(seq))
    if not os.path.exists('att_blast'):
        os.mkdir('att_blast')
    genomes = os.listdir(genomes_dir)

    for genome in genomes:
        with open(f'./att_blast/{genome[:-3]}_hits.txt', 'w') as ouf:
            subprocess.run(['blastn', '-query', seq,
                            '-subject', f'{os.path.join(genomes_dir, genome)}',
                            '-outfmt', '7'], stdout=ouf)


def get_hit_positions_by_id(blast_tab):
    blast_results = SearchIO.read(blast_tab, 'blast-tab', comments=True)
    hit_positions_by_id = {}

    for entry in blast_results:
        hit_positions = [sorted([hit.hit_start, hit.hit_end]) for hit in entry]
        hit_positions_by_id[entry.id] = sorted(hit_positions)

    return hit_positions_by_id
