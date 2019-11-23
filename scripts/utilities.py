import os
import subprocess
from Bio import SearchIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def check_dir(out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
        try:
            len(os.listdir(out_dir)) == 0
        except RuntimeError:
            print('The specified output directory is not empty!')


def ncbi_acc_download(acc_ids):
    for acc_id in acc_ids:
        print(f'Downloading {acc_id}...')
        # TODO: fix local path!
        subprocess.run(['/home/liubov/.local/bin/ncbi-acc-download', '-F', 'fasta', acc_id])


def blast_genomes_against_seq(genomes_dir, seq, output_dir):
    check_dir(output_dir)
    genomes = os.listdir(genomes_dir)

    for genome in genomes:
        with open(os.path.join(output_dir, f'{genome[:-3]}_fmt7.txt'), 'w') as ouf:
            subprocess.run(['blastn', '-query', seq,
                            '-subject', f'{os.path.join(genomes_dir, genome)}',
                            '-outfmt', '7'], stdout=ouf)


def read_blasttab(blast_tab):
    return SearchIO.read(blast_tab, 'blast-tab', comments=True)


def save_as_csv(Pipolin, pipolins, out_file):
    with open(out_file, 'w') as ouf:
        print(','.join(Pipolin._fields), file=ouf)
        for i_p, _ in enumerate(pipolins):
            words = ['None' if i is None else str(i) for i in pipolins[i_p]]
            print(','.join(words), file=ouf)