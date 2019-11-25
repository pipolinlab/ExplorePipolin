import os
import subprocess
from Bio import SearchIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class Feature:
    def __init__(self, start, end, node):
        self.start = start
        self.end = end
        self.node = node

    def normalize(self):
        if self.start > self.end:
            return Feature(self.end, self.start, self.node)
        else:
            return self


class Pipolin:
    def __init__(self, strain_id):
        self.strain_id = strain_id
        self.polymerases = []
        self.atts = []

    def get_pipolin_bounds(self):
        if not self.is_complete_genome():
            raise AssertionError('Should be complete!')
        polymerases = sorted((i.normalize() for i in self.polymerases), key=lambda p: p.start)
        atts = sorted((i.normalize() for i in self.atts), key=lambda p: p.start)

        if not self._is_polymerase_inside(atts, polymerases):
            raise AssertionError('The polymerases are not within att bounds!')

        return atts[0].start, atts[-1].end

    def _is_polymerase_inside(self, atts, polymerases):
        return atts[0].start < polymerases[0].start and polymerases[0].end < atts[-1].end

    def get_contigs(self):
        contigs = set(i.node for i in self.polymerases)
        contigs.update(i.node for i in self.atts)
        return contigs

    def is_complete_genome(self):
        return self.strain_id == self.polymerases[0].node


def check_dir(out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


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
    # TODO: rewrite this!
    with open(out_file, 'w') as ouf:
        print(','.join(Pipolin._fields), file=ouf)
        for i_p, _ in enumerate(pipolins):
            words = ['None' if i is None else str(i) for i in pipolins[i_p]]
            print(','.join(words), file=ouf)