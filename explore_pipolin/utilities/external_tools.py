import os
import subprocess

from explore_pipolin.utilities import define_gquery_id


def run_prokka(gquery_id, pipolins_dir, proteins, prokka_dir):
    os.makedirs(prokka_dir, exist_ok=True)
    subprocess.run(['prokka', '--outdir', prokka_dir,
                    '--prefix', gquery_id,
                    '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--rfam', '--proteins', proteins, '--force',
                    '--locustag', gquery_id,
                    os.path.join(pipolins_dir, gquery_id + '.fa')])


def blast_for_repeats(gquery_id, repeats_dir):
    with open(os.path.join(repeats_dir, gquery_id + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', os.path.join(repeats_dir, gquery_id + '.left'),
                        '-subject', os.path.join(repeats_dir, gquery_id + '.right'),
                        '-outfmt', '5', '-perc_identity', '100', '-word_size', '6',
                        '-strand', 'plus'], stdout=ouf)


def run_aragorn(genome, aragorn_results):
    os.makedirs(aragorn_results, exist_ok=True)
    with open(os.path.join(aragorn_results, define_gquery_id(genome=genome) + '.batch'), 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome], stdout=ouf)


def blast_genome_against_seq(genome, seq, seq_type, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, define_gquery_id(genome=genome) + '.fmt5'), 'w') as ouf:
        if seq_type == 'protein':
            subprocess.run(['tblastn', '-query', seq, '-subject', genome, '-evalue', '0.1', '-outfmt', '5'],
                           stdout=ouf)
        elif seq_type == 'nucleotide':
            subprocess.run(['blastn', '-query', seq, '-subject', genome, '-evalue', '0.1', '-outfmt', '5'],
                           stdout=ouf)
