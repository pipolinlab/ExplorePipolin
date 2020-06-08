import os
import subprocess

from explore_pipolin.common import define_genome_id


def run_prokka(gquery_id, pipolins_dir, proteins, prokka_dir):
    os.makedirs(prokka_dir, exist_ok=True)
    subprocess.run(['prokka', '--outdir', prokka_dir,
                    '--prefix', gquery_id,
                    # TODO: number of CPUs is hardcoded. To pass it as an argument?
                    '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--rfam', '--proteins', proteins, '--force',
                    '--locustag', gquery_id,
                    os.path.join(pipolins_dir, gquery_id + '.fa')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def blast_for_repeats(gquery_id, repeats_dir):
    with open(os.path.join(repeats_dir, gquery_id + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', os.path.join(repeats_dir, gquery_id + '.left'),
                        '-subject', os.path.join(repeats_dir, gquery_id + '.right'),
                        '-outfmt', '5', '-perc_identity', '100', '-word_size', '6',
                        '-strand', 'plus'], stdout=ouf)


def run_aragorn(genome, aragorn_results_dir):
    with open(os.path.join(aragorn_results_dir, define_genome_id(genome_path=genome) + '.batch'), 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome], stdout=ouf)


def tblastn_against_ref_pipolb(genome, ref_pipolb, out_dir):
    with open(os.path.join(out_dir, define_genome_id(genome_path=genome) + '.fmt5'), 'w') as ouf:
        subprocess.run(['tblastn', '-query', ref_pipolb, '-subject', genome, '-evalue', '0.1', '-outfmt', '5'],
                       stdout=ouf)


def blastn_against_ref_att(genome, ref_att, out_dir):
    with open(os.path.join(out_dir, define_genome_id(genome_path=genome) + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', ref_att, '-subject', genome, '-evalue', '0.1', '-outfmt', '5'],
                       stdout=ouf)
