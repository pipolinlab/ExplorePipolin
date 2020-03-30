#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
from Bio import SeqIO
import subprocess
from itertools import combinations
from utilities import CONTEXT_SETTINGS
from utilities import blast_genomes_against_seq
from utilities import Feature, Contig, GQuery
from utilities import save_to_shelve
from utilities import read_blastxml
from utilities import run_aragorn
from utilities import read_seqio_records


def feature_from_blasthit(hit, entry_id):
    start = hit.hit_start
    end = hit.hit_end
    return Feature(start=start, end=end, frame=hit.hit_frame, contig=entry_id)


@task
def create_gqueries(genomes):
    gqueries = []
    genomes_dict = read_seqio_records(files=genomes, file_format='fasta')
    for f_key, f_value in genomes_dict.items():
        gqueries.append(GQuery(gquery_id=f_key))
        for s_key, s_value in f_value.items():
            contig = Contig(contig_id=s_key, contig_length=len(s_value.seq))
            gqueries[-1].contigs.append(contig)

    return gqueries


@task
def add_polb_features(gqueries, polbs_blast_dir):
    for i_q, gquery in enumerate(gqueries):
        polbs = read_blastxml(os.path.join(polbs_blast_dir, f'{gquery.gquery_id}.fmt5'))
        for entry in polbs:
            for hit in entry:
                polb_feature = feature_from_blasthit(hit=hit, entry_id=entry.id)
                gqueries[i_q].polymerases.append(polb_feature)

    return gqueries


@task
def run_blast_against_polb(genomes, root_dir, ref_polb):
    polbs_blast_path = os.path.join(root_dir, 'polb_blast')
    blast_genomes_against_seq(genomes, ref_polb, polbs_blast_path)
    return polbs_blast_path


@task
def detect_trnas(genomes, out_dir):
    aragorn_results = os.path.join(out_dir, 'aragorn_results')
    run_aragorn(genomes, aragorn_results)
    return aragorn_results


def get_kmer_set(kmer_size, seq):
    kmer_set = set()
    for i in range(0, len(seq) - kmer_size):
        kmer_set.add(seq[i:i + kmer_size])

    return kmer_set


def remove_overlapping(repeated_kmers):
    kmer_pairs = list(combinations(repeated_kmers, 2))
    return repeated_kmers


def find_repeated_kmers(left_seq, right_seq):
    repeated_kmers = []
    for kmer_size in range(80, 11, -1):
        left_set = get_kmer_set(kmer_size=kmer_size, seq=left_seq)
        right_set = get_kmer_set(kmer_size=kmer_size, seq=right_seq)
        repeated_kmers.append(left_set.intersection(right_set))
    
    return remove_overlapping(repeated_kmers=repeated_kmers)


def save_left_right_subsequences(genome_contigs_dict, gquery, repeats_dir):
    if gquery.is_single_contig():
        contig_id = next(iter(genome_contigs_dict.keys()))
        left_window, right_window = gquery.get_left_right_windows()
        left_seq = genome_contigs_dict[contig_id][left_window[0]:left_window[1]]
        SeqIO.write(sequences=left_seq, handle=os.path.join(repeats_dir, gquery.gquery_id + '.left'), format='fasta')
        right_seq = genome_contigs_dict[contig_id][right_window[0]:right_window[1]]
        SeqIO.write(sequences=right_seq, handle=os.path.join(repeats_dir, gquery.gquery_id + '.right'), format='fasta')


def blast_for_identical(gquery_id, repeats_dir):
    with open(os.path.join(repeats_dir, gquery_id + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', os.path.join(repeats_dir, gquery_id + '.left'),
                        '-subject', os.path.join(repeats_dir, gquery_id + '.right'),
                        '-outfmt', '5', '-perc_identity', '100', '-word_size', '12', '-strand', 'plus'], stdout=ouf)


def filter_repeats(repeats):
    for entry in repeats:
        for hit in entry:
            print(hit)


@task
def find_atts_denovo(genomes, gqueries, root_dir) -> str:
    repeats_dir = os.path.join(root_dir, 'repeats')
    os.makedirs(repeats_dir, exist_ok=True)
    genomes_dict = read_seqio_records(files=genomes, file_format='fasta')
    for i_q, gquery in enumerate(gqueries):
        save_left_right_subsequences(genome_contigs_dict=genomes_dict[gquery.gquery_id],
                                     gquery=gquery, repeats_dir=repeats_dir)
        blast_for_identical(gquery_id=gquery.gquery_id, repeats_dir=repeats_dir)
        repeats = read_blastxml(os.path.join(repeats_dir, gquery.gquery_id + '.fmt5'))
        repeats_filtered = filter_repeats(repeats=repeats)

    return repeats_dir


@task   # TODO: replace with de-novo atts search finction!
def run_blast_against_att(genomes, out_dir, ref_att):
    atts_blast_path = os.path.join(out_dir, 'att_blast')
    blast_genomes_against_seq(genomes, ref_att, atts_blast_path)
    return atts_blast_path


@task   # TODO: replace with aragorn search function!
def run_blast_against_trna(genomes, out_dir, ref_trna):
    trna_blast_path = os.path.join(out_dir, 'trna_blast')
    blast_genomes_against_seq(genomes, ref_trna, trna_blast_path)
    return trna_blast_path


@task
def identify_pipolins_roughly(genomes, out_dir, polbs_blast, atts_blast):
    out_file = os.path.join(out_dir, 'shelve.db')
    save_to_shelve(out_file, pipolins, 'pipolins')
    return 'pipolins'


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ref-polb', type=click.Path(exists=True))
@click.argument('ref-att', type=click.Path(exists=True))
@click.argument('ref-trna', type=click.Path(exists=True))
@click.argument('genomes', nargs=-1, type=click.Path(exists=True))
@click.argument('out-dir')
def main(ref_polb, ref_att, ref_trna, genomes, out_dir):
    """
    GENOMES_DIR contains each genome in a separate FASTA file (strain_id.fa).
    If there are several contigs in the genome, each contig should have unique name.
    If OUT_DIR exists, it should be empty.
    """
    identify_pipolins_roughly(genomes, out_dir, ref_polb, ref_att, ref_trna)


if __name__ == '__main__':
    main()
