#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import prefect
from prefect import task
from Bio import SeqIO
import subprocess
from collections import defaultdict
from utilities import CONTEXT_SETTINGS
from utilities import blast_genome_against_seq
from utilities import Orientation, Feature, Contig, GQuery
from utilities import save_to_shelve
from utilities import read_blastxml
from utilities import run_aragorn
from utilities import read_seqio_records


def feature_from_blasthit(hit, entry_id):
    return Feature(start=hit.hit_start, end=hit.hit_end,
                   frame=Orientation.orientation_from_blast(hit.hit_frame), contig=entry_id)


@task
def create_gquery(genome):
    gquery = GQuery(gquery_id=os.path.splitext(os.path.basename(genome))[0])
    genome_dict = read_seqio_records(file=genome, file_format='fasta')
    for key, value in genome_dict.items():
        contig = Contig(contig_id=key, contig_length=len(value.seq))
        gquery.contigs.append(contig)

    return gquery


@task
def run_blast_against_ref(genome, root_dir, reference, dir_name):
    blast_path = os.path.join(root_dir, dir_name)
    blast_genome_against_seq(genome=genome, seq=reference, output_dir=blast_path)
    return blast_path


@task
def add_features_from_blast(gquery, blast_dir, feature_type):
    entries = read_blastxml(blast_xml=os.path.join(blast_dir, f'{gquery.gquery_id}.fmt5'))
    for entry in entries:
        for hit in entry:
            feature = feature_from_blasthit(hit=hit, entry_id=entry.id)
            gquery.get_features_by_type(feature_type).append(feature)


@task
def detect_trnas(genome, root_dir):
    aragorn_results = os.path.join(root_dir, 'aragorn_results')
    run_aragorn(genome, aragorn_results)
    return aragorn_results


def read_aragorn_batch(aragorn_batch):
    entries = defaultdict(list)
    with open(aragorn_batch) as inf:
        for line in inf:
            if line[0] == '>':
                entry = line.split(sep=' ')[0][1:]
            else:
                hit = line.split(sep='\t')
                if len(hit) < 2:
                    continue
                else:
                    coordinates = hit[0].split(sep=' ')[-1]
                    if coordinates[0] == 'c':
                        start, end = (int(i) for i in coordinates[2:-1].split(sep=','))
                        entries[entry].append((start, end, Orientation.REVERSE))
                    else:
                        start, end = (int(i) for i in coordinates[1:-1].split(sep=','))
                        entries[entry].append((start, end, Orientation.FORWARD))

    return entries


@task
def add_features_from_aragorn(gquery, aragorn_dir):
    entries = read_aragorn_batch(aragorn_batch=os.path.join(aragorn_dir, f'{gquery.gquery_id}.batch'))
    for contig, hits in entries.items():
        for hit in hits:
            feature = Feature(start=hit[0], end=hit[1], frame=hit[2], contig=contig)
            gquery.trnas.append(feature)


def save_left_right_subsequences(genome_dict, gquery, repeats_dir):
    contig_id = next(iter(genome_dict.keys()))
    left_window, right_window = gquery.get_left_right_windows()
    left_seq = genome_dict[contig_id][left_window[0]:left_window[1]]
    SeqIO.write(sequences=left_seq, handle=os.path.join(repeats_dir, gquery.gquery_id + '.left'), format='fasta')
    right_seq = genome_dict[contig_id][right_window[0]:right_window[1]]
    SeqIO.write(sequences=right_seq, handle=os.path.join(repeats_dir, gquery.gquery_id + '.right'), format='fasta')


def blast_for_identical(gquery_id, repeats_dir):
    with open(os.path.join(repeats_dir, gquery_id + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', os.path.join(repeats_dir, gquery_id + '.left'),
                        '-subject', os.path.join(repeats_dir, gquery_id + '.right'),
                        '-outfmt', '5', '-perc_identity', '100', '-word_size', '12',
                        '-strand', 'plus'], stdout=ouf)


def get_repeats_location(repeats, gquery):
    left_window, right_window = gquery.get_left_right_windows()
    qrepeats_location = []
    srepeats_location = []
    for entry in repeats:
        for hit in entry:
            qrepeats_location.append((hit.query_start + left_window[0], hit.query_end + left_window[0]))
            srepeats_location.append((hit.hit_start + right_window[0], hit.hit_end + right_window[0]))

    return qrepeats_location, srepeats_location


def remove_overlapping_atts(gquery, qrepeats_location, srepeats_location):
    atts_location = [(i.start, i.end) for i in gquery.atts]
    leave_this = []
    for i, i_rep in enumerate(qrepeats_location):
        # TODO: FIX THIS!
        for i_att in atts_location:
            max_start = max(i_rep[0], i_att[0])
            min_end = min(i_rep[1], i_att[1])
            if max_start <= min_end:
                continue
            else:
                leave_this.append(i)
                break

    q_filtered = [qrepeats_location[i] for i in leave_this]
    s_filtered = [srepeats_location[i] for i in leave_this]
    return q_filtered, s_filtered


def leave_overlapping_trnas(gquery, qrepeats_location, srepeats_location):
    trnas_location = [(i.start, i.end) for i in gquery.trnas]
    leave_this = []
    for i, i_rep in enumerate(zip(qrepeats_location, srepeats_location)):
        for i_trna in trnas_location:
            q_max_start = max(i_rep[0][0], i_trna[0])
            q_min_end = min(i_rep[0][1], i_trna[1])
            s_max_start = max(i_rep[1][0], i_trna[0])
            s_min_end = min(i_rep[1][1], i_trna[1])
            if q_max_start <= q_min_end or s_max_start <= s_min_end:
                leave_this.append(i)

    q_filtered = [qrepeats_location[i] for i in leave_this]
    s_filtered = [srepeats_location[i] for i in leave_this]
    return q_filtered, s_filtered


@task
def find_atts_denovo(genome, gquery, root_dir):
    logger = prefect.context.get('logger')

    if not gquery.is_single_contig():
        logger.warning('This step is only for complete genomes. Pass...')
        return

    repeats_dir = os.path.join(root_dir, 'repeats_denovo')
    os.makedirs(repeats_dir, exist_ok=True)

    genome_dict = read_seqio_records(file=genome, file_format='fasta')
    save_left_right_subsequences(genome_dict=genome_dict, gquery=gquery, repeats_dir=repeats_dir)

    blast_for_identical(gquery_id=gquery.gquery_id, repeats_dir=repeats_dir)
    repeats = read_blastxml(os.path.join(repeats_dir, gquery.gquery_id + '.fmt5'))
    qrepeats_location, srepeats_location = get_repeats_location(repeats=repeats, gquery=gquery)
    q_filtered1, s_filtered1 = remove_overlapping_atts(gquery=gquery,
                                                       qrepeats_location=qrepeats_location,
                                                       srepeats_location=srepeats_location)
    q_filtered, s_filtered = leave_overlapping_trnas(gquery=gquery,
                                                     qrepeats_location=q_filtered1,
                                                     srepeats_location=s_filtered1)
    with open(os.path.join(repeats_dir, gquery.gquery_id + '.loc'), 'w') as ouf:
        print('lrepeat_loc', 'rrepeat_loc', sep='\t', file=ouf)
        for q, s in zip(q_filtered, s_filtered):
            print(q, s, sep='\t', file=ouf)


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
