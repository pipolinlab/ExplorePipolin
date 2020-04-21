#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import prefect
from prefect import task
from prefect.engine import signals
from pipolin_finder.utilities import CONTEXT_SETTINGS
from pipolin_finder.utilities import save_left_right_subsequences
from pipolin_finder.utilities import blast_for_identical
from pipolin_finder.utilities import extract_repeats
from pipolin_finder.utilities import set_proper_location
from pipolin_finder.utilities import read_aragorn_batch
from pipolin_finder.utilities import blast_genome_against_seq
from pipolin_finder.utilities import Orientation, Feature, Contig, GQuery
from pipolin_finder.utilities import read_blastxml
from pipolin_finder.utilities import run_aragorn
from pipolin_finder.utilities import read_seqio_records
from pipolin_finder.utilities import define_gquery_id


@task
def create_gquery(genome):
    gquery = GQuery(gquery_id=define_gquery_id(genome=genome))
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
            feature = gquery.feature_from_blasthit(hit=hit, contig_id=entry.id)
            gquery.get_features_by_type(feature_type).append(feature)


@task
def detect_trnas_with_aragorn(genome, root_dir):
    aragorn_results = os.path.join(root_dir, 'aragorn_results')
    run_aragorn(genome, aragorn_results)
    return aragorn_results


@task
def add_features_from_aragorn(gquery, aragorn_dir):
    entries = read_aragorn_batch(aragorn_batch=os.path.join(aragorn_dir, f'{gquery.gquery_id}.batch'))
    for contig_id, hits in entries.items():
        for hit in hits:
            feature = Feature(start=hit[0], end=hit[1], frame=hit[2], contig=gquery.get_contig_by_id(contig_id))
            gquery.trnas.append(feature)

    for att in gquery.atts:
        target_trna = gquery.define_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas.append(target_trna)


@task()
def find_atts_denovo(genome, gquery, root_dir):
    logger = prefect.context.get('logger')

    if not gquery.is_single_contig():
        logger.warning('This step is only for complete genomes. Pass...')
        raise signals.SKIP()

    repeats_dir = os.path.join(root_dir, 'atts_denovo')

    left_window, right_window = gquery.get_left_right_windows()
    save_left_right_subsequences(genome=genome, left_window=left_window, right_window=right_window,
                                 repeats_dir=repeats_dir)

    blast_for_identical(gquery_id=gquery.gquery_id, repeats_dir=repeats_dir)
    left_repeats, right_repeats = extract_repeats(file=os.path.join(repeats_dir, gquery.gquery_id + '.fmt5'))
    left_repeats = [set_proper_location(seq_range=i, shift=left_window[0]) for i in left_repeats]
    right_repeats = [set_proper_location(seq_range=i, shift=right_window[0]) for i in right_repeats]
    atts_denovo = [(i, j) for i, j in zip(left_repeats, right_repeats) if gquery.is_att_denovo(i, j)]

    with open(os.path.join(repeats_dir, gquery.gquery_id + '.atts'), 'w') as ouf:
        print('attL_start', 'attL_end', 'attR_start', 'attR_end', sep='\t', file=ouf)
        for att_pair in atts_denovo:
            print(att_pair[0][0], att_pair[0][1], att_pair[1][0], att_pair[1][1], sep='\t', file=ouf)

    return repeats_dir


@task
def add_features_atts_denovo(gquery, atts_denovo_dir):
    with open(os.path.join(atts_denovo_dir, gquery.gquery_id + '.atts')) as inf:
        _ = inf.readline()
        for line in inf:
            att_pair = line.strip().split(sep='\t')
            gquery.denovo_atts.append(Feature(start=att_pair[0], end=att_pair[1],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))
            gquery.denovo_atts.append(Feature(start=att_pair[2], end=att_pair[3],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))

    for att in gquery.denovo_atts:
        target_trna = gquery.define_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas_denovo.append(target_trna)


@task
def are_polbs_present(gquery):
    logger = prefect.context.get('logger')

    if len(gquery.polbs) == 0:
        logger.warning('No piPolB! => No pipolin!')
        raise signals.FAIL()


@task(skip_on_upstream_skip=False)
def are_atts_present(gquery):
    logger = prefect.context.get('logger')

    if len(gquery.atts) == 0 and len(gquery.denovo_atts) == 0:
        logger.warning('There is piPolB, but no atts were found! Not able to define pipolin bounds!')
        # TODO: probably, it makes sense to output piPolB(s) alone
        raise signals.SKIP()

    if len(gquery.atts) == 0:
        logger.warning(f'No "usual" atts were found, but some atts were found by denovo search!'
                       f'For more details, check the {gquery.gquery_id}.atts file in the atts_denovo directory!')
        # TODO: check that it's only one repeat!
        gquery.atts.extend(gquery.denovo_atts)
        gquery.target_trnas.extend(gquery.target_trnas_denovo)
    if len(gquery.denovo_atts) != 0:
        logger.warning(f'Some atts were found by denovo search, but we are not going to use them!'
                       f'For more details, check the {gquery.gquery_id}.atts file in the atts_denovo directory!')


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
