#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS
from utilities import blast_genomes_against_seq
from utilities import Feature, Contig, GQuery
from utilities import save_to_shelve
from utilities import read_blastxml
from utilities import run_aragorn
from utilities import read_seqio_records


def feature_from_blasthit(hit, entry_id):
    start = hit.hit_start if hit.hit_frame == 1 else hit.hit_end
    end = hit.hit_end if hit.hit_frame == 1 else hit.hit_start
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


def find_repeats_in_genome(genome_records, pipolin):
    if pipolin.is_complete_genome():
        pass


@task
def find_att_repeats(genomes, pipolins, root_dir):
    att_repeats_dir = os.path.join(root_dir, 'att_repeats')
    os.makedirs(att_repeats_dir, exist_ok=True)
    genomes_dict = read_seqio_records(files=genomes, file_format='fasta')
    for i_p, pipolin in enumerate(pipolins):
        with open(os.path.join(att_repeats_dir, f'{os.path.basename(pipolin.strain_id)[:-3]}.batch'), 'w') as ouf:
            repeats = find_repeats_in_genome(genome_records=genomes_dict[pipolin.strain_id], pipolin=pipolin)
            # for repeat in repeats:
            #     print(repeat, file=ouf)
    return att_repeats_dir


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
