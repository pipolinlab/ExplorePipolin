"""
All the tasks should be defined in this file!
"""
import os

import pkg_resources
from prefect import task
from prefect import context

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import Any, Optional, Sequence

from explore_pipolin.tasks_related.easyfig_coloring import add_colours
from explore_pipolin.common import Feature, FeatureType, Pipolin, Genome, \
    define_genome_id, FeaturesContainer, Strand, Range, ContigID
from explore_pipolin.utilities.logging import genome_specific_logging
from explore_pipolin.tasks_related.misc import join_it
from explore_pipolin.utilities.io import create_pipolb_entries
from explore_pipolin.utilities.io import create_seqio_records_dict
from explore_pipolin.utilities.io import read_aragorn_batch
from explore_pipolin.utilities.external_tools import ExternalTools, RealExternalTools
from explore_pipolin.utilities.external_tools import run_prokka, run_aragorn
from explore_pipolin.tasks_related.misc import create_fragment_record
from explore_pipolin.utilities.io import read_gff_records
from explore_pipolin.tasks_related.including_atts import include_atts_into_gb
from explore_pipolin.tasks_related.including_atts import include_atts_into_gff
from explore_pipolin.utilities.io import write_genbank_records
from explore_pipolin.utilities.io import write_gff_records
from explore_pipolin.utilities.io import read_genome_contigs_from_file

_REF_PIPOLB = pkg_resources.resource_filename('explore_pipolin', 'data/pi-polB_Firmicutes.faa')


# Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
# https://github.com/biopython/biopython/issues/1755


@task()
def create_genome(genome_file) -> Genome:
    contigs = read_genome_contigs_from_file(genome_file=genome_file)
    genome_id = define_genome_id(genome_file)
    return Genome(genome_id=genome_id, genome_file=genome_file, contigs=contigs)


@task()
@genome_specific_logging
def find_pipolbs(genome: Genome, out_dir: str, ext: ExternalTools = RealExternalTools) -> Genome:
    results_dir = os.path.join(out_dir, 'pipolbs_search')
    os.makedirs(results_dir, exist_ok=True)

    prodigal_output_file = os.path.join(results_dir, genome.id + '.faa')
    ext.find_cdss(genome_file=genome.file, output_file=prodigal_output_file)
    hmm_output_file = os.path.join(results_dir, genome.id + '.tbl')
    ext.find_pipolbs(prodigal_output_file, hmm_output_file)
    entries_pipolb = create_pipolb_entries(hmmsearch_table=hmm_output_file)

    for entry in entries_pipolb:
        feature = Feature(location=Range(start=entry[1], end=entry[2]),
                          strand=Strand.from_pm_one_encoding(entry[3]),
                          contig_id=ContigID(entry[0]), genome=genome)
        genome.features.add_features(feature, feature_type=FeatureType.PIPOLB)

    return genome


@task()
@genome_specific_logging
def find_trnas(genome: Genome, out_dir) -> Genome:
    aragorn_results_dir = os.path.join(out_dir, 'trnas_search')
    os.makedirs(aragorn_results_dir, exist_ok=True)

    output_file = os.path.join(aragorn_results_dir, genome.id + '.batch')
    run_aragorn(genome_file=genome.file, output_file=output_file)
    entries = read_aragorn_batch(aragorn_batch=output_file)

    add_trna_features_from_aragorn_entries(entries=entries, genome=genome)
    find_and_add_target_trnas_features(genome.features)

    return genome


def add_trna_features_from_aragorn_entries(entries, genome: Genome):
    for contig_id, hits in entries.items():
        for hit in hits:
            # "correct strange coordinates in -l mode" as in Prokka
            start = max(hit[0], 1)
            end = min(hit[1], genome.get_contig_by_id(contig_id=contig_id).length)
            trna_feature = Feature(location=Range(start=start, end=end),
                                   strand=hit[2], contig_id=contig_id, genome=genome)
            genome.features.add_features(trna_feature, feature_type=FeatureType.TRNA)


def find_and_add_target_trnas_features(features: FeaturesContainer):
    for att in features.get_features(FeatureType.ATT):
        target_trna = features.get_features(FeatureType.TRNA).get_overlapping(att)
        if target_trna is not None:
            features.add_features(target_trna, feature_type=FeatureType.TARGET_TRNA)


@task()
@genome_specific_logging
def are_pipolbs_present(genome: Genome):
    logger = context.get('logger')

    if len(genome.features.get_features(FeatureType.PIPOLB)) == 0:
        logger.warning('No piPolBs were found!')
        return False

    return True


@task()
def return_result_if_true_else_none(result_to_filter: Any, filter_by: bool) -> Optional[Any]:
    if filter_by:
        return result_to_filter

    return None


@task()
@genome_specific_logging
def extract_pipolin_regions(genome: Genome, pipolins: Sequence[Pipolin], out_dir: str):
    genome_dict = create_seqio_records_dict(file=genome.file, file_format='fasta')

    pipolins_dir = os.path.join(out_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    logger = context.get('logger')

    with open(os.path.join(pipolins_dir, genome.id + '.fa'), 'w') as ouf:
        fragment_records = [create_fragment_record(fragment=f, genome_dict=genome_dict) for f in pipolins.fragments]

        for fragment_record, fragment in zip(fragment_records, pipolins.fragments):
            logger.info(f'@pipolin fragment length {len(fragment_record)} from {fragment.contig_id}')

        record = sum(join_it(fragment_records, SeqRecord(seq='N' * 100)), SeqRecord(seq=''))

        logger.info(f'@@@pipolin record total length {len(record)}')

        record.id = genome.id
        record.name = genome.id
        record.description = genome.id
        SeqIO.write(sequences=record, handle=ouf, format='fasta')

    return pipolins_dir


@task()
@genome_specific_logging
def annotate_pipolins(genome: Genome, pipolins_dir, out_dir):
    prokka_results_dir = os.path.join(out_dir, 'prokka_results')
    os.makedirs(prokka_results_dir, exist_ok=True)
    run_prokka(genome_id=genome.id, pipolins_dir=pipolins_dir, prokka_results_dir=prokka_results_dir)
    return prokka_results_dir


@task()
@genome_specific_logging
def include_atts(genome: Genome, prokka_dir, out_dir, pipolins: Sequence[Pipolin]):
    gb_records = create_seqio_records_dict(file=os.path.join(prokka_dir, genome.id + '.gbk'),
                                           file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, genome.id + '.gff'))

    include_atts_into_gb(gb_records=gb_records, genome=genome, pipolins=pipolins)
    include_atts_into_gff(gff_records=gff_records, genome=genome, pipolins=pipolins)

    prokka_atts_dir = os.path.join(out_dir, 'results')
    os.makedirs(prokka_atts_dir, exist_ok=True)

    write_genbank_records(gb_records=gb_records, out_dir=prokka_atts_dir, genome=genome)
    write_gff_records(gff_records=gff_records, out_dir=prokka_atts_dir, genome=genome)

    return prokka_atts_dir


@task()
@genome_specific_logging
def easyfig_add_colours(genome: Genome, in_dir):
    gb_records = create_seqio_records_dict(file=os.path.join(in_dir, genome.id + '.gbk'),
                                           file_format='genbank')
    add_colours(gb_records[genome.id])
    write_genbank_records(gb_records=gb_records, out_dir=in_dir, genome=genome)
