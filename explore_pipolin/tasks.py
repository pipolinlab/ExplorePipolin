"""
All the tasks should be defined in this file!
"""
import os

import pkg_resources
from prefect import task
from prefect import context

from Bio import SeqIO
from typing import Any, Optional, Sequence

from explore_pipolin.tasks_related.easyfig_coloring import add_colours
from explore_pipolin.common import Feature, FeatureType, Pipolin, Genome, \
    define_genome_id, FeaturesContainer, Strand, Range, ContigID
from explore_pipolin.utilities.logging import genome_specific_logging
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
def save_pipolin_sequences(genome: Genome, pipolins: Sequence[Pipolin], out_dir: str):
    pipolins_dir = os.path.join(out_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    genome_dict = create_seqio_records_dict(file=genome.file, file_format='fasta')

    for i, pipolin in enumerate(pipolins):
        with open(os.path.join(pipolins_dir, genome.id + f'_{i}.fa'), 'w') as ouf:
            records = [create_fragment_record(f, genome_dict) for f in pipolin.fragments]
            SeqIO.write(records, ouf, 'fasta')

    return pipolins_dir


@task()
@genome_specific_logging
def annotate_pipolins(genome: Genome, pipolins_dir, out_dir):
    prokka_results_dir = os.path.join(out_dir, 'prokka_results')
    os.makedirs(prokka_results_dir, exist_ok=True)

    for fasta_file in os.listdir(pipolins_dir):
        if fasta_file.startswith(genome.id):
            prefix = os.path.splitext(fasta_file)
            input_file = os.path.join(pipolins_dir, fasta_file)
            run_prokka(prefix=prefix, input_file=input_file, prokka_results_dir=prokka_results_dir)

    return prokka_results_dir


@task()
@genome_specific_logging
def include_atts(genome: Genome, prokka_dir, out_dir, pipolins: Sequence[Pipolin]):
    results_dir = os.path.join(out_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    for prokka_file in os.listdir(prokka_dir):

        pipolin_index = int(os.path.splitext(prokka_file)[0].split(sep='_')[-1])
        cur_pipolin = pipolins[pipolin_index]

        if prokka_file.startswith(genome.id) and prokka_file.endswith('.gbk'):
            gb_records = create_seqio_records_dict(file=os.path.join(prokka_dir, prokka_file),
                                                   file_format='genbank')

            include_atts_into_gb(gb_records=gb_records, pipolin=cur_pipolin)

            output_file = os.path.join(results_dir, f'{prokka_file}.gbk')
            write_genbank_records(gb_records=gb_records, output_file=output_file)

        if prokka_file.startswith(genome.id) and prokka_file.endswith('.gff'):
            gff_records = read_gff_records(file=os.path.join(prokka_dir, prokka_file))
            include_atts_into_gff(gff_records=gff_records, pipolin=cur_pipolin)

            output_file = os.path.join(results_dir, f'{prokka_file}.gff')
            write_gff_records(gff_records=gff_records, output_file=output_file)

    return results_dir


@task()
@genome_specific_logging
def easyfig_add_colours(genome: Genome, in_dir):
    for gbk_file in os.listdir(in_dir):
        if gbk_file.startswith(genome.id):
            gb_records = create_seqio_records_dict(file=os.path.join(in_dir, gbk_file),
                                                   file_format='genbank')
            for record in gb_records.values():
                add_colours(record)

            write_genbank_records(gb_records=gb_records, output_file=os.path.join(in_dir, gbk_file))
