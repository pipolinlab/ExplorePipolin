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
from explore_pipolin.common import Feature, FeatureType, RepeatPair, Pipolin, Genome, \
    define_genome_id, FeaturesContainer, Strand, Range
from explore_pipolin.utilities.logging import genome_specific_logging
from explore_pipolin.tasks_related.misc import join_it, add_features_from_blast_entries
from explore_pipolin.utilities.io import read_blastxml, write_repeats, write_atts_denovo, create_pipolb_entries
from explore_pipolin.utilities.io import create_seqio_records_dict
from explore_pipolin.utilities.io import read_aragorn_batch
from explore_pipolin.tasks_related.atts_denovo_search import find_repeats, is_att_denovo
from explore_pipolin.utilities.external_tools import ExternalTools, RealExternalTools
from explore_pipolin.utilities.external_tools import run_prokka, run_aragorn
from explore_pipolin.tasks_related.misc import create_fragment_record
from explore_pipolin.utilities.io import read_gff_records
from explore_pipolin.tasks_related.including_atts import include_atts_into_gb
from explore_pipolin.tasks_related.including_atts import include_atts_into_gff
from explore_pipolin.utilities.io import write_genbank_records
from explore_pipolin.utilities.io import write_gff_records
from explore_pipolin.utilities.io import read_genome_contigs_from_file
from explore_pipolin.tasks_related.scaffolding import Scaffolder, create_pipolin_fragments_single_contig

_REF_PIPOLB = pkg_resources.resource_filename('explore_pipolin', 'data/pi-polB_Firmicutes.faa')


@task
def create_genome(genome_file) -> Genome:
    contigs = read_genome_contigs_from_file(genome_file=genome_file)
    genome_id = define_genome_id(genome_file)
    return Genome(genome_id=genome_id, genome_file=genome_file, contigs=contigs)


@task
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
                          contig_id=entry[0], genome=genome)
        genome.features.add_feature(feature=feature, feature_type=FeatureType.PIPOLB)

    return genome


@task
@genome_specific_logging
def find_atts(genome: Genome, out_dir: str, ext: ExternalTools = RealExternalTools) -> Genome:
    blast_results_dir = os.path.join(out_dir, 'atts_search')
    os.makedirs(blast_results_dir, exist_ok=True)

    output_file = os.path.join(blast_results_dir, genome.id) + '.fmt5'
    ext.blastn_against_ref_att(genome_file=genome.file, output_file=output_file)
    entries = read_blastxml(blast_xml=output_file)

    add_features_from_blast_entries(entries=entries, feature_type=FeatureType.ATT, genome=genome)

    return genome


@task
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
            genome.features.add_feature(feature=trna_feature, feature_type=FeatureType.TRNA)


def find_and_add_target_trnas_features(features: FeaturesContainer):
    for att in features.get_features(FeatureType.ATT):
        target_trna = features.get_features(FeatureType.TRNA).get_overlapping(att)
        if target_trna is not None:
            features.add_feature(feature=target_trna, feature_type=FeatureType.TARGET_TRNA)


@task()
@genome_specific_logging
def are_pipolbs_present(genome: Genome):
    logger = context.get('logger')

    if len(genome.features.get_features(FeatureType.PIPOLB)) == 0:
        logger.warning('No piPolBs were found!')
        return False

    return True


@task
def return_result_if_true_else_none(result_to_filter: Any, filter_by: bool) -> Optional[Any]:
    if filter_by:
        return result_to_filter

    return None


@task()
@genome_specific_logging
def find_atts_denovo(genome: Genome, out_dir) -> Genome:
    atts_denovo_dir = os.path.join(out_dir, 'atts_denovo_search')
    os.makedirs(atts_denovo_dir, exist_ok=True)

    repeats: Sequence[RepeatPair] = find_repeats(genome=genome, repeats_dir=atts_denovo_dir)
    write_repeats(genome=genome, repeats=repeats, out_dir=atts_denovo_dir)

    atts_denovo: Sequence[RepeatPair] = [rep for rep in repeats if is_att_denovo(genome, rep)]
    for atts_pair in atts_denovo:
        genome.features.add_feature(feature=Feature(
            location=atts_pair.left_range,
            strand=Strand.FORWARD,
            contig_id=genome.contigs[0].id,
            genome=genome), feature_type=FeatureType.ATT_DENOVO)
        genome.features.add_feature(feature=Feature(
            location=atts_pair.right_range,
            strand=Strand.FORWARD,
            contig_id=genome.contigs[0].id,
            genome=genome), feature_type=FeatureType.ATT_DENOVO)

    for att in genome.features.get_features(FeatureType.ATT_DENOVO):
        target_trna = genome.features.get_features(FeatureType.TRNA).get_overlapping(att)
        if target_trna is not None:
            genome.features.add_feature(feature=target_trna, feature_type=FeatureType.TARGET_TRNA_DENOVO)

    write_atts_denovo(genome.features.get_features(FeatureType.ATT_DENOVO), genome, atts_denovo_dir)

    return genome


@task()
@genome_specific_logging
def are_atts_present(genome: Genome) -> Genome:
    logger = context.get('logger')

    num_atts = len(genome.features.get_features(FeatureType.ATT))
    num_atts_denovo = len(genome.features.get_features(FeatureType.ATT_DENOVO))
    if num_atts == 0 and num_atts_denovo == 0:
        logger.warning('\n\n>>>There is piPolB, but no atts were found! Not able to define pipolin bounds!\n')
        # TODO: probably, it makes sense to output_ecoli piPolB(s) alone

    elif len(genome.features.get_features(FeatureType.ATT)) == 0:
        logger.warning(f'\n\n>>>No "usual" atts were found, but some atts were found by denovo search!'
                       f'For more details, check the {genome.id}.atts_denovo file '
                       f'in the atts_denovo_search directory!\n')
        # TODO: check that it's only one repeat! Although, this shouldn't be a problem.
        atts_frames = [att.strand for att in genome.features.get_features(FeatureType.ATT_DENOVO)]
        if len(set(atts_frames)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_frames).pop() == genome.features.get_features(FeatureType.TARGET_TRNA_DENOVO).first.strand:
            reverse_denovo_atts = []
            for att in genome.features.get_features(FeatureType.ATT_DENOVO):
                reverse_denovo_atts.append(Feature(location=Range(start=att.start, end=att.end),
                                                   strand=-att.strand, contig_id=att.contig_id, genome=genome))
            [genome.features.add_feature(feature=i, feature_type=FeatureType.ATT) for i in reverse_denovo_atts]
        else:
            atts_denovo = genome.features.get_features(FeatureType.ATT_DENOVO)
            [genome.features.add_feature(feature=i, feature_type=FeatureType.ATT) for i in atts_denovo]
        target_trnas_denovo = genome.features.get_features(FeatureType.TARGET_TRNA_DENOVO)
        [genome.features.add_feature(feature=i, feature_type=FeatureType.TARGET_TRNA) for i in target_trnas_denovo]

    elif len(genome.features.get_features(FeatureType.ATT_DENOVO)) != 0:
        logger.warning(f'\n\n>>>Some atts were found by denovo search, but we are not going to use them!'
                       f'For more details, check the {genome.id}.atts file '
                       f'in the atts_denovo directory!\n')

    return genome


@task()
@genome_specific_logging
def scaffold_pipolins(genome: Genome) -> Pipolin:
    # Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
    # https://github.com/biopython/biopython/issues/1755
    logger = context.get('logger')

    # TODO: make this inside Scaffolder!
    if genome.features.is_on_the_same_contig(FeatureType.PIPOLB, FeatureType.ATT, FeatureType.TARGET_TRNA):
        logger.warning('>>> Scaffolding is not required!')
        return create_pipolin_fragments_single_contig(genome)
    else:
        logger.warning('>>> Scaffolding is required!')
        scaffolder = Scaffolder(genome=genome)
        return scaffolder.scaffold()


@task()
@genome_specific_logging
def extract_pipolin_regions(genome: Genome, pipolin: Pipolin, out_dir: str):
    genome_dict = create_seqio_records_dict(file=genome.file, file_format='fasta')

    pipolins_dir = os.path.join(out_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    logger = context.get('logger')

    with open(os.path.join(pipolins_dir, genome.id + '.fa'), 'w') as ouf:
        fragment_records = [create_fragment_record(fragment=f, genome_dict=genome_dict) for f in pipolin.fragments]

        for fragment_record, fragment in zip(fragment_records, pipolin.fragments):
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
def annotate_pipolins(genome: Genome, pipolins_dir, proteins, out_dir):
    prokka_results_dir = os.path.join(out_dir, 'prokka_results')
    os.makedirs(prokka_results_dir, exist_ok=True)
    run_prokka(genome_id=genome.id, pipolins_dir=pipolins_dir,
               proteins=proteins, prokka_results_dir=prokka_results_dir)
    return prokka_results_dir


@task()
@genome_specific_logging
def include_atts(genome: Genome, prokka_dir, out_dir, pipolin: Pipolin):
    gb_records = create_seqio_records_dict(file=os.path.join(prokka_dir, genome.id + '.gbk'),
                                           file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, genome.id + '.gff'))

    include_atts_into_gb(gb_records=gb_records, genome=genome, pipolin=pipolin)
    include_atts_into_gff(gff_records=gff_records, genome=genome, pipolin=pipolin)

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
