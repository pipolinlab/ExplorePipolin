"""
All the tasks should be defined in this file!
"""
import os

import pkg_resources
from prefect import task
from prefect import context
from prefect.engine import signals

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import Any, Optional, Sequence

from explore_pipolin.utilities.easyfig_coloring import add_colours
from explore_pipolin.common import Feature, FeatureType, RepeatPair, Pipolin, Genome, \
    define_genome_id, FeaturesContainer
from explore_pipolin.utilities.logging import genome_specific_logging
from explore_pipolin.utilities.misc import join_it, get_contig_orientation, \
    is_single_target_trna_per_contig, add_features_from_blast_entries
from explore_pipolin.utilities.io import read_blastxml, write_repeats, write_atts_denovo
from explore_pipolin.utilities.io import read_seqio_records
from explore_pipolin.utilities.io import read_aragorn_batch
from explore_pipolin.utilities.atts_denovo_search import find_repeats, is_att_denovo
from explore_pipolin.utilities.external_tools import tblastn_against_ref_pipolb, blastn_against_ref_att
from explore_pipolin.utilities.external_tools import run_prokka, run_aragorn
from explore_pipolin.utilities.misc import create_fragment_record
from explore_pipolin.utilities.io import read_gff_records
from explore_pipolin.utilities.including_atts_into_annotation import include_atts_into_gb
from explore_pipolin.utilities.including_atts_into_annotation import include_atts_into_gff
from explore_pipolin.utilities.io import write_genbank_records
from explore_pipolin.utilities.io import write_gff_records
from explore_pipolin.utilities.io import read_genome_contigs_from_file
from explore_pipolin.utilities.scaffolding import Scaffolder, create_pipolin_fragments_single_contig

_REF_PIPOLB = pkg_resources.resource_filename('explore_pipolin', 'data/pi-polB.faa')
_REF_ATT = pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa')


@task
def create_genome(genome_file) -> Genome:
    contigs = read_genome_contigs_from_file(genome_file=genome_file)
    genome_id = define_genome_id(genome_file)
    return Genome(genome_id=genome_id, genome_file=genome_file, contigs=contigs)


@task
@genome_specific_logging
def find_pipolbs(genome: Genome, out_dir) -> Genome:
    blast_results_dir = os.path.join(out_dir, 'polb_blast')
    os.makedirs(blast_results_dir, exist_ok=True)

    output_file = os.path.join(blast_results_dir, genome.genome_id) + '.fmt5'
    tblastn_against_ref_pipolb(genome_file=genome.genome_file, ref_pipolb=_REF_PIPOLB, output_file=output_file)
    entries = read_blastxml(blast_xml=output_file)

    add_features_from_blast_entries(entries=entries, feature_type=FeatureType.PIPOLB, genome=genome)

    return genome


@task
@genome_specific_logging
def find_atts(genome: Genome, out_dir) -> Genome:
    blast_results_dir = os.path.join(out_dir, 'att_blast')
    os.makedirs(blast_results_dir, exist_ok=True)

    output_file = os.path.join(blast_results_dir, genome.genome_id) + '.fmt5'
    blastn_against_ref_att(genome_file=genome.genome_file, ref_att=_REF_ATT, output_file=output_file)
    entries = read_blastxml(blast_xml=output_file)

    add_features_from_blast_entries(entries=entries, feature_type=FeatureType.ATT, genome=genome)

    return genome


@task
@genome_specific_logging
def find_trnas(genome: Genome, out_dir) -> Genome:
    aragorn_results_dir = os.path.join(out_dir, 'aragorn_results')
    os.makedirs(aragorn_results_dir, exist_ok=True)

    output_file = os.path.join(aragorn_results_dir, genome.genome_id + '.batch')
    run_aragorn(genome_file=genome.genome_file, output_file=output_file)
    entries = read_aragorn_batch(aragorn_batch=output_file)

    add_trna_features_from_aragorn_entries(entries=entries, genome=genome)
    find_and_add_target_trnas_features(genome.features)

    return genome


def add_trna_features_from_aragorn_entries(entries, genome: Genome):
    for contig_id, hits in entries.items():
        for hit in hits:
            trna_feature = Feature(start=hit[0], end=hit[1], strand=hit[2], contig_id=contig_id, genome=genome)
            genome.features.add_feature(feature=trna_feature, feature_type=FeatureType.TRNA)


def find_and_add_target_trnas_features(features: FeaturesContainer):
    for att in features.get_features(FeatureType.ATT):
        target_trna = features.find_overlapping_feature(att, FeatureType.TRNA)
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
def find_atts_denovo(genome: Genome, out_dir):
    logger = context.get('logger')

    if not genome.is_single_contig():
        logger.warning('This step is only for complete genomes. Skip...')
        raise signals.SKIP()

    atts_denovo_dir = os.path.join(out_dir, 'atts_denovo')
    os.makedirs(atts_denovo_dir, exist_ok=True)

    repeats: Sequence[RepeatPair] = find_repeats(genome, atts_denovo_dir)
    write_repeats(genome=genome, repeats=repeats, out_dir=atts_denovo_dir)

    atts_denovo: Sequence[RepeatPair] = [rep for rep in repeats if is_att_denovo(genome, rep)]
    for atts_pair in atts_denovo:
        genome.features.add_feature(feature=atts_pair.left, feature_type=FeatureType.ATT_DENOVO)
        genome.features.add_feature(feature=atts_pair.right, feature_type=FeatureType.ATT_DENOVO)

    for att in genome.features.get_features(FeatureType.ATT_DENOVO):
        target_trna = genome.features.find_overlapping_feature(att, FeatureType.TRNA)
        if target_trna is not None:
            genome.features.add_feature(feature=target_trna, feature_type=FeatureType.TARGET_TRNA_DENOVO)

    write_atts_denovo(genome.features.get_features(FeatureType.ATT_DENOVO), genome, atts_denovo_dir)

    return atts_denovo_dir


@task(skip_on_upstream_skip=False)
@genome_specific_logging
def are_atts_present(genome: Genome) -> Genome:
    logger = context.get('logger')

    num_atts = len(genome.features.get_features(FeatureType.ATT))
    num_atts_denovo = len(genome.features.get_features(FeatureType.ATT_DENOVO))
    if num_atts == 0 and num_atts_denovo == 0:
        logger.warning('\n\n>>>There is piPolB, but no atts were found! Not able to define pipolin bounds!\n')
        # TODO: probably, it makes sense to output piPolB(s) alone
        # raise signals.SKIP() # let's try cutting from both sides and proceed with annotation

    elif len(genome.features.get_features(FeatureType.ATT)) == 0:
        logger.warning(f'\n\n>>>No "usual" atts were found, but some atts were found by denovo search!'
                       f'For more details, check the {genome.genome_id}.atts file '
                       f'in the atts_denovo directory!\n')
        # TODO: check that it's only one repeat! Although, this shouldn't be a problem.
        atts_frames = [att.strand for att in genome.features.get_features(FeatureType.ATT_DENOVO)]
        if len(set(atts_frames)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_frames).pop() == genome.features.get_features(FeatureType.TARGET_TRNA_DENOVO)[0].strand:
            reverse_denovo_atts = []
            for att in genome.features.get_features(FeatureType.ATT_DENOVO):
                reverse_denovo_atts.append(Feature(start=att.start, end=att.end, strand=-att.strand,
                                                   contig_id=att.contig, genome=genome))
            [genome.features.add_feature(feature=i, feature_type=FeatureType.ATT) for i in reverse_denovo_atts]
        else:
            atts_denovo = genome.features.get_features(FeatureType.ATT_DENOVO)
            [genome.features.add_feature(feature=i, feature_type=FeatureType.ATT) for i in atts_denovo]
        target_trnas_denovo = genome.features.get_features(FeatureType.TARGET_TRNA_DENOVO)
        [genome.features.add_feature(feature=i, feature_type=FeatureType.TARGET_TRNA) for i in target_trnas_denovo]

    elif len(genome.features.get_features(FeatureType.ATT_DENOVO)) != 0:
        logger.warning(f'\n\n>>>Some atts were found by denovo search, but we are not going to use them!'
                       f'For more details, check the {genome.genome_id}.atts file '
                       f'in the atts_denovo directory!\n')

    return genome


@task()
@genome_specific_logging
def analyse_pipolin_orientation(genome: Genome) -> Genome:
    is_single_target_trna_per_contig(genome=genome)
    for contig in genome.contigs:
        contig.contig_orientation = get_contig_orientation(contig=contig, genome=genome)

    return genome


@task()
@genome_specific_logging
def scaffold_pipolins(genome: Genome) -> Pipolin:
    # Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
    # https://github.com/biopython/biopython/issues/1755
    logger = context.get('logger')

    if genome.features.is_on_the_same_contig(FeatureType.PIPOLB, FeatureType.ATT, FeatureType.TARGET_TRNA):
        logger.warning('>>> Scaffolding is not required!')
        return create_pipolin_fragments_single_contig(genome)
    else:
        logger.warning('>>> Scaffolding is required!')
        scaffolder = Scaffolder(genome=genome)
        return scaffolder.try_creating_single_record()


@task()
@genome_specific_logging
def extract_pipolin_regions(genome: Genome, pipolin: Pipolin, out_dir: str):
    genome_dict = read_seqio_records(file=genome.genome_file, file_format='fasta')

    pipolins_dir = os.path.join(out_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    logger = context.get('logger')

    with open(os.path.join(pipolins_dir, genome.genome_id + '.fa'), 'w') as ouf:
        fragment_records = [create_fragment_record(fragment=f, genome_dict=genome_dict) for f in pipolin.fragments]

        for fragment_record, fragment in zip(fragment_records, pipolin.fragments):
            logger.info(f'@pipolin fragment length {len(fragment_record)} from {fragment.contig.contig_id}')

        record = sum(join_it(fragment_records, SeqRecord(seq='N' * 100)), SeqRecord(seq=''))

        logger.info(f'@@@pipolin record total length {len(record)}')

        record.id = genome.genome_id
        record.name = genome.genome_id
        record.description = genome.genome_id
        SeqIO.write(sequences=record, handle=ouf, format='fasta')

    return pipolins_dir


@task()
@genome_specific_logging
def annotate_pipolins(genome: Genome, pipolins_dir, proteins, out_dir):
    prokka_results_dir = os.path.join(out_dir, 'prokka')
    os.makedirs(prokka_results_dir, exist_ok=True)
    run_prokka(genome_id=genome.genome_id, pipolins_dir=pipolins_dir,
               proteins=proteins, prokka_results_dir=prokka_results_dir)
    return prokka_results_dir


@task()
@genome_specific_logging
def include_atts_into_annotation(genome: Genome, prokka_dir, out_dir, pipolin: Pipolin):
    gb_records = read_seqio_records(file=os.path.join(prokka_dir, genome.genome_id + '.gbk'),
                                    file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, genome.genome_id + '.gff'))

    include_atts_into_gb(gb_records=gb_records, genome=genome, pipolin=pipolin)
    include_atts_into_gff(gff_records=gff_records, genome=genome, pipolin=pipolin)

    prokka_atts_dir = os.path.join(out_dir, 'prokka_atts')
    os.makedirs(prokka_atts_dir, exist_ok=True)

    write_genbank_records(gb_records=gb_records, out_dir=prokka_atts_dir, genome=genome)
    write_gff_records(gff_records=gff_records, out_dir=prokka_atts_dir, genome=genome)

    return prokka_atts_dir


@task()
@genome_specific_logging
def easyfig_add_colours(genome: Genome, in_dir):
    gb_records = read_seqio_records(file=os.path.join(in_dir, genome.genome_id + '.gbk'),
                                    file_format='genbank')
    add_colours(gb_records[genome.genome_id])
    write_genbank_records(gb_records=gb_records, out_dir=in_dir, genome=genome)
