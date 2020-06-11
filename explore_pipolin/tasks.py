"""
All the tasks should be defined in this file!
"""
import os
from prefect import task
from prefect import context
from prefect.engine import signals

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import Any, Optional, Sequence

from explore_pipolin.utilities.easyfig_coloring import add_colours
from explore_pipolin.common import Feature, FeatureType, RepeatPair, Pipolin, GenomeFeatures
from explore_pipolin.utilities.logging import genome_specific_logging
from explore_pipolin.utilities.misc import feature_from_blasthit, join_it, get_contig_orientation, \
    is_single_target_trna_per_contig
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
from explore_pipolin.utilities.io import create_genome_from_file
from explore_pipolin.utilities.scaffolding import Scaffolder, create_pipolin_fragments_single_contig


@task
def create_features_container(genome_file) -> GenomeFeatures:
    return GenomeFeatures(genome=create_genome_from_file(genome_file))


@task()
@genome_specific_logging
def run_blast_against_pipolb(features_container, ref_pipolb, out_dir):
    blast_results_dir = os.path.join(out_dir, 'polb_blast')
    os.makedirs(blast_results_dir, exist_ok=True)
    tblastn_against_ref_pipolb(genome=features_container.genome.genome_file, ref_pipolb=ref_pipolb,
                               out_dir=blast_results_dir)
    return blast_results_dir


@task()
@genome_specific_logging
def run_blast_against_att(features_container, ref_att, out_dir):
    blast_results_dir = os.path.join(out_dir, 'att_blast')
    os.makedirs(blast_results_dir, exist_ok=True)
    blastn_against_ref_att(genome=features_container.genome.genome_file, ref_att=ref_att, out_dir=blast_results_dir)
    return blast_results_dir


@task()
@genome_specific_logging
def add_features_from_blast(features_container: GenomeFeatures, blast_dir, feature_type: FeatureType) -> GenomeFeatures:
    entries = read_blastxml(blast_xml=os.path.join(blast_dir, f'{features_container.genome.genome_id}.fmt5'))
    for entry in entries:
        for hit in entry:
            feature = feature_from_blasthit(hit=hit, contig_id=entry.id, genome=features_container.genome)
            features_container.get_features(feature_type).append(feature)

    return features_container


@task()
@genome_specific_logging
def detect_trnas_with_aragorn(features_container, out_dir):
    aragorn_results_dir = os.path.join(out_dir, 'aragorn_results')
    os.makedirs(aragorn_results_dir, exist_ok=True)
    run_aragorn(genome=features_container.genome.genome_file, aragorn_results_dir=aragorn_results_dir)
    return aragorn_results_dir


@task()
@genome_specific_logging
def add_features_from_aragorn(features_container: GenomeFeatures, aragorn_results_dir) -> GenomeFeatures:
    entries = read_aragorn_batch(aragorn_batch=os.path.join(aragorn_results_dir,
                                                            f'{features_container.genome.genome_id}.batch'))
    for contig_id, hits in entries.items():
        for hit in hits:
            feature = Feature(start=hit[0], end=hit[1], strand=hit[2], contig_id=contig_id,
                              genome=features_container.genome)
            features_container.get_features(FeatureType.TRNA).append(feature)
    for att in features_container.get_features(FeatureType.ATT):
        target_trna = features_container.find_overlapping_feature(att, FeatureType.TRNA)
        if target_trna is not None:
            features_container.get_features(FeatureType.TARGET_TRNA).append(target_trna)

    return features_container


@task()
@genome_specific_logging
def are_pipolbs_present(features_container: GenomeFeatures):
    logger = context.get('logger')

    if len(features_container.get_features(FeatureType.PIPOLB)) == 0:
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
def find_atts_denovo(features_container: GenomeFeatures, out_dir):
    logger = context.get('logger')

    if not features_container.genome.is_single_contig():
        logger.warning('This step is only for complete genomes. Skip...')
        raise signals.SKIP()

    atts_denovo_dir = os.path.join(out_dir, 'atts_denovo')
    os.makedirs(atts_denovo_dir, exist_ok=True)

    repeats: Sequence[RepeatPair] = find_repeats(features_container, atts_denovo_dir)
    write_repeats(features_container=features_container, repeats=repeats, out_dir=atts_denovo_dir)

    atts_denovo: Sequence[RepeatPair] = [rep for rep in repeats if is_att_denovo(features_container, rep)]
    for atts_pair in atts_denovo:
        features_container.get_features(FeatureType.ATT_DENOVO).append(atts_pair.left)
        features_container.get_features(FeatureType.ATT_DENOVO).append(atts_pair.right)

    for att in features_container.get_features(FeatureType.ATT_DENOVO):
        target_trna = features_container.find_overlapping_feature(att, FeatureType.TRNA)
        if target_trna is not None:
            features_container.get_features(FeatureType.TARGET_TRNA_DENOVO).append(target_trna)

    write_atts_denovo(features_container.get_features(FeatureType.ATT_DENOVO),
                      features_container.genome, atts_denovo_dir)

    return atts_denovo_dir


@task(skip_on_upstream_skip=False)
@genome_specific_logging
def are_atts_present(features_container: GenomeFeatures):
    logger = context.get('logger')

    if len(features_container.get_features(
            FeatureType.ATT)) == 0 and len(features_container.get_features(FeatureType.ATT_DENOVO)) == 0:
        logger.warning('\n\n>>>There is piPolB, but no atts were found! Not able to define pipolin bounds!\n')
        # TODO: probably, it makes sense to output piPolB(s) alone
        # raise signals.SKIP() # let's try cutting from both sides and proceed with annotation

    elif len(features_container.get_features(FeatureType.ATT)) == 0:
        logger.warning(f'\n\n>>>No "usual" atts were found, but some atts were found by denovo search!'
                       f'For more details, check the {features_container.genome.genome_id}.atts file '
                       f'in the atts_denovo directory!\n')
        # TODO: check that it's only one repeat! Although, this shouldn't be a problem.
        atts_frames = [att.strand for att in features_container.get_features(FeatureType.ATT_DENOVO)]
        if len(set(atts_frames)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_frames).pop() == features_container.get_features(FeatureType.TARGET_TRNA_DENOVO)[0].strand:
            reverse_denovo_atts = []
            for att in features_container.get_features(FeatureType.ATT_DENOVO):
                reverse_denovo_atts.append(Feature(start=att.start, end=att.end, strand=-att.strand,
                                                   contig_id=att.contig, genome=features_container.genome))
            features_container.get_features(FeatureType.ATT).extend(reverse_denovo_atts)
        else:
            features_container.get_features(FeatureType.ATT).extend(
                features_container.get_features(FeatureType.ATT_DENOVO))
        features_container.get_features(FeatureType.TARGET_TRNA).extend(
            features_container.get_features(FeatureType.TARGET_TRNA_DENOVO))

    elif len(features_container.get_features(FeatureType.ATT_DENOVO)) != 0:
        logger.warning(f'\n\n>>>Some atts were found by denovo search, but we are not going to use them!'
                       f'For more details, check the {features_container.genome.genome_id}.atts file '
                       f'in the atts_denovo directory!\n')

    return features_container


@task()
@genome_specific_logging
def analyse_pipolin_orientation(features_container: GenomeFeatures):
    is_single_target_trna_per_contig(features_container=features_container)
    for contig in features_container.genome.contigs:
        contig.contig_orientation = get_contig_orientation(contig=contig, features_container=features_container)

    return features_container


@task()
@genome_specific_logging
def scaffold_pipolins(features_container: GenomeFeatures) -> Pipolin:
    # Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
    # https://github.com/biopython/biopython/issues/1755
    logger = context.get('logger')

    if features_container.is_on_the_same_contig(FeatureType.PIPOLB, FeatureType.ATT, FeatureType.TARGET_TRNA):
        logger.warning('>>> Scaffolding is not required!')
        return create_pipolin_fragments_single_contig(features_container)
    else:
        logger.warning('>>> Scaffolding is required!')
        scaffolder = Scaffolder(features_container=features_container)
        return scaffolder.try_creating_single_record()


@task()
@genome_specific_logging
def extract_pipolin_regions(features_container: GenomeFeatures, pipolin: Pipolin, out_dir: str):
    genome_dict = read_seqio_records(file=features_container.genome.genome_file, file_format='fasta')

    pipolins_dir = os.path.join(out_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    logger = context.get('logger')

    with open(os.path.join(pipolins_dir, features_container.genome.genome_id + '.fa'), 'w') as ouf:
        fragment_records = [create_fragment_record(fragment=f, genome_dict=genome_dict) for f in pipolin.fragments]

        for fragment_record, fragment in zip(fragment_records, pipolin.fragments):
            logger.info(f'@pipolin fragment length {len(fragment_record)} from {fragment.contig.contig_id}')

        record = sum(join_it(fragment_records, SeqRecord(seq='N' * 100)), SeqRecord(seq=''))

        logger.info(f'@@@pipolin record total length {len(record)}')

        record.id = features_container.genome.genome_id
        record.name = features_container.genome.genome_id
        record.description = features_container.genome.genome_id
        SeqIO.write(sequences=record, handle=ouf, format='fasta')

    return pipolins_dir


@task()
@genome_specific_logging
def annotate_pipolins(features_container, pipolins_dir, proteins, out_dir):
    prokka_results_dir = os.path.join(out_dir, 'prokka')
    os.makedirs(prokka_results_dir, exist_ok=True)
    run_prokka(genome_id=features_container.genome.genome_id, pipolins_dir=pipolins_dir,
               proteins=proteins, prokka_results_dir=prokka_results_dir)
    return prokka_results_dir


@task()
@genome_specific_logging
def include_atts_into_annotation(features_container, prokka_dir, out_dir, pipolin: Pipolin):
    gb_records = read_seqio_records(file=os.path.join(prokka_dir, features_container.genome.genome_id + '.gbk'),
                                    file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, features_container.genome.genome_id + '.gff'))

    include_atts_into_gb(gb_records=gb_records, genome=features_container.genome, pipolin=pipolin)
    include_atts_into_gff(gff_records=gff_records, genome=features_container.genome, pipolin=pipolin)

    prokka_atts_dir = os.path.join(out_dir, 'prokka_atts')
    os.makedirs(prokka_atts_dir, exist_ok=True)

    write_genbank_records(gb_records=gb_records, out_dir=prokka_atts_dir, genome=features_container.genome)
    write_gff_records(gff_records=gff_records, out_dir=prokka_atts_dir, genome=features_container.genome)

    return prokka_atts_dir


@task()
@genome_specific_logging
def easyfig_add_colours(features_container: GenomeFeatures, in_dir):
    gb_records = read_seqio_records(file=os.path.join(in_dir, features_container.genome.genome_id + '.gbk'),
                                    file_format='genbank')
    add_colours(gb_records[features_container.genome.genome_id])
    write_genbank_records(gb_records=gb_records, out_dir=in_dir, genome=features_container.genome)
