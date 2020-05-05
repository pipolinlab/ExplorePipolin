import os
import prefect
from prefect import task
from prefect.engine import signals

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
def create_gquery(genome) -> GQuery:
    gquery = GQuery(gquery_id=define_gquery_id(genome=genome))
    genome_dict = read_seqio_records(file=genome, file_format='fasta')
    for key, value in genome_dict.items():
        contig = Contig(contig_id=key, contig_length=len(value.seq))
        gquery.contigs.append(contig)

    return gquery


@task
def run_blast_against_polb(genome, root_dir, reference):
    blast_path = os.path.join(root_dir, 'polb_blast')
    blast_genome_against_seq(genome=genome, seq=reference, seq_type='protein', output_dir=blast_path)
    return blast_path


@task
def run_blast_against_att(genome, root_dir, reference):
    blast_path = os.path.join(root_dir, 'att_blast')
    blast_genome_against_seq(genome=genome, seq=reference, seq_type='nucleotide', output_dir=blast_path)
    return blast_path


@task
def add_features_from_blast(gquery: GQuery, blast_dir, feature_type):
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
def add_features_from_aragorn(gquery: GQuery, aragorn_dir):
    entries = read_aragorn_batch(aragorn_batch=os.path.join(aragorn_dir, f'{gquery.gquery_id}.batch'))
    for contig_id, hits in entries.items():
        for hit in hits:
            feature = Feature(start=hit[0], end=hit[1], frame=hit[2], contig=gquery.get_contig_by_id(contig_id))
            gquery.trnas.append(feature)

    for att in gquery.atts:
        target_trna = gquery.find_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas.append(target_trna)


@task()
def find_atts_denovo(genome, gquery: GQuery, root_dir):
    logger = prefect.context.get('logger')

    if not gquery.is_single_contig():
        logger.warning('This step is only for complete genomes. Pass...')
        raise signals.SKIP()

    repeats_dir = os.path.join(root_dir, 'atts_denovo')

    left_window, right_window = gquery.get_left_right_windows()
    save_left_right_subsequences(genome=genome, left_window=left_window, right_window=right_window,
                                 repeats_dir=repeats_dir)

    blast_for_identical(gquery_id=gquery.gquery_id, repeats_dir=repeats_dir)
    left_repeats, right_repeats, sequences = extract_repeats(file=os.path.join(repeats_dir, gquery.gquery_id + '.fmt5'))
    left_repeats = [set_proper_location(seq_range=i, shift=left_window[0]) for i in left_repeats]
    right_repeats = [set_proper_location(seq_range=i, shift=right_window[0]) for i in right_repeats]

    with open(os.path.join(repeats_dir, gquery.gquery_id + '.repeats'), 'w') as ouf:
        polbs_locations = sorted([(x.start, x.end) for x in gquery.polbs], key=lambda x: x[0])
        print('left_repeat', 'right_repeat', 'length', 'polbs',
              'd_to_the_left', 'd_to_the_right', 'sequence', sep='\t', file=ouf)
        for repeat in zip(left_repeats, right_repeats, sequences):
            print(repeat[0], repeat[1], repeat[0][1] - repeat[0][0], ','.join([str(i) for i in polbs_locations]),
                  polbs_locations[0][0] - repeat[0][1], repeat[1][0] - polbs_locations[-1][-1], repeat[2],
                  sep='\t', file=ouf)

    atts_denovo = [(i, j) for i, j in zip(left_repeats, right_repeats) if gquery.is_att_denovo(i, j)]

    with open(os.path.join(repeats_dir, gquery.gquery_id + '.atts'), 'w') as ouf:
        print('attL_start', 'attL_end', 'attR_start', 'attR_end', sep='\t', file=ouf)
        for att_pair in atts_denovo:
            print(att_pair[0][0], att_pair[0][1], att_pair[1][0], att_pair[1][1], sep='\t', file=ouf)

    return repeats_dir


@task
def add_features_atts_denovo(gquery: GQuery, atts_denovo_dir):
    with open(os.path.join(atts_denovo_dir, gquery.gquery_id + '.atts')) as inf:
        _ = inf.readline()
        for line in inf:
            att_pair = [int(i) for i in line.strip().split(sep='\t')]
            gquery.denovo_atts.append(Feature(start=att_pair[0], end=att_pair[1],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))
            gquery.denovo_atts.append(Feature(start=att_pair[2], end=att_pair[3],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))

    for att in gquery.denovo_atts:
        print(type(att.start), type(att.end))
        target_trna = gquery.find_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas_denovo.append(target_trna)


@task
def are_polbs_present(gquery: GQuery):
    logger = prefect.context.get('logger')

    if len(gquery.polbs) == 0:
        logger.warning('No piPolB! => No pipolins!')
        raise signals.FAIL()


@task(skip_on_upstream_skip=False)
def are_atts_present(gquery):
    logger = prefect.context.get('logger')

    if len(gquery.atts) == 0 and len(gquery.denovo_atts) == 0:
        logger.warning('\n\n>>>There is piPolB, but no atts were found! Not able to define pipolin bounds!\n')
        # TODO: probably, it makes sense to output piPolB(s) alone
        # raise signals.SKIP() # let's try cutting from both sides and proceed with annotation

    elif len(gquery.atts) == 0:
        logger.warning(f'\n\n>>>No "usual" atts were found, but some atts were found by denovo search!'
                       f'For more details, check the {gquery.gquery_id}.atts file in the atts_denovo directory!\n')
        # TODO: check that it's only one repeat! Although, this shouldn't be a problem.
        gquery.atts.extend(gquery.denovo_atts)
        gquery.target_trnas.extend(gquery.target_trnas_denovo)

    elif len(gquery.denovo_atts) != 0:
        logger.warning(f'\n\n>>>Some atts were found by denovo search, but we are not going to use them!'
                       f'For more details, check the {gquery.gquery_id}.atts file in the atts_denovo directory!\n')
