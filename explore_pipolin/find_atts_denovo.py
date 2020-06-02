import os
from prefect import task
from prefect import context
from prefect.engine import signals
from explore_pipolin.utilities import GQuery
from explore_pipolin.utilities import save_left_right_subsequences
from explore_pipolin.utilities import blast_for_identical
from explore_pipolin.utilities import extract_repeats
from explore_pipolin.utilities import set_proper_location


@task()
def find_atts_denovo(genome, gquery: GQuery, root_dir):
    logger = context.get('logger')

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
