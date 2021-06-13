import logging
import os

import click

from explore_pipolin.flow import get_flow
import explore_pipolin.settings as settings
from explore_pipolin.settings import GlobalSettings, _DEFAULT_OUT_DIR_PREFIX, _NO_BORDER_INFLATE

from explore_pipolin.utilities.external_tools import check_external_dependencies
from explore_pipolin.common import CONTEXT_SETTINGS


def check_genome_file_names(genome):
    if len(genome) > 1:
        if len(genome) != len(set(os.path.basename(i) for i in genome)):
            logging.fatal('GENOME files should have different names!')
            exit(1)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genome', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir-prefix', type=str,
              help=f'Use this prefix for the output directory, '
                   f'instead of the default "{_DEFAULT_OUT_DIR_PREFIX}" prefix.')
@click.option('--out-dir', type=click.Path(),
              help='Use this output directory instead.')
@click.option('--pipolb-hmm-profile', type=click.Path(exists=True),
              help="piPolB's HMM profile to use as 1st priority."
                   'If not provided, the default profile will be used instead.')
@click.option('--ref-att', type=click.Path(exists=True),
              help='Att sequence in FASTA file to use as 1st priority. '
                   'If not provided, the default file will be used instead.')
@click.option('--percent-identity', type=int, default=85, show_default=True,
              help='Minimum percent identity in direct repeats search')
@click.option('--max-inflate', type=int, default=_NO_BORDER_INFLATE, show_default=True,
              help='If no borders of pipolin are found (no ATTs), '
                   'inflate the analysed region from both sides of piPolB.')
@click.option('--no-annotation', is_flag=True, help='Do not run the annotation step (i.e. Prokka).')
@click.option('--proteins', type=click.Path(exists=True),
              help='Prokka param: FASTA or GBK file to use as 1st priority. '
                   'If not provided, the default file will be used instead.')
@click.option('--skip-colours', is_flag=True,
              help='Do not add an Easyfig-compatible colouring scheme to the final Genbank file.')
@click.option('--cpus', default=8, type=int, show_default=True,
              help='Prokka param: Number of CPUs to use [0=all]')
def main(
        genome,
        out_dir_prefix,
        out_dir,
        pipolb_hmm_profile,
        ref_att,
        percent_identity,
        max_inflate,
        no_annotation,
        proteins,
        skip_colours,
        cpus,
):
    """
    ExplorePipolin is a search tool for prediction and analysis of pipolins, bacterial mobile genetic elements.
    """

    # from graphviz import Digraph
    # graph: Digraph = get_flow().visualize(filename='/Users/liubov/Documents/ExplorePipolin_data/testtest')
    # graph.node_attr.update(fontsize='18', fontname='DejaVuSansMono')
    # graph.render('/Users/liubov/Documents/ExplorePipolin_data/test')
    # # modify test in text editor
    # # dot -Tpdf test > test.pdf
    # # check the result

    check_genome_file_names(genome)

    check_external_dependencies()

    settings.set_instance(GlobalSettings.create_instance(
        out_dir_prefix, out_dir, pipolb_hmm_profile, ref_att, percent_identity, max_inflate, proteins, cpus
    ))
    os.makedirs(settings.get_instance().out_dir, exist_ok=True)

    state = get_flow().run(
        genome_file=genome,
        no_annotation=no_annotation,
        skip_colours=skip_colours,
    )
    assert state.is_successful()


if __name__ == '__main__':
    main()
