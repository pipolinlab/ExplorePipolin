import datetime
import logging
import os
from typing import Optional

import click

from explore_pipolin.flow import get_flow

from explore_pipolin.utilities.external_tools import check_external_dependencies
from explore_pipolin.utilities.logging import set_logging_dir
from explore_pipolin.common import CONTEXT_SETTINGS

from explore_pipolin.tasks.reconstruct_pipolins import NO_BORDER_INFLATE


def check_genome_file_names(genome):
    if len(genome) > 1:
        if len(genome) != len(set(os.path.basename(i) for i in genome)):
            logging.fatal('GENOME files should have different names!')
            exit(1)


_DEFAULT_OUT_DIR_PREFIX = 'results'
_SUFFIX = datetime.datetime.now().strftime('_%H%M%S')


def get_out_dir_name(out_dir_prefix: Optional[str], out_dir: Optional[str]) -> str:
    if out_dir_prefix and out_dir:
        logging.fatal('Options --out-dir-prefix and --out-dir are mutually exclusive!')
        exit(1)
    elif out_dir_prefix:
        return os.path.join(os.getcwd(), out_dir_prefix + _SUFFIX)
    elif out_dir:
        return out_dir
    else:
        return os.path.join(os.getcwd(), _DEFAULT_OUT_DIR_PREFIX + _SUFFIX)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genome', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir-prefix', type=str,
              help=f'Use this prefix for the output directory, '
                   f'instead of the default "{_DEFAULT_OUT_DIR_PREFIX}" prefix.')
@click.option('--out-dir', type=click.Path(exists=True),
              help='Use the existing output directory. If the directory contains results of a previous run, '
                   'such as found piPolBs, ATTs and tRNAs, the program will reuse them, unless '
                   '--do-not-reuse option is specified.')
@click.option('--pipolb-hmm-profile', type=click.Path(exists=True),
              help='If not provided, the default profile will be used instead.')
@click.option('--ref-att', type=click.Path(exists=True),
              help='Att sequence in FASTA file to use as 1st priority. '
                   'If not provided, the default file will be used instead.')
@click.option('--percent-identity', type=int, default=85, show_default=True,
              help='Minimum percent identity in direct repeats search')
@click.option('--max-inflate', type=int, default=NO_BORDER_INFLATE, show_default=True,
              help='If no borders of pipolin are found (no ATTs), '
                   'inflate the analysed region from both sides of piPolB.')
@click.option('--no-annotation', is_flag=True, help='Do not run the annotation step (i.e. Prokka).')
@click.option('--proteins', type=click.Path(exists=True),
              help='Prokka param: FASTA or GBK file to use as 1st priority. '
                   'If not provided, the default file will be used instead.')
@click.option('--skip-colours', is_flag=True,
              help='Do not add an Easyfig-compatible colouring scheme to the final Genbank file features.')
@click.option('--cpus', default=8, type=int, show_default=True,
              help='Prokka param: Number of CPUs to use [0=all]')
@click.option('--do-not-reuse', is_flag=True,
              help='Do not reuse information about piPolBs, ATTs and tRNAs found in a previous run for '
                   'the same genome. I.e. it will run all the analysis from scratch for the given genome.')
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
        do_not_reuse,
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

    out_dir_name = get_out_dir_name(out_dir_prefix, out_dir)
    set_logging_dir(out_dir_name)

    state = get_flow().run(
        genome_file=genome,
        out_dir=out_dir_name,
        pipolb_hmm_profile=pipolb_hmm_profile,
        ref_att=ref_att,
        percent_identity=percent_identity,
        max_inflate=max_inflate,
        no_annotation=no_annotation,
        proteins=proteins,
        skip_colours=skip_colours,
        cpus=cpus,
        do_not_reuse=do_not_reuse,
    )
    assert state.is_successful()


if __name__ == '__main__':
    main()
