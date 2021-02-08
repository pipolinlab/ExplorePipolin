import logging
import os

import click

from explore_pipolin.flow import get_flow

from explore_pipolin.utilities.external_tools import check_external_dependencies
from explore_pipolin.utilities.logging import set_logging_dir
from explore_pipolin.common import CONTEXT_SETTINGS


def check_genome_file_names(genome):
    if len(genome) > 1:
        if len(genome) != len(set(os.path.basename(i) for i in genome)):
            logging.fatal('GENOME files should have different names!')
            exit(1)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genome', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir', type=click.Path(), required=True, help='Output directory')
@click.option('--add-colours', is_flag=True,
              help='Add colours to the final Genbank file features. '
                   'The genomic structure can be visualized further using Easyfig.')
@click.option('--pipolb-hmm-profile', type=click.Path(exists=True),
              help='If not provided, the default profile will be used instead.')
@click.option('--ref-att', type=click.Path(exists=True),
              help='Att sequence in FASTA file to use as 1st priority. '
                   'If not provided, the default file will be used instead.')
@click.option('--perc-identity', type=int, default=85, show_default=True,
              help='Minimum percent identity in direct repeats search')
@click.option('--proteins', type=click.Path(exists=True),
              help='Prokka param: FASTA or GBK file to use as 1st priority. '
                   'If not provided, the default file will be used instead.')
@click.option('--cpus', default=8, type=int, show_default=True,
              help='Prokka param: Number of CPUs to use [0=all]')
def explore_pipolin(
        genome,
        out_dir,
        add_colours,
        pipolb_hmm_profile,
        ref_att,
        perc_identity,
        proteins,
        cpus,
):
    """
    ExplorePipolin is a search tool for prediction and analysis of pipolins, bacterial mobile genetic elements.
    """

    check_genome_file_names(genome)

    check_external_dependencies()

    logging_dir = os.path.join(out_dir, 'logs')
    os.makedirs(logging_dir, exist_ok=True)
    set_logging_dir(logging_dir)

    state = get_flow().run(
        genome_file=genome,
        out_dir=out_dir,
        add_colours=add_colours,
        pipolb_hmm_profile=pipolb_hmm_profile,
        ref_att=ref_att,
        perc_identity=perc_identity,
        proteins=proteins,
        cpus=cpus,
    )
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolin()
