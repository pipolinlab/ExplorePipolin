import logging
import os

import click

from explore_pipolin.flow import get_flow

from explore_pipolin.utilities.external_tools import check_aragorn, check_blast, check_prokka
from explore_pipolin.utilities.logging import set_logging_dir
from explore_pipolin.common import CONTEXT_SETTINGS


def check_external_dependencies():
    check_blast()
    check_aragorn()
    check_prokka()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genome', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir', type=click.Path(), required=True, help='Specify the output directory!')
@click.option('--add-colours', is_flag=True,
              help='Add colours to the final Genbank file features. '
                   'The genomic structure can be visualized further using Easyfig.')
def explore_pipolin(genome, out_dir, add_colours):
    """
    ExplorePipolin is a search tool that identifies and analyses pipolin elements within bacterial genome(s).
    """
    # from graphviz import Digraph
    # graph: Digraph = get_flow().visualize(filename='/dev/null')
    # graph.node_attr.update(fontsize='18', fontname='DejaVuSansMono')
    # graph.render('data/test')
    # # modify test in text editor
    # # dot -Tpdf test > test.pdf
    # # check the result

    if len(genome) > 1:
        if len(genome) != len(set(os.path.basename(i) for i in genome)):
            logging.fatal('GENOME files should have different names!')
            exit(1)

    check_external_dependencies()

    logging_dir = os.path.join(out_dir, 'results')
    os.makedirs(logging_dir, exist_ok=True)
    set_logging_dir(logging_dir)

    state = get_flow().run(genome_file=genome, out_dir=out_dir, add_colours=add_colours)
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolin()
