import click
from explore_pipolin.flow import get_flow

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir', type=click.Path(), required=True, help='Specify the output directory!')
@click.option('--add-colours', is_flag=True,
              help='Add colours to the final Genbank file features. '
                   'The genomic structure can be visualized further using Easyfig.')
def explore_pipolin(genomes, out_dir, add_colours):
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

    state = get_flow().run(genomes=genomes, out_dir=out_dir, add_colours=add_colours)
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolin()
