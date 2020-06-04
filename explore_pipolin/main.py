import click
from prefect.engine.state import State
from explore_pipolin.flow import get_flow

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# TODO: is it possible to write better?
def task_state_handler(_obj, old_state: State, _new_state):
    if 'gquery' in old_state.cached_inputs:
        gquery = old_state.cached_inputs['gquery'].value
        gquery_id = gquery.genome.id
        print(f'>>> It was {gquery_id}...')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir', type=click.Path(), required=True)
@click.option('--abricate_dir', type=click.Path(exists=True))
def explore_pipolin(genomes, out_dir, abricate_dir):
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

    state = get_flow().run(genomes=genomes, out_dir=out_dir, abricate_dir=abricate_dir,
                           task_runner_state_handlers=[task_state_handler])
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolin()
