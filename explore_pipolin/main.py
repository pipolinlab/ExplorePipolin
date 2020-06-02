import click
import pkg_resources
from prefect import Flow, Parameter, unmapped
from prefect.engine.state import State
from prefect.tasks.core.constants import Constant
from explore_pipolin.utilities import CONTEXT_SETTINGS
from explore_pipolin.identify_pipolins_roughly import create_gquery
from explore_pipolin.identify_pipolins_roughly import run_blast_against_polb
from explore_pipolin.identify_pipolins_roughly import run_blast_against_att
from explore_pipolin.identify_pipolins_roughly import add_features_from_blast
from explore_pipolin.identify_pipolins_roughly import detect_trnas_with_aragorn
from explore_pipolin.identify_pipolins_roughly import add_features_from_aragorn
from explore_pipolin.identify_pipolins_roughly import are_polbs_present
from explore_pipolin.identify_pipolins_roughly import find_atts_denovo
from explore_pipolin.identify_pipolins_roughly import add_features_atts_denovo
from explore_pipolin.identify_pipolins_roughly import are_atts_present
from explore_pipolin.analyse_pipolin_orientation import analyse_pipolin_orientation
from explore_pipolin.scaffold_gapped_pipolins import scaffold_pipolins
from explore_pipolin.extract_pipolin_regions import extract_pipolin_regions
from explore_pipolin.annotate_pipolins import annotate_pipolins
from explore_pipolin.include_atts_into_annotation import include_atts_into_annotation
from explore_pipolin.easyfig_add_colours import easyfig_add_colours
from explore_pipolin.include_atts_into_annotation import set_correct_positions

REF_POLB = Constant(pkg_resources.resource_filename('explore_pipolin', 'data/pi-polB.faa'))
REF_ATT = Constant(pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa'))
PROTEINS = Constant(pkg_resources.resource_filename('explore_pipolin', '/data/HHpred_proteins.faa'))


def get_flow():
    with Flow('MAIN') as flow:
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')
        abricate_dir = Parameter('abricate_dir')

        gquery = create_gquery.map(genome=genomes)

        polbs_blast = run_blast_against_polb.map(genome=genomes, root_dir=unmapped(out_dir),
                                                 reference=unmapped(REF_POLB))
        t_add_polbs = add_features_from_blast.map(gquery=gquery, blast_dir=polbs_blast,
                                                  feature_type=unmapped(Constant('polbs')))

        atts_blast = run_blast_against_att.map(genome=genomes, root_dir=unmapped(out_dir),
                                               reference=unmapped(REF_ATT))
        t_add_atts = add_features_from_blast.map(gquery=gquery, blast_dir=atts_blast,
                                                 feature_type=unmapped(Constant('atts')))

        aragorn_results = detect_trnas_with_aragorn.map(genome=genomes, root_dir=unmapped(out_dir))
        t_add_trnas = add_features_from_aragorn.map(gquery=gquery, aragorn_dir=aragorn_results,
                                                    upstream_tasks=[t_add_polbs, t_add_atts])

        t_check_polbs = are_polbs_present.map(gquery=gquery, upstream_tasks=[t_add_trnas])

        atts_denovo = find_atts_denovo.map(genome=genomes, gquery=gquery, root_dir=unmapped(out_dir),
                                           upstream_tasks=[t_check_polbs])
        t_add_denovo_atts = add_features_atts_denovo.map(gquery=gquery, atts_denovo_dir=atts_denovo)

        t_check_features = are_atts_present.map(gquery=gquery, upstream_tasks=[t_add_denovo_atts])

        t_set_orientations = analyse_pipolin_orientation.map(gquery=gquery, upstream_tasks=[t_check_features])
        t_scaffolding = scaffold_pipolins.map(gquery=gquery, upstream_tasks=[t_set_orientations])

        pipolin_sequences = extract_pipolin_regions.map(genome=genomes, gquery=gquery,
                                                        root_dir=unmapped(out_dir), upstream_tasks=[t_scaffolding])
        prokka = annotate_pipolins.map(gquery=gquery, pipolins_dir=pipolin_sequences,
                                       proteins=unmapped(PROTEINS), root_dir=unmapped(out_dir))
        prokka_atts = include_atts_into_annotation.map(gquery=gquery, prokka_dir=prokka,
                                                       root_dir=unmapped(out_dir))
        # TODO: make this task optional
        t_easyfig = easyfig_add_colours.map(gquery=gquery, in_dir=prokka_atts,
                                            abricate_dir=unmapped(abricate_dir),
                                            upstream_tasks=[t_scaffolding])
        # TODO: when easyfig skipped, set (skip_on_upstream_skip=False)
        # prokka_atts_positions = set_correct_positions.map(gquery=gquery, prokka_atts_dir=prokka_atts,
        #                                                   root_dir=unmapped(out_dir), upstream_tasks=[t_easyfig])

    return flow


# TODO: is it possible to write better?
def task_state_handler(obj, old_state: State, new_state):
    if 'gquery' in old_state.cached_inputs:
        gquery_id = old_state.cached_inputs['gquery'].value.gquery_id
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