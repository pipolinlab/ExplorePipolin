import click
import pkg_resources
from prefect import Flow, Parameter, unmapped
from prefect.tasks.core.constants import Constant
from pipolin_finder.utilities import CONTEXT_SETTINGS
from pipolin_finder.identify_pipolins_roughly import create_gquery
from pipolin_finder.identify_pipolins_roughly import run_blast_against_ref
from pipolin_finder.identify_pipolins_roughly import add_features_from_blast
from pipolin_finder.identify_pipolins_roughly import detect_trnas_with_aragorn
from pipolin_finder.identify_pipolins_roughly import add_features_from_aragorn
from pipolin_finder.identify_pipolins_roughly import are_polbs_present
from pipolin_finder.identify_pipolins_roughly import find_atts_denovo
from pipolin_finder.identify_pipolins_roughly import add_features_atts_denovo
from pipolin_finder.identify_pipolins_roughly import are_atts_present
from pipolin_finder.analyse_pipolin_orientation import analyse_pipolin_orientation
from pipolin_finder.scaffold_gapped_pipolins import scaffold_pipolins
from pipolin_finder.extract_pipolin_regions import extract_pipolin_regions
from pipolin_finder.annotate_pipolins import annotate_pipolins
from pipolin_finder.include_atts_into_annotation import include_atts_into_annotation
from pipolin_finder.easyfig_add_colours import easyfig_add_colours

REF_POLB = Constant(pkg_resources.resource_filename('pipolin_finder', 'data/pi-polB.faa'))
REF_ATT = Constant(pkg_resources.resource_filename('pipolin_finder', 'data/attL.fa'))
PROTEINS = Constant(pkg_resources.resource_filename('pipolin_finder', '/data/HHpred_proteins.faa'))


def get_flow():
    with Flow('MAIN') as flow:
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')

        gquery = create_gquery.map(genome=genomes)

        polbs_blast = run_blast_against_ref.map(genome=genomes, root_dir=unmapped(out_dir),
                                                reference=unmapped(REF_POLB),
                                                dir_name=unmapped(Constant('polb_blast')))
        t_add_polbs = add_features_from_blast.map(gquery=gquery, blast_dir=polbs_blast,
                                                  feature_type=unmapped(Constant('polbs')))

        atts_blast = run_blast_against_ref.map(genome=genomes, root_dir=unmapped(out_dir),
                                               reference=unmapped(REF_ATT),
                                               dir_name=unmapped(Constant('att_blast')))
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
        easyfig_add_colours.map(gquery=gquery, in_dir=prokka_atts, abricate_dir=unmapped(Constant(None)))

    return flow


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--out-dir', type=click.Path(), required=True)
def explore_pipolins(genomes, out_dir):
    """
    PipolinFinder is a search tool that identifies, extracts and annotates pipolin elements within bacterial genome(s).
    """

    # get_flow().visualize()

    state = get_flow().run(genomes=genomes, out_dir=out_dir)
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolins()
