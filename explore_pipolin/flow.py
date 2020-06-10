import pkg_resources
from prefect import Flow, Parameter, unmapped, case
from prefect.tasks.control_flow import FilterTask
from prefect.tasks.core.constants import Constant

from explore_pipolin import tasks
from explore_pipolin.common import FeatureType

_REF_PIPOLB = Constant(pkg_resources.resource_filename('explore_pipolin', 'data/pi-polB.faa'))
_REF_ATT = Constant(pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa'))
_PROTEINS = Constant(pkg_resources.resource_filename('explore_pipolin', '/data/HHpred_proteins.faa'))


_DEFAULT_FILTER = FilterTask()


def get_flow():
    with Flow('MAIN') as flow:
        genome = Parameter('genome')
        out_dir = Parameter('out_dir')
        add_colours = Parameter('add_colours')

        gquery = tasks.create_gquery.map(genome_file=genome)

        pipolb_blast_dir = tasks.run_blast_against_pipolb.map(gquery=gquery, ref_pipolb=unmapped(_REF_PIPOLB),
                                                              out_dir=unmapped(out_dir))
        gquery = tasks.add_features_from_blast.map(gquery=gquery, blast_dir=pipolb_blast_dir,
                                                   feature_type=unmapped(Constant(FeatureType.PIPOLB)))

        t_check_pipolbs = tasks.are_pipolbs_present.map(gquery=gquery)
        gquery = _DEFAULT_FILTER(tasks.return_result_if_true_else_none.map(result_to_filter=gquery,
                                                                           filter_by=t_check_pipolbs))

        att_blast_dir = tasks.run_blast_against_att.map(gquery=gquery, ref_att=unmapped(_REF_ATT),
                                                        out_dir=unmapped(out_dir))
        gquery = tasks.add_features_from_blast.map(gquery=gquery, blast_dir=att_blast_dir,
                                                   feature_type=unmapped(Constant(FeatureType.ATT)))

        aragorn_results_dir = tasks.detect_trnas_with_aragorn.map(gquery=gquery, out_dir=unmapped(out_dir))
        gquery = tasks.add_features_from_aragorn.map(gquery=gquery, aragorn_results_dir=aragorn_results_dir)

        atts_denovo = tasks.find_atts_denovo.map(gquery=gquery, out_dir=unmapped(out_dir))

        gquery = tasks.are_atts_present.map(gquery=gquery, upstream_tasks=[atts_denovo])

        gquery = tasks.analyse_pipolin_orientation.map(gquery=gquery)
        pipolin = tasks.scaffold_pipolins.map(gquery=gquery)

        pipolin_sequences = tasks.extract_pipolin_regions.map(gquery=gquery, pipolin=pipolin, out_dir=unmapped(out_dir))
        prokka = tasks.annotate_pipolins.map(gquery=gquery, pipolins_dir=pipolin_sequences,
                                             proteins=unmapped(_PROTEINS), out_dir=unmapped(out_dir))
        prokka_atts = tasks.include_atts_into_annotation.map(gquery=gquery, prokka_dir=prokka, pipolin=pipolin,
                                                             out_dir=unmapped(out_dir))
        with case(add_colours, True):
            tasks.easyfig_add_colours.map(gquery=gquery, in_dir=prokka_atts)

    flow.set_reference_tasks([prokka_atts])
    return flow
