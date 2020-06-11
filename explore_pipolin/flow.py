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

        features_container = tasks.create_features_container.map(genome_file=genome)

        pipolb_blast_dir = tasks.run_blast_against_pipolb.map(features_container=features_container,
                                                              ref_pipolb=unmapped(_REF_PIPOLB),
                                                              out_dir=unmapped(out_dir))
        features_container = tasks.add_features_from_blast.map(features_container=features_container,
                                                               blast_dir=pipolb_blast_dir,
                                                               feature_type=unmapped(Constant(FeatureType.PIPOLB)))

        t_check_pipolbs = tasks.are_pipolbs_present.map(features_container=features_container)
        features_container = _DEFAULT_FILTER(tasks.return_result_if_true_else_none.map(
            result_to_filter=features_container, filter_by=t_check_pipolbs))

        att_blast_dir = tasks.run_blast_against_att.map(features_container=features_container,
                                                        ref_att=unmapped(_REF_ATT),
                                                        out_dir=unmapped(out_dir))
        features_container = tasks.add_features_from_blast.map(features_container=features_container,
                                                               blast_dir=att_blast_dir,
                                                               feature_type=unmapped(Constant(FeatureType.ATT)))

        aragorn_results_dir = tasks.detect_trnas_with_aragorn.map(features_container=features_container,
                                                                  out_dir=unmapped(out_dir))
        features_container = tasks.add_features_from_aragorn.map(features_container=features_container,
                                                                 aragorn_results_dir=aragorn_results_dir)

        atts_denovo = tasks.find_atts_denovo.map(features_container=features_container, out_dir=unmapped(out_dir))

        features_container = tasks.are_atts_present.map(features_container=features_container,
                                                        upstream_tasks=[atts_denovo])

        features_container = tasks.analyse_pipolin_orientation.map(features_container=features_container)
        pipolin = tasks.scaffold_pipolins.map(features_container=features_container)

        pipolin_sequences = tasks.extract_pipolin_regions.map(features_container=features_container, pipolin=pipolin,
                                                              out_dir=unmapped(out_dir))
        prokka = tasks.annotate_pipolins.map(features_container=features_container, pipolins_dir=pipolin_sequences,
                                             proteins=unmapped(_PROTEINS), out_dir=unmapped(out_dir))
        prokka_atts = tasks.include_atts_into_annotation.map(features_container=features_container, prokka_dir=prokka,
                                                             pipolin=pipolin, out_dir=unmapped(out_dir))
        with case(add_colours, True):
            tasks.easyfig_add_colours.map(features_container=features_container, in_dir=prokka_atts)

    flow.set_reference_tasks([prokka_atts])
    return flow
