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

        genome = tasks.create_genome.map(genome_file=genome)

        pipolb_blast_dir = tasks.run_blast_against_pipolb.map(genome=genome, ref_pipolb=unmapped(_REF_PIPOLB),
                                                              out_dir=unmapped(out_dir))
        genome = tasks.add_features_from_blast.map(genome=genome, blast_dir=pipolb_blast_dir,
                                                   feature_type=unmapped(Constant(FeatureType.PIPOLB)))

        t_check_pipolbs = tasks.are_pipolbs_present.map(genome=genome)
        genome = _DEFAULT_FILTER(tasks.return_result_if_true_else_none.map(
            result_to_filter=genome, filter_by=t_check_pipolbs)
        )

        att_blast_dir = tasks.run_blast_against_att.map(genome=genome, ref_att=unmapped(_REF_ATT),
                                                        out_dir=unmapped(out_dir))
        genome = tasks.add_features_from_blast.map(genome=genome, blast_dir=att_blast_dir,
                                                   feature_type=unmapped(Constant(FeatureType.ATT)))

        aragorn_results_dir = tasks.detect_trnas_with_aragorn.map(genome=genome, out_dir=unmapped(out_dir))
        genome = tasks.add_features_from_aragorn.map(genome=genome, aragorn_results_dir=aragorn_results_dir)

        atts_denovo = tasks.find_atts_denovo.map(genome=genome, out_dir=unmapped(out_dir))

        genome = tasks.are_atts_present.map(genome=genome, upstream_tasks=[atts_denovo])

        genome = tasks.analyse_pipolin_orientation.map(genome=genome)
        pipolin = tasks.scaffold_pipolins.map(genome=genome)

        pipolin_sequences = tasks.extract_pipolin_regions.map(genome=genome, pipolin=pipolin,
                                                              out_dir=unmapped(out_dir))
        prokka = tasks.annotate_pipolins.map(genome=genome, pipolins_dir=pipolin_sequences,
                                             proteins=unmapped(_PROTEINS), out_dir=unmapped(out_dir))
        prokka_atts = tasks.include_atts_into_annotation.map(genome=genome, prokka_dir=prokka,
                                                             pipolin=pipolin, out_dir=unmapped(out_dir))
        with case(add_colours, True):
            tasks.easyfig_add_colours.map(genome=genome, in_dir=prokka_atts)

    flow.set_reference_tasks([prokka_atts])
    return flow
