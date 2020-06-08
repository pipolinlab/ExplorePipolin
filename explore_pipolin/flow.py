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
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')
        add_colours = Parameter('add_colours')

        gquery = tasks.create_gquery.map(genome=genomes)

        pipolb_blast_dir = tasks.run_blast_against_pipolb.map(genome=genomes, out_dir=unmapped(out_dir),
                                                              ref_pipolb=unmapped(_REF_PIPOLB))
        t_add_pipolbs = tasks.add_features_from_blast.map(gquery=gquery, blast_dir=pipolb_blast_dir,
                                                          feature_type=unmapped(Constant(FeatureType.PIPOLB)))

        att_blast_dir = tasks.run_blast_against_att.map(genome=genomes, out_dir=unmapped(out_dir),
                                                        ref_att=unmapped(_REF_ATT))
        t_add_atts = tasks.add_features_from_blast.map(gquery=gquery, blast_dir=att_blast_dir,
                                                       feature_type=unmapped(Constant(FeatureType.ATT)))

        aragorn_results_dir = tasks.detect_trnas_with_aragorn.map(genome=genomes, out_dir=unmapped(out_dir))
        t_add_trnas = tasks.add_features_from_aragorn.map(gquery=gquery, aragorn_results_dir=aragorn_results_dir,
                                                          upstream_tasks=[t_add_pipolbs, t_add_atts])
        # TODO: move it before att blast and aragorn!
        t_check_pipolbs = tasks.are_pipolbs_present.map(gquery=gquery, upstream_tasks=[t_add_trnas])

        gquery = _DEFAULT_FILTER(tasks.filter_no_pipolbs.map(task_to_filter=gquery, filter_by=t_check_pipolbs))
        genomes = _DEFAULT_FILTER(tasks.filter_no_pipolbs.map(task_to_filter=genomes, filter_by=t_check_pipolbs))

        atts_denovo = tasks.find_atts_denovo.map(genome=genomes, gquery=gquery, out_dir=unmapped(out_dir))
        t_add_denovo_atts = tasks.add_features_atts_denovo.map(gquery=gquery, atts_denovo_dir=atts_denovo)

        t_check_atts = tasks.are_atts_present.map(gquery=gquery, upstream_tasks=[t_add_denovo_atts])

        t_set_orientations = tasks.analyse_pipolin_orientation.map(gquery=gquery, upstream_tasks=[t_check_atts])
        t_scaffolding = tasks.scaffold_pipolins.map(gquery=gquery, upstream_tasks=[t_set_orientations])

        pipolin_sequences = tasks.extract_pipolin_regions.map(genome=genomes, gquery=gquery,
                                                              out_dir=unmapped(out_dir),
                                                              upstream_tasks=[t_scaffolding])
        prokka = tasks.annotate_pipolins.map(gquery=gquery, pipolins_dir=pipolin_sequences,
                                             proteins=unmapped(_PROTEINS), out_dir=unmapped(out_dir))
        prokka_atts = tasks.include_atts_into_annotation.map(gquery=gquery, prokka_dir=prokka,
                                                             out_dir=unmapped(out_dir))
        with case(add_colours, True):
            tasks.easyfig_add_colours.map(gquery=gquery, in_dir=prokka_atts)

    return flow
