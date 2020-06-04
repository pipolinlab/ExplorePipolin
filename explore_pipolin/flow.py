import pkg_resources
from prefect import Flow, Parameter, unmapped
from prefect.tasks.core.constants import Constant

from explore_pipolin import tasks
from explore_pipolin.utilities import FeatureType

_REF_POLB = Constant(pkg_resources.resource_filename('explore_pipolin', 'data/pi-polB.faa'))
_REF_ATT = Constant(pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa'))
_PROTEINS = Constant(pkg_resources.resource_filename('explore_pipolin', '/data/HHpred_proteins.faa'))


def get_flow():
    with Flow('MAIN') as flow:
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')
        abricate_dir = Parameter('abricate_dir')

        gquery = tasks.create_gquery.map(genome=genomes)

        polbs_blast = tasks.run_blast_against_polb.map(genome=genomes, root_dir=unmapped(out_dir),
                                                       reference=unmapped(_REF_POLB))
        t_add_polbs = tasks.add_features_from_blast.map(gquery=gquery, blast_dir=polbs_blast,
                                                        feature_type=unmapped(Constant(FeatureType.POLB)))

        atts_blast = tasks.run_blast_against_att.map(genome=genomes, root_dir=unmapped(out_dir),
                                                     reference=unmapped(_REF_ATT))
        t_add_atts = tasks.add_features_from_blast.map(gquery=gquery, blast_dir=atts_blast,
                                                       feature_type=unmapped(Constant(FeatureType.ATT)))

        aragorn_results = tasks.detect_trnas_with_aragorn.map(genome=genomes, root_dir=unmapped(out_dir))
        t_add_trnas = tasks.add_features_from_aragorn.map(gquery=gquery, aragorn_dir=aragorn_results,
                                                          upstream_tasks=[t_add_polbs, t_add_atts])

        t_check_polbs = tasks.are_polbs_present.map(gquery=gquery, upstream_tasks=[t_add_trnas])

        atts_denovo = tasks.find_atts_denovo.map(genome=genomes, gquery=gquery, root_dir=unmapped(out_dir),
                                                 upstream_tasks=[t_check_polbs])
        t_add_denovo_atts = tasks.add_features_atts_denovo.map(gquery=gquery, atts_denovo_dir=atts_denovo)

        t_check_features = tasks.are_atts_present.map(gquery=gquery, upstream_tasks=[t_add_denovo_atts])

        t_set_orientations = tasks.analyse_pipolin_orientation.map(gquery=gquery, upstream_tasks=[t_check_features])
        t_scaffolding = tasks.scaffold_pipolins.map(gquery=gquery, upstream_tasks=[t_set_orientations])

        pipolin_sequences = tasks.extract_pipolin_regions.map(genome=genomes, gquery=gquery,
                                                              root_dir=unmapped(out_dir),
                                                              upstream_tasks=[t_scaffolding])
        prokka = tasks.annotate_pipolins.map(gquery=gquery, pipolins_dir=pipolin_sequences,
                                             proteins=unmapped(_PROTEINS), root_dir=unmapped(out_dir))
        prokka_atts = tasks.include_atts_into_annotation.map(gquery=gquery, prokka_dir=prokka,
                                                             root_dir=unmapped(out_dir))
        # TODO: make this task optional
        t_easyfig = tasks.easyfig_add_colours.map(gquery=gquery, in_dir=prokka_atts,
                                                  abricate_dir=unmapped(abricate_dir),
                                                  upstream_tasks=[t_scaffolding])
        # TODO: when easyfig skipped, set (skip_on_upstream_skip=False)
        # prokka_atts_positions = tasks.set_correct_positions.map(gquery=gquery, prokka_atts_dir=prokka_atts,
        #                                                         root_dir=unmapped(out_dir),
        #                                                         upstream_tasks=[t_easyfig])

    return flow
