from prefect import Flow, Parameter, unmapped, case
from prefect.tasks.control_flow import FilterTask

from explore_pipolin.tasks_related.atts_search import find_atts, find_atts_denovo, are_atts_present
from explore_pipolin.tasks_related.find_pipolins import find_pipolins
from explore_pipolin.tasks_related.scaffolding import refine_pipolins
from explore_pipolin import tasks

_DEFAULT_FILTER = FilterTask()


def get_flow():
    with Flow('MAIN') as flow:
        genome_file = Parameter('genome_file')
        out_dir = Parameter('out_dir')
        add_colours = Parameter('add_colours')

        genome = tasks.create_genome.map(genome_file=genome_file)

        genome = tasks.find_pipolbs.map(genome=genome, out_dir=unmapped(out_dir))

        t_check_pipolbs = tasks.are_pipolbs_present.map(genome=genome)
        genome = _DEFAULT_FILTER(tasks.return_result_if_true_else_none.map(
            result_to_filter=genome, filter_by=t_check_pipolbs)
        )

        genome = tasks.find_trnas.map(genome=genome, out_dir=unmapped(out_dir))

        genome = find_atts.map(genome=genome, out_dir=unmapped(out_dir))
        genome = find_atts_denovo.map(genome=genome, out_dir=unmapped(out_dir))

        genome = are_atts_present.map(genome=genome)

        pipolins = find_pipolins.map(genome=genome)
        pipolins = refine_pipolins.map(genome=genome, pipolins=pipolins)

        pipolin_sequences = tasks.extract_pipolin_regions.map(genome=genome, pipolins=pipolins,
                                                              out_dir=unmapped(out_dir))
        prokka = tasks.annotate_pipolins.map(genome=genome, pipolins_dir=pipolin_sequences,
                                             out_dir=unmapped(out_dir))
        prokka_atts = tasks.include_atts.map(genome=genome, prokka_dir=prokka,
                                             pipolins=pipolins, out_dir=unmapped(out_dir))
        with case(add_colours, True):
            tasks.easyfig_add_colours.map(genome=genome, in_dir=prokka_atts)

    flow.set_reference_tasks([prokka_atts])
    return flow
