from prefect import Flow, Parameter, unmapped, case
from prefect.tasks.control_flow import FilterTask

from explore_pipolin.tasks.misc import create_genome
from explore_pipolin.tasks.find_pipolbs import find_pipolbs, are_pipolbs_present, \
    return_result_if_true_else_none
from explore_pipolin.tasks.find_trnas import find_trnas
from explore_pipolin.tasks.find_atts import find_atts, find_atts_denovo, are_atts_present
from explore_pipolin.tasks.find_pipolins import find_pipolins
from explore_pipolin.tasks.scaffold_pipolins import scaffold_pipolins
from explore_pipolin.tasks.annotate_pipolins import save_pipolin_sequences, annotate_pipolins
from explore_pipolin.tasks.include_atts import include_atts
from explore_pipolin.tasks.easyfig_coloring import easyfig_add_colours

_DEFAULT_FILTER = FilterTask()


def get_flow():
    with Flow('MAIN') as flow:
        genome_file = Parameter('genome_file')
        out_dir = Parameter('out_dir')
        add_colours = Parameter('add_colours')

        genome = create_genome.map(genome_file=genome_file, output_dir=unmapped(out_dir))

        genome = find_pipolbs.map(genome=genome)

        t_check_pipolbs = are_pipolbs_present.map(genome=genome)
        genome = _DEFAULT_FILTER(return_result_if_true_else_none.map(
            result_to_filter=genome, filter_by=t_check_pipolbs)
        )

        genome = find_trnas.map(genome=genome)

        genome = find_atts.map(genome=genome)
        genome = find_atts_denovo.map(genome=genome)

        genome = are_atts_present.map(genome=genome)

        pipolins = find_pipolins.map(genome=genome)
        scaffolded_pipolins = scaffold_pipolins.map(genome=genome, pipolins=pipolins)

        pipolin_seqs_dir = save_pipolin_sequences.map(genome=genome, pipolins=scaffolded_pipolins)
        prokka_dir = annotate_pipolins.map(genome=genome, pipolins_dir=pipolin_seqs_dir)
        results_dir = include_atts.map(genome=genome, prokka_dir=prokka_dir, pipolins=pipolins)
        with case(add_colours, True):
            easyfig_add_colours.map(genome=genome, in_dir=results_dir)

    flow.set_reference_tasks([results_dir])
    return flow
