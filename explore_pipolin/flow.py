from prefect import Flow, Parameter, case, unmapped, triggers
from prefect.tasks.control_flow import FilterTask

from explore_pipolin.tasks.prepare_for_the_analysis import prepare_for_the_analysis
from explore_pipolin.tasks.find_pipolbs import find_pipolbs, are_pipolbs_present, \
    continue_if_true_else_finished
from explore_pipolin.tasks.find_trnas import find_trnas
from explore_pipolin.tasks.find_atts import find_atts, find_atts_denovo, are_atts_present
from explore_pipolin.tasks.find_pipolins import find_pipolins
from explore_pipolin.tasks.reconstruct_pipolins import reconstruct_pipolins
from explore_pipolin.tasks.annotate_pipolins import save_pipolin_sequences, annotate_pipolins
from explore_pipolin.tasks.generate_results import generate_results

_DEFAULT_FILTER = FilterTask(trigger=triggers.all_successful)


def get_flow():
    with Flow('MAIN') as flow:
        genome_file = Parameter('genome_file')
        just_find_pipolbs = Parameter('just_find_pipolbs')
        no_annotation = Parameter('no_annotation')

        genome = prepare_for_the_analysis.map(original_file=genome_file)

        genome = find_pipolbs.map(genome=genome)

        t_check_pipolbs = are_pipolbs_present.map(genome=genome)
        genome = _DEFAULT_FILTER(continue_if_true_else_finished.map(
            result_to_filter=genome, filter_1=t_check_pipolbs, filter_2=unmapped(just_find_pipolbs))
        )

        genome = find_trnas.map(genome=genome)

        genome = find_atts.map(genome=genome)
        genome = find_atts_denovo.map(genome=genome)
        genome = are_atts_present.map(genome=genome)

        pipolins = find_pipolins.map(genome=genome)
        reconstructed_pipolins = reconstruct_pipolins.map(
            genome=genome, pipolins=pipolins
        )

        pipolin_seqs_dir = save_pipolin_sequences.map(
            genome=genome, pipolins=reconstructed_pipolins
        )

        with case(no_annotation, False):
            prokka_dir = annotate_pipolins.map(genome=genome, pipolins_dir=pipolin_seqs_dir)
            results_dir = generate_results.map(
                genome=genome, prokka_dir=prokka_dir, pipolins=reconstructed_pipolins
            )

    flow.set_reference_tasks([results_dir])
    return flow
