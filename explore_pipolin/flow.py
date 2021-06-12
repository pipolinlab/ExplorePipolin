from prefect import Flow, Parameter, unmapped, case
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
from explore_pipolin.tasks.easyfig_coloring import easyfig_add_colours

_DEFAULT_FILTER = FilterTask()


def get_flow():
    with Flow('MAIN') as flow:
        genome_file = Parameter('genome_file')
        out_dir = Parameter('out_dir')
        pipolb_hmm_profile = Parameter('pipolb_hmm_profile')
        ref_att = Parameter('ref_att')
        percent_identity = Parameter('percent_identity')
        max_inflate = Parameter('max_inflate')
        no_annotation = Parameter('no_annotation')
        proteins = Parameter('proteins')
        skip_colours = Parameter('skip_colours')
        cpus = Parameter('cpus')
        do_not_reuse = Parameter('do_not_reuse')

        genome = prepare_for_the_analysis.map(original_file=genome_file, out_dir=unmapped(out_dir))

        genome = find_pipolbs.map(
            genome=genome,
            pipolb_hmm_profile=unmapped(pipolb_hmm_profile),
            do_not_reuse=unmapped(do_not_reuse),
        )

        t_check_pipolbs = are_pipolbs_present.map(genome=genome)
        genome = _DEFAULT_FILTER(continue_if_true_else_finished.map(
            result_to_filter=genome, filter_by=t_check_pipolbs)
        )

        genome = find_trnas.map(
            genome=genome, do_not_reuse=unmapped(do_not_reuse)
        )

        genome = find_atts.map(
            genome=genome,
            ref_att=unmapped(ref_att),
            do_not_reuse=unmapped(do_not_reuse),
        )
        genome = find_atts_denovo.map(
            genome=genome,
            percent_identity=unmapped(percent_identity),
            do_not_reuse=unmapped(do_not_reuse),
        )
        genome = are_atts_present.map(genome=genome)

        pipolins = find_pipolins.map(genome=genome)
        reconstructed_pipolins = reconstruct_pipolins.map(
            genome=genome, pipolins=pipolins, no_border_inflate=unmapped(max_inflate)
        )

        pipolin_seqs_dir = save_pipolin_sequences.map(
            genome=genome, pipolins=reconstructed_pipolins
        )

        with case(no_annotation, False):
            prokka_dir = annotate_pipolins.map(
                genome=genome,
                pipolins_dir=pipolin_seqs_dir,
                proteins=unmapped(proteins),
                cpus=unmapped(cpus)
            )
            results_dir = generate_results.map(
                genome=genome, prokka_dir=prokka_dir, pipolins=reconstructed_pipolins
            )

            with case(skip_colours, False):
                easyfig_add_colours.map(genome=genome, in_dir=results_dir)

    flow.set_reference_tasks([results_dir])
    return flow
