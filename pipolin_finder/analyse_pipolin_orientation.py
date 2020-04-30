from prefect import task


@task
def analyse_pipolin_orientation(gquery):
    gquery.is_single_target_trna_per_contig()
    for contig in gquery.contigs:
        contig.contig_orientation = gquery.get_contig_orientation(contig)
