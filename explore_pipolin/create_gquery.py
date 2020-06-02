from prefect import task
from explore_pipolin.utilities import GQuery, Contig
from explore_pipolin.utilities import define_gquery_id
from explore_pipolin.utilities import read_seqio_records


@task
def create_gquery(genome) -> GQuery:
    gquery = GQuery(gquery_id=define_gquery_id(genome=genome))
    genome_dict = read_seqio_records(file=genome, file_format='fasta')
    for key, value in genome_dict.items():
        contig = Contig(contig_id=key, contig_length=len(value.seq))
        gquery.contigs.append(contig)

    return gquery
