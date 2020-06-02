import os
from prefect import task
from explore_pipolin.utilities import GQuery
from explore_pipolin.utilities import read_blastxml


@task
def add_features_from_blast(gquery: GQuery, blast_dir, feature_type):
    entries = read_blastxml(blast_xml=os.path.join(blast_dir, f'{gquery.gquery_id}.fmt5'))
    for entry in entries:
        for hit in entry:
            feature = gquery.feature_from_blasthit(hit=hit, contig_id=entry.id)
            gquery.get_features_by_type(feature_type).append(feature)
