import os
from prefect import task
from explore_pipolin.utilities import GQuery, Feature
from explore_pipolin.utilities import read_aragorn_batch


@task
def add_features_from_aragorn(gquery: GQuery, aragorn_dir):
    entries = read_aragorn_batch(aragorn_batch=os.path.join(aragorn_dir, f'{gquery.gquery_id}.batch'))
    for contig_id, hits in entries.items():
        for hit in hits:
            feature = Feature(start=hit[0], end=hit[1], frame=hit[2], contig=gquery.get_contig_by_id(contig_id))
            gquery.trnas.append(feature)

    for att in gquery.atts:
        target_trna = gquery.find_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas.append(target_trna)
