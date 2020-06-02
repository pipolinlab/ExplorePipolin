import os
from prefect import task
from explore_pipolin.utilities import GQuery, Feature, Orientation


@task
def add_features_atts_denovo(gquery: GQuery, atts_denovo_dir):
    with open(os.path.join(atts_denovo_dir, gquery.gquery_id + '.atts')) as inf:
        _ = inf.readline()
        for line in inf:
            att_pair = [int(i) for i in line.strip().split(sep='\t')]
            gquery.denovo_atts.append(Feature(start=att_pair[0], end=att_pair[1],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))
            gquery.denovo_atts.append(Feature(start=att_pair[2], end=att_pair[3],
                                              frame=Orientation.FORWARD, contig=gquery.contigs[0]))

    for att in gquery.denovo_atts:
        target_trna = gquery.find_target_trna(att)
        if target_trna is not None:
            gquery.target_trnas_denovo.append(target_trna)
