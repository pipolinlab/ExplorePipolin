#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

from typing import Mapping, MutableMapping
from itertools import *

import csv
import click
from glob import glob
import os.path
from BCBio import GFF
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from explore_pipolin.utilities import CONTEXT_SETTINGS, get_roary_groups


class FeatureSet:
    def __init__(self, seq_record: SeqRecord):
        self.seq_record = seq_record
        self._features = {feature.id : feature for feature in seq_record.features}

    def __getitem__(self, item) -> SeqFeature:
        return self._features[item]

    def get(self, key, default=None) -> SeqFeature:
        return self._features.get(key, default)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('roary-dir', type=click.Path(exists=True))
@click.argument('annot-dir', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(roary_dir, annot_dir, output):
    roary_groups = get_roary_groups(roary_dir)
    gffs = read_gff_features(annot_dir)
    entries = []
    for group_name, genome in islice(roary_groups.items(), 10):
        for genome_name, genes in genome.items():
            for gene in genes:
                feature = gffs[genome_name].get(gene)
                entries.append({
                    'genome': genome_name,
                    'gene': group_name,
                    'start': feature.location.start,
                    'end': feature.location.end})

    with open(output, 'w') as outf:
        writer = csv.DictWriter(outf, ('genome', 'gene', 'start', 'end'))
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)


def read_gffs(annot_dir) -> Mapping[str, SeqRecord]:
    gff_records: MutableMapping[str, SeqRecord] = {}
    for file in glob(os.path.join(annot_dir, "*.[13].gff")):
        with open(file) as inf:
            entry, = GFF.parse(inf)
            gff_records[entry.id] = entry

            # gff_records[os.path.splitext(file)[0]] = records
    return gff_records


def read_gff_features(annot_dir) -> Mapping[str, FeatureSet]:
    return {name: FeatureSet(sr) for name, sr in read_gffs(annot_dir).items()}


if __name__ == '__main__':
    main()
