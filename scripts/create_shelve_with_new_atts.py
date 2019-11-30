#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from utilities import CONTEXT_SETTINGS
from utilities import Feature, Pipolin
from utilities import read_from_shelve
from utilities import save_to_shelve


def parse_hmmer_tbl(file):
    hit_names = []
    hit_bounds = []
    with open(file) as inf:
        for line in inf:
            if line[0] != '#':
                entry = line.strip().split()
                hit_names.append(entry[2])
                hit_bounds.append((int(entry[6]), int(entry[7])))
    return hit_names, hit_bounds


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('old-shelve-file', type=click.Path(exists=True))
@click.argument('att-hmm-dir', type=click.Path(exists=True))
@click.argument('out-file')
def create_shelve_with_new_atts(old_shelve_file, att_hmm_dir, out_file):
    """
    Creates a new shelve file with pipolins, containing new att bounds defined by HMMER.
    """
    pipolins = read_from_shelve(old_shelve_file, 'pipolins')

    for pipolin in pipolins:
        pipolin.atts = []
        hit_names, hit_bounds = parse_hmmer_tbl(os.path.join(att_hmm_dir, f'{pipolin.strain_id}-atts.tbl'))
        for name, bounds in zip(hit_names, hit_bounds):
            pipolin.atts.append(Feature(bounds[0], bounds[1], name))

    save_to_shelve(out_file, pipolins, 'pipolins')


if __name__ == '__main__':
    create_shelve_with_new_atts()
