#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from explore_pipolin.utilities import CONTEXT_SETTINGS
from explore_pipolin.utilities import Feature, GQuery
from explore_pipolin.utilities import read_from_shelve
from explore_pipolin.utilities import save_to_shelve


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


def store_new_att_bounds(shelve_in_dir, object_name, att_hmm_dir):
    pipolins = read_from_shelve(os.path.join(shelve_in_dir, 'shelve.db'), 'pipolins')
    for pipolin in pipolins:
        pipolin.atts = []
        hit_names, hit_bounds = parse_hmmer_tbl(os.path.join(att_hmm_dir, f'{pipolin.strain_id}-atts.tbl'))
        for name, bounds in zip(hit_names, hit_bounds):
            pipolin.atts.append(Feature(bounds[0], bounds[1], name))
    save_to_shelve(os.path.join(shelve_in_dir, 'shelve.db'), pipolins, object_name)

    return object_name


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('att-hmm-dir', type=click.Path(exists=True))
@click.option('--object-name', required=True)
def main(shelve_file, object_name, att_hmm_dir):
    """
    Creates a new shelve object ("short-pipolins" or "long-pipolins") with att bounds, defined by HMMER.
    """
    # TODO: This is not a proper way to store new att positions, because pi-polB positions are not updated!
    store_new_att_bounds(shelve_file, object_name, att_hmm_dir)


if __name__ == '__main__':
    main()
