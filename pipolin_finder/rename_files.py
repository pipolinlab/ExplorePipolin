#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import shutil
import click
from pipolin_finder.utilities import CONTEXT_SETTINGS


def create_dict_w_strainnames(names_file):
    accs_to_strainnames = {}
    with open(names_file) as inf:
        for line in inf:
            names = line.strip().split(sep='\t')
            accs_to_strainnames[names[1]] = names[0]
    return accs_to_strainnames


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in_dir', type=click.Path(exists=True))
@click.argument('out_dir', type=click.Path())
@click.argument('names_file', type=click.Path(exists=True))
def rename(in_dir, out_dir, names_file):
    """
    TODO
    """
    accs_to_strainnames = create_dict_w_strainnames(names_file)
    os.makedirs(out_dir, exist_ok=True)
    for root, _, files in os.walk(in_dir):
        for filename in files:
            old_path = os.path.join(root, filename)
            base, extension = os.path.splitext(filename)
            try:
                new_filename = accs_to_strainnames[base] + extension
            except KeyError:
                new_filename = base + extension

            new_path = os.path.join(out_dir, new_filename)
            shutil.copy(old_path, new_path)


if __name__ == '__main__':
    rename()