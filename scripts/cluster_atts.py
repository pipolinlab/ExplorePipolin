#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('msa', type=click.Path(exists=True))
@click.argument('sequences', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path(exists=True))
def cluster_msa(msa, sequences, out_dir):
    """
    TODO
    """
    alignment = AlignIO.read(msa, 'fasta')
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'upgma')
    tree = constructor.build_tree(alignment)
    # tree.ladderize()
    # Phylo.draw(tree)

    # TODO: think if it possible to get rid of hard-coded clade names ?
    inner_156 = tree.find_clades(name='Inner156', order='level')
    subtree1 = tree.from_clade(list(inner_156)[0])
    terminals1 = [i.name for i in subtree1.get_terminals()]
    terminals1.sort()

    inner_155 = tree.find_clades(name='Inner155', order='level')
    subtree2 = tree.from_clade(list(inner_155)[0])
    terminals2 = [i.name for i in subtree2.get_terminals()]
    terminals2.sort()

    seqs = SeqIO.to_dict(SeqIO.parse(sequences, 'fasta'))
    seqs1 = [seqs[name] for name in terminals1]
    seqs2 = [seqs[name] for name in terminals2]
    seqs2_rev = [
        SeqRecord(seq=i.seq.reverse_complement(), id=i.id, name=i.name, description=i.description) for i in seqs2
    ]
    unified_seqs = seqs1 + seqs2_rev
    with open(os.path.join(out_dir, 'att_sequences.unified.fa'), 'w') as ouf:
        SeqIO.write(seqs1, ouf, 'fasta')


if __name__ == '__main__':
    cluster_msa()
