#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from prefect import task
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import PipolinFragment, GQuery

# Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
# https://github.com/biopython/biopython/issues/1755


def create_assembly_gap_record(record):
    source_feature = SeqFeature(type='source', location=FeatureLocation(1, 100, strand=+1),
                                qualifiers={'mol_type': record.features[0].qualifiers['mol_type'],
                                            'organism': record.features[0].qualifiers['organism'],
                                            'strain': record.features[0].qualifiers['strain']})
    assembly_gap_seq = Seq('N' * 100, alphabet=IUPACAmbiguousDNA())
    assembly_gap_qualifiers = {'estimated_length': ['unknown'],
                               'gap_type': ['within_scaffolds'],
                               'linkage_evidence': ['pipolin_structure']}
    assembly_gap_feature = SeqFeature(type='assembly_gap',
                                      location=FeatureLocation(1, 100, strand=+1),
                                      qualifiers=assembly_gap_qualifiers)
    assembly_gap_record = SeqRecord(seq=assembly_gap_seq, id=record.id, name=record.name,
                                    description=record.description, features=[source_feature, assembly_gap_feature],
                                    annotations=record.annotations)

    return assembly_gap_record


@task
def scaffold_pipolins(gquery: GQuery):
    if gquery.is_single_contig() or gquery.is_on_the_same_contig():
        print('>>>Scaffolding is not required!')
        start, end = gquery.get_pipolin_bounds()
        pipolin = PipolinFragment(contig=gquery.get_contig_by_id(gquery.polbs[0].contig.contig_id),
                                  start=start, end=end)

        pipolin.atts.extend(gquery.atts)
        gquery.pipolin_fragments.append(pipolin)
    else:
        print('>>>Scaffolding is required!')
        gquery.try_creating_single_record()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
@click.option('--long', is_flag=True, help='If long, pipolin regions are identified by leftmost and rightmost '
                                           'atts. If it is not specified, the closest atts to pi-polB will determine '
                                           'the pipolin bounds.')
def main(in_dir, out_dir, long):
    """
    The script takes IN_DIR with *.gbk files (generated after annotation and with included atts),
    detects not assembled pipolins and tries to order the contigs into one sequence,
    filling in gaps with NNs and adding the /assembly_gap feature to the *.gbk files.
    """
    scaffold_pipolins(in_dir, out_dir, long)


if __name__ == '__main__':
    main()
