#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
from collections import defaultdict
from typing import Any

import click
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import read_genbank_records, write_genbank_records
from utilities import GenBankRecords
from utilities import check_dir


def get_unchangeable_contig(record_set):
    """
    If there are a polB and an att, as well as there is only polB on the contig,
    the whole contig should be included into the assembly.
    """
    unchangeable_contig = set()

    for record_id, record in record_set.items():
        atts_polbs = []
        for feature in record.features:
            if feature.type == 'repeat_region':
                atts_polbs.append('att')
            if 'gene' in feature.qualifiers:
                if feature.qualifiers['gene'] == ['pi-polB']:
                    atts_polbs.append('polB')
        if 'att' in atts_polbs and 'polB' in atts_polbs:
            unchangeable_contig.add(record_id)
        elif 'polB' in atts_polbs:
            unchangeable_contig.add(record_id)

    return unchangeable_contig


def create_assembly_gap_record(record):
    assembly_gap_seq = Seq('N' * 100, alphabet=IUPACAmbiguousDNA())
    assembly_gap_qualifiers = {'estimated_length': ['unknown'],
                               'gap_type': ['within_scaffolds'],
                               'linkage_evidence': ['pipolin_structure']}
    assembly_gap_feature = SeqFeature(type='assembly_gap',
                                      location=FeatureLocation(1, 100, strand=+1),
                                      qualifiers=assembly_gap_qualifiers)
    assembly_gap_record = SeqRecord(seq=assembly_gap_seq, id=record.id, name=record.name,
                                    description=record.description, features=[assembly_gap_feature],
                                    annotations=record.annotations)

    return assembly_gap_record


def add_assembly_gap_to_unchangeable(record):

    att_index = []
    polB_index = []
    for i_f, feature in enumerate(record.features):
        if feature.type == 'repeat_region':
            att_index.append(i_f)
        if 'gene' in feature.qualifiers:
            if feature.qualifiers['gene'] == ['pi-polB']:
                polB_index.append(i_f)

    if len(att_index) > 1:
        raise AssertionError('Only one att near the pi-polB is expected!')
    if len(polB_index) > 1:
        raise AssertionError('Only one pi-polB per record is expected!')

    if len(att_index) != 0:
        assembly_gap_record = create_assembly_gap_record(record)
        modified_record = assembly_gap_record + record if att_index[0] > polB_index[0] \
            else record + assembly_gap_record
    # TODO: continue from here next time!
    return modified_record


def create_single_record(record_set):
    unchangeable_contig = get_unchangeable_contig(record_set)
    if len(unchangeable_contig) > 1:
        raise AssertionError('Only a single pipolin region per genome is expected, but got more!')
    elif len(unchangeable_contig) == 0:
        raise AssertionError('At least one unchangeable contig is expected, but got 0!')
    unchangeable_contig = unchangeable_contig.pop()
    record_set[unchangeable_contig] = add_assembly_gap_to_unchangeable(record_set[unchangeable_contig])
    print(record_set[unchangeable_contig])


    # return single_record


def assemble_gapped_pipolins(gb_records: GenBankRecords):
    for strain_id, record_set in gb_records.items():
        if len(record_set) > 1:
            gb_records[strain_id] = create_single_record(record_set)


def get_unique_pipolin_features(gb_records: GenBankRecords):
    """
    Create a set of the unique pipolin features (excluding the hypothetical proteins)
    grouped by a feature type. Only ungapped pipolins will be traversed for features.
    """
    pipolin_features = {'repeat_region': set(), 'tRNA': set(), 'CDS': set()}
    for record_set in gb_records.values():
        if len(record_set) == 1:
            for record in record_set.values():
                for feature in record.features:
                    if feature.type in pipolin_features:
                        if 'product' in feature.qualifiers:
                            if len(feature.qualifiers['product']) > 1:
                                raise AssertionError('Only a single product is expected!')
                            if feature.qualifiers['product'][0] != 'hypothetical protein':
                                pipolin_features[feature.type].add(feature.qualifiers['product'][0])

    return pipolin_features


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def main(in_dir, out_dir):
    """
    The script takes IN_DIR with *.gbk files (generated after annotation and including atts),
    detects not assembled pipolins and tries to order the contigs into one sequence,
    filling in gaps with NNs and adding the /assembly_gap feature to the *.gbk files.
    """
    gb_records = read_genbank_records(in_dir)
    pipolin_features = get_unique_pipolin_features(gb_records)
    assemble_gapped_pipolins(gb_records)
    
    check_dir(out_dir)
    write_genbank_records(gb_records, out_dir)


if __name__ == '__main__':
    main()
