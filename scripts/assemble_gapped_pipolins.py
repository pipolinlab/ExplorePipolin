#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import read_genbank_records, write_genbank_records
from utilities import GenBankRecords
from utilities import check_dir


def get_atts_from_record(record):
    atts = []
    for i_f, feature in enumerate(record.features):
        if feature.type == 'repeat_region':
            atts.append((feature.location.start, feature.location.end, i_f))
    return atts


def get_polbs_from_record(record):
    polbs = []
    for i_f, feature in enumerate(record.features):
        if 'gene' in feature.qualifiers:
            if feature.qualifiers['gene'] == ['pi-polB']:
                polbs.append((feature.location.start, feature.location.end, i_f))
    return polbs


def get_unchangeable_contigs(record_set):
    """
    If there are a polB and an att, the whole contig is included into the assembly.
    """
    unchangeable_contigs = []
    for record_id, record in record_set.items():
        atts = get_atts_from_record(record)
        polbs = get_polbs_from_record(record)

        if len(atts) != 0 and len(polbs) != 0:
            unchangeable_contigs.append(record_id)
    return unchangeable_contigs


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
    assembly_gap_record = create_assembly_gap_record(record)
    atts = get_atts_from_record(record)
    polbs = get_polbs_from_record(record)

    modified_record = assembly_gap_record + record if atts[0][2] > polbs[0][2] else record + assembly_gap_record
    return modified_record


def get_att_only_contigs(record_set):
    att_only_contigs = []
    for record_id, record in record_set.items():
        atts = get_atts_from_record(record)
        polbs = get_polbs_from_record(record)

        if len(polbs) == 0 and len(atts) != 0:
            att_only_contigs.append(record_id)

    return att_only_contigs


def get_assembly_gap_position(unchangeable_record):
    for i_f, feature in enumerate(unchangeable_record.features):
        if feature.type == 'assembly_gap':
            return i_f


def glue_unchangeable_and_att(unchangeable_record, att_record):
    gap_position = get_assembly_gap_position(unchangeable_record)
    print(gap_position)
    glued_record = att_record + unchangeable_record if gap_position == 0 else unchangeable_record + att_record
    return glued_record


def create_single_record(record_set, pipolin_features):
    unchangeable_contigs = get_unchangeable_contigs(record_set)
    print(f'The unchangeable contigs: {unchangeable_contigs}!')
    for contig in unchangeable_contigs:
        record_set[contig] = add_assembly_gap_to_unchangeable(record_set[contig])

    if len(unchangeable_contigs) == 1:
        att_only_contigs = get_att_only_contigs(record_set)
        if len(att_only_contigs) == 1:
            single_record = glue_unchangeable_and_att(record_set[unchangeable_contigs[0]],
                                                     record_set[att_only_contigs[0]])
        else:
            # TODO: traverse the features to infer the order
    elif len(unchangeable_contigs) == 0:
        # TODO: not trivial here...

    # return single_record


def assemble_gapped_pipolins(gb_records: GenBankRecords, pipolin_features):
    for strain_id, record_set in gb_records.items():
        if len(record_set) > 1:
            print(f'Assembling pipolin region for {strain_id}...')
            gb_records[strain_id] = create_single_record(record_set, pipolin_features)


def get_unique_pipolin_features(gb_records: GenBankRecords):
    """
    Create a set of the unique pipolin features (excluding hypothetical proteins)
    grouped by a feature type. Only ungapped pipolins will be traversed for features.
    """
    pipolin_features = {'tRNA': set(), 'CDS': set()}
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
    assemble_gapped_pipolins(gb_records, pipolin_features)
    
    check_dir(out_dir)
    write_genbank_records(gb_records, out_dir)


if __name__ == '__main__':
    main()
