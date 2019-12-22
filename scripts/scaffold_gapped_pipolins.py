#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import read_genbank_records, write_genbank_records
from utilities import write_gff_records
from utilities import GenBankRecords
from utilities import check_dir

# Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
# https://github.com/biopython/biopython/issues/1755


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


def add_assembly_gap_to_unchangeable(record):
    assembly_gap_record = create_assembly_gap_record(record)
    atts = get_atts_from_record(record)
    polbs = get_polbs_from_record(record)

    modified_record = assembly_gap_record + record if atts[-1][2] > polbs[-1][2] else record + assembly_gap_record
    return modified_record


def get_att_only_contigs(record_set):
    att_only_contigs = []
    for record_id, record in record_set.items():
        atts = get_atts_from_record(record)
        polbs = get_polbs_from_record(record)

        if len(polbs) == 0 and len(atts) != 0:
            att_only_contigs.append(record_id)

    return att_only_contigs


def get_polb_only_contigs(record_set):
    polb_only_contigs = []
    for record_id, record in record_set.items():
        atts = get_atts_from_record(record)
        polbs = get_polbs_from_record(record)

        if len(polbs) != 0 and len(atts) == 0:
            polb_only_contigs.append(record_id)

    return polb_only_contigs


def get_assembly_gap_position(unchangeable_record):
    for i_f, feature in enumerate(unchangeable_record.features):
        if feature.type == 'assembly_gap':
            return i_f


def glue_unchangeable_and_att(unchangeable_record, att_record) -> SeqRecord:
    gap_position = get_assembly_gap_position(unchangeable_record)
    if gap_position == 1:
        att_record = cut_att_contig(att_record, 'left')
        glued_record = att_record + unchangeable_record
    else:
        att_record = cut_att_contig(att_record, 'right')
        glued_record = unchangeable_record + att_record
    return glued_record


def check_trna(att_record):
    atts = get_atts_from_record(att_record)
    for att in atts:
        left_feature = att_record.features[att[2] - 1] if att[2] - 1 >= 0 else None
        right_feature = att_record.features[att[2] + 1] if att[2] + 1 < len(att_record.features) else None
        if left_feature is not None and left_feature.type == 'tRNA':
            if 'product' in left_feature.qualifiers:
                if left_feature.qualifiers['product'] == ['tRNA-Leu']:
                    return True
        if right_feature is not None and right_feature.type == 'tRNA':
            if 'product' in right_feature.qualifiers:
                if right_feature.qualifiers['product'] == ['tRNA-Leu']:
                    return True


def get_trna_contig(record_set, att_only_contigs):
    for att_contig in att_only_contigs:
        is_trna = check_trna(record_set[att_contig])
        if is_trna:
            return att_contig
    return None


def cut_att_contig(att_record, direction):
    atts = get_atts_from_record(att_record)
    if direction == 'right':
        new_att_record = att_record[:atts[0][1] + 50]
    else:
        new_att_record = att_record[atts[0][0] - 50:]
    source_feature = SeqFeature(type='source', location=FeatureLocation(1, len(new_att_record), strand=+1),
                                qualifiers={'mol_type': att_record.features[0].qualifiers['mol_type'],
                                            'organism': att_record.features[0].qualifiers['organism'],
                                            'strain': att_record.features[0].qualifiers['strain']})
    new_att_record.features.insert(0, source_feature)

    return new_att_record


def finish_all_separate_contigs(record_set, pipolin_features) -> SeqRecord:
    polb_contigs = get_polb_only_contigs(record_set)
    modify_polb_only_record(polb_contigs, record_set)

    att_contigs = get_att_only_contigs(record_set)
    right_atts = []
    left_atts = []
    trna_contig = get_trna_contig(record_set, att_contigs)
    if trna_contig is None:
        raise NotImplementedError
    else:
        att_contigs.remove(trna_contig)
        right_atts.append(trna_contig)

    if len(att_contigs) == 1:
        print('>>>The only right and left atts are found!')
        left_atts.append(att_contigs[0])
        right_contig = cut_att_contig(record_set[right_atts[0]], 'right')
        left_contig = cut_att_contig(record_set[left_atts[0]], 'left')
        return left_contig + record_set[polb_contigs[0]] + right_contig
    else:
        assert False


def modify_polb_only_record(polb_contigs, record_set):
    if len(polb_contigs) != 1:
        raise NotImplementedError
    polbs = get_polbs_from_record(record_set[polb_contigs[0]])
    if len(polbs) != 1:
        raise NotImplementedError
    if polbs[0][0] > polbs[0][1]:
        raise AssertionError(f'The contig {polb_contigs[0]} probably must be reverted!')
    else:
        print(f'>>> The contig {polb_contigs[0]} probably has the correct direction!')
        assembly_gap = create_assembly_gap_record(record_set[polb_contigs[0]])
        record_set[polb_contigs[0]] = assembly_gap + record_set[polb_contigs[0]] + assembly_gap


def create_single_record(record_set, pipolin_features) -> SeqRecord:
    unchangeable_contigs = get_unchangeable_contigs(record_set)
    print(f'The unchangeable contigs: {unchangeable_contigs}!')
    for contig in unchangeable_contigs:
        record_set[contig] = add_assembly_gap_to_unchangeable(record_set[contig])

    if len(unchangeable_contigs) == 1:
        return finish_one_unchangeable_contig(record_set, unchangeable_contigs)

    elif len(unchangeable_contigs) == 0:
        # TODO: not trivial here...
        return finish_all_separate_contigs(record_set, pipolin_features)
    else:
        raise AssertionError('Only a single pipolin region is expected per genome!')


def finish_one_unchangeable_contig(record_set, unchangeable_contigs) -> SeqRecord:
    att_only_contigs = get_att_only_contigs(record_set)

    if len(att_only_contigs) == 1:
        print('The single record was assembled!!!\n')
        return glue_unchangeable_and_att(record_set[unchangeable_contigs[0]], record_set[att_only_contigs[0]])

    else:
        gap_position = get_assembly_gap_position(record_set[unchangeable_contigs[0]])
        direction = 'left' if gap_position == 1 else 'right'
        if direction == 'right':
            trna_contig = get_trna_contig(record_set, att_only_contigs)
            if trna_contig is None:
                raise NotImplementedError
            else:
                att_only_contigs.remove(trna_contig)
                if len(att_only_contigs) == 1:
                    print('The single record was assembled!!!\n')
                    att_record = cut_att_contig(record_set[att_only_contigs[0]], 'right')
                    return record_set[unchangeable_contigs[0]] + att_record
                else:
                    raise NotImplementedError
        else:
            raise NotImplementedError


def assemble_gapped_pipolins(gb_records: GenBankRecords, pipolin_features):
    for strain_id, record_set in gb_records.items():
        if len(record_set) > 1:
            print(f'Assembling pipolin region for {strain_id}...')
            try:
                gb_records[strain_id][strain_id] = create_single_record(record_set, pipolin_features)
                gb_records[strain_id][strain_id].id = strain_id
            except NotImplementedError:
                print(f'FAILED: {strain_id}')
                continue

            keys_to_delele = []
            for key in gb_records[strain_id].keys():
                if key != strain_id:
                    keys_to_delele.append(key)
            for key in keys_to_delele:
                del gb_records[strain_id][key]


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
    The script takes IN_DIR with *.gbk files (generated after annotation and with included atts),
    detects not assembled pipolins and tries to order the contigs into one sequence,
    filling in gaps with NNs and adding the /assembly_gap feature to the *.gbk files.
    NOTE: WORKS ONLY FOR SHORT PIPOLINS !!!
    """
    gb_records = read_genbank_records(in_dir)
    pipolin_features = get_unique_pipolin_features(gb_records)
    assemble_gapped_pipolins(gb_records, pipolin_features)
    
    check_dir(out_dir)
    write_genbank_records(gb_records, out_dir)
    write_gff_records(gb_records, out_dir)


if __name__ == '__main__':
    main()
