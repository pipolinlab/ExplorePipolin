from Bio.SeqIO import SeqRecord

from explore_pipolin.common import Pipolin, FeatureType, Range, AttFeature, AttType, get_rec_id_by_contig_id
from explore_pipolin.utilities.io import SeqIORecords
import explore_pipolin.settings as settings


_PRODUCTS_TO_COLOUR = {'default': '255 250 240'}


def read_colors(colors_tsv: str) -> None:
    with open(colors_tsv) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            values = line.strip().split(sep='\t')
            if len(values) == 0:   # skip empty lines
                continue
            elif len(values) != 3:
                raise AssertionError(f'{len(values)} columns in the line {line.strip()}.\n'
                                     f'3 columns are expected.')
            else:
                _PRODUCTS_TO_COLOUR[values[0]] = values[2]


def easyfig_add_colours(gb_records: SeqIORecords, pipolin: Pipolin):
    colours = settings.get_instance().colours
    read_colors(colours)

    for record in gb_records.values():
        add_colours(record)

    for fragment in pipolin.fragments:
        fragment_shift = fragment.start

        for att in [f for f in fragment.features if f.ftype == FeatureType.ATT]:
            att_start, att_end = (att.start - fragment_shift), (att.end - fragment_shift)
            att: AttFeature
            for feature in gb_records[get_rec_id_by_contig_id(gb_records, fragment.contig_id)].features:
                if feature.location.start == att_start and feature.location.end == att_end:
                    if att.att_type == AttType.CONSERVED:
                        if 'FeatureType.ATT.CONSERVED' in _PRODUCTS_TO_COLOUR:
                            feature.qualifiers['colour'] = _PRODUCTS_TO_COLOUR['FeatureType.ATT.CONSERVED']
                    else:
                        if 'FeatureType.ATT.DENOVO' in _PRODUCTS_TO_COLOUR:
                            feature.qualifiers['colour'] = _PRODUCTS_TO_COLOUR['FeatureType.ATT.DENOVO']

        for ttrna in [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]:
            ttrna_start, ttrna_end = (ttrna.start - fragment_shift), (ttrna.end - fragment_shift)

            for feature in gb_records[get_rec_id_by_contig_id(gb_records, fragment.contig_id)].features:
                feature_range = Range(start=feature.location.start, end=feature.location.end)
                if feature.type == 'tRNA' and feature_range.is_overlapping(Range(start=ttrna_start, end=ttrna_end)):
                    if 'FeatureType.TARGET_TRNA' in _PRODUCTS_TO_COLOUR:
                        feature.qualifiers['colour'] = _PRODUCTS_TO_COLOUR['FeatureType.TARGET_TRNA']


def add_colours(record: SeqRecord):
    for feature in record.features:
        _colour_feature(feature.qualifiers)


def _colour_feature(qualifiers):
    if 'product' in qualifiers:
        for product in qualifiers['product']:
            if product in _PRODUCTS_TO_COLOUR:
                qualifiers['colour'] = [_PRODUCTS_TO_COLOUR[product]]
            else:
                qualifiers['colour'] = [_PRODUCTS_TO_COLOUR['default']]
    elif 'linkage_evidence' in qualifiers:
        is_paired_ends = qualifiers['linkage_evidence'] == ['paired-ends']
        is_pipolin_structure = qualifiers['linkage_evidence'] == ['pipolin_structure']

        if is_paired_ends and ('paired-ends' in _PRODUCTS_TO_COLOUR):
            qualifiers['colour'] = [_PRODUCTS_TO_COLOUR['paired-ends']]

        elif is_pipolin_structure and ('pipolin_structure' in _PRODUCTS_TO_COLOUR):
            qualifiers['colour'] = [_PRODUCTS_TO_COLOUR['pipolin_structure']]

        else:
            qualifiers['colour'] = [_PRODUCTS_TO_COLOUR['default']]
    else:
        qualifiers['colour'] = [_PRODUCTS_TO_COLOUR['default']]
