import csv
import os
import shelve
from typing import Sequence, MutableSequence, Mapping
from itertools import groupby
import subprocess
from Bio import SearchIO, SeqIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class Feature:
    def __init__(self, start, end, node):
        self.start: int = start
        self.end: int = end
        self.node: str = node

    def normalize(self):
        # TODO: doesn't work as there is no info about direction! Do not use 'blast-tab'!
        if self.start > self.end:
            return Feature(self.end, self.start, self.node)
        else:
            return self


class Pipolin:
    def __init__(self, strain_id):
        self.strain_id: str = strain_id
        self.polymerases: MutableSequence[Feature] = []
        self.atts: MutableSequence[Feature] = []

    def get_pipolin_bounds(self, long):
        if not self.is_complete_genome():
            raise AssertionError('Should be complete!')

        polymerases = sorted((i.normalize() for i in self.polymerases), key=lambda p: p.start)
        atts = sorted((i.normalize() for i in self.atts), key=lambda p: p.start)

        if not self._is_polymerase_inside(atts, polymerases):
            raise AssertionError('The polymerases are not within att bounds!')

        if len(atts) > 3:
            raise AssertionError('There are more than 3 atts!')

        if long:
            if len(atts) == 3:
                return atts[0].start - 50, atts[2].end + 50
            else:
                return atts[0].start - 50, atts[1].end + 50
        else:
            if polymerases[-1].end < atts[1].end:
                return atts[0].start - 50, atts[1].end + 50
            else:
                return atts[1].start - 50, atts[2].end + 50

    @staticmethod
    def _is_polymerase_inside(atts, polymerases):
        return atts[0].start < polymerases[0].start and polymerases[-1].end < atts[-1].end

    @staticmethod
    def _dict_by_node_normalized(features):
        return {node: sorted(list(ps), key=lambda p: p.start) for node, ps
                in groupby((i.normalize() for i in features), key=lambda x: x.node)}

    @classmethod
    def _add_dangling_atts(cls, atts, things_to_return, length_by_contig):
        for node, features in atts.items():
            if node not in things_to_return:
                things_to_return[node] = cls._get_dangling_feature(features, node, length_by_contig)

    @classmethod
    def _get_dangling_feature(cls, features, node, length_by_contig):
        pos_left = features[0].start - 50000
        left = pos_left if pos_left > 0 else 0
        pos_right = features[-1].end + 50000
        right = pos_right if pos_right < length_by_contig[node] else -1
        return left, right

    @classmethod
    def _get_pipolin_two_atts(cls, atts, polymerases, long):
        if cls._is_polymerase_inside(atts, polymerases):
            if long:
                return atts[0].start - 50, atts[2].end + 50
            else:
                return atts[0].start - 50, atts[1].end + 50
        else:
            raise AssertionError('The polymerases are not within att bounds!')

    @classmethod
    def _get_pipolin_single_att(cls, att: Feature, polymerases, length_by_contig):
        if att.end < polymerases[0].start:
            left = att.start - 50
            pos_right = polymerases[-1].end + 50000
            right = pos_right if pos_right < length_by_contig[att.node] else -1
        else:
            pos_left = polymerases[0].start - 50000
            left = pos_left if pos_left > 0 else 0
            right = att.end + 50
        return left, right

    def get_contigs_with_bounds(self, length_by_contig, long):
        # TODO: check boarders when +/-50 nt !
        # TODO: check indexes carefully: 0-based and 1-based !!!
        polymerases = self._dict_by_node_normalized(self.polymerases)
        atts = self._dict_by_node_normalized(self.atts)

        things_to_return = {}
        for node, features in polymerases.items():
            self.add_pipolin_for_node(atts, features, node, things_to_return, length_by_contig, long)

        return things_to_return

    def add_pipolin_for_node(self, atts, features, node, things_to_return, length_by_contig, long):
        if node in atts:
            if long:
                if len(atts[node]) == 3:
                    things_to_return[node] = self._get_pipolin_two_atts(atts[node], features, long)
            if len(atts[node]) == 2:
                things_to_return[node] = self._get_pipolin_two_atts(atts[node], features, long=False)
            elif len(atts[node]) == 1:
                things_to_return[node] = self._get_pipolin_single_att(atts[node][0], features, length_by_contig)
                self._add_dangling_atts(atts, things_to_return, length_by_contig)
            else:
                raise AssertionError('More than two atts on one contig!')
        else:
            things_to_return[node] = self._get_dangling_feature(features, node, length_by_contig)
            self._add_dangling_atts(atts, things_to_return, length_by_contig)

    def is_complete_genome(self):
        return self.strain_id == self.polymerases[0].node


def check_dir(out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


def ncbi_acc_download(acc_ids):
    for acc_id in acc_ids:
        print(f'Downloading {acc_id}...')
        # TODO: fix local path!
        subprocess.run(['/home/liubov/.local/bin/ncbi-acc-download', '-F', 'fasta', acc_id])


# def read_blasttab(blast_tab):
#     return SearchIO.read(blast_tab, 'blast-tab', comments=True)


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')


def save_to_shelve(out_file, object_to_save, key_name):
    shelve_db = shelve.open(os.path.splitext(out_file)[0])
    shelve_db[key_name] = object_to_save
    shelve_db.close()


def read_from_shelve(shelve_file, key_name):
    shelve_db = shelve.open(os.path.splitext(shelve_file)[0])
    shelve_object = shelve_db[key_name]
    shelve_db.close()
    return shelve_object


# def save_as_csv(Pipolin, pipolins, out_file):
#     # TODO: rewrite this!
#     with open(out_file, 'w') as ouf:
#         print(','.join(Pipolin._fields), file=ouf)
#         for i_p, _ in enumerate(pipolins):
#             words = ['None' if i is None else str(i) for i in pipolins[i_p]]
#             print(','.join(words), file=ouf)


def get_roary_groups(roary_dir) -> Mapping[str, Mapping[str, Sequence[str]]]:
    roary_groups = {}
    with open(os.path.join(roary_dir, 'gene_presence_absence.csv')) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        header = next(reader)
        for entry in reader:
            group_name = entry[0]
            genes = {genome: genes.split('\t') if genes else [] for genome, genes in zip(header[14:],entry[14:])}
            roary_groups[group_name] = genes
    return roary_groups


def blast_seqs_against_seq(dir_with_seqs, seq, output_dir):
    check_dir(output_dir)
    genomes = os.listdir(dir_with_seqs)

    for genome in genomes:
        with open(os.path.join(output_dir, f'{genome[:-3]}_fmt5.txt'), 'w') as ouf:
            subprocess.run(['blastn', '-query', seq,
                            '-subject', f'{os.path.join(dir_with_seqs, genome)}',
                            '-outfmt', '5'], stdout=ouf)


GenBankRecords = Mapping[str, Mapping[str, SeqIO.SeqRecord]]


def read_genbank_records(annot_dir) -> GenBankRecords:
    genbank_records: GenBankRecords = {}
    for file in os.listdir(annot_dir):
        if file.endswith('.gbk'):
            record = SeqIO.to_dict(SeqIO.parse(os.path.join(annot_dir, file), format='genbank'))
            genbank_records[os.path.splitext(file)[0]] = record
    return genbank_records


def write_genbank_records(genbank_records, new_annot_dir):
    for key, value in genbank_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(new_annot_dir, f'{key}.gbk'), 'w') as ouf:
            SeqIO.write(records, ouf, 'genbank')