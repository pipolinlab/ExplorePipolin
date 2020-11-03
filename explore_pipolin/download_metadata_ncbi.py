import os
from typing import Tuple, Sequence

import click
import subprocess


from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.utilities.external_tools import subprocess_with_retries, EmptyResult, NoAssembly


def read_gb_ids(input_file):
    with open(input_file) as inf:
        return [line.strip() for line in inf if line != '']


def get_assembly_info(gb_id, id_type) -> Tuple[str, bool, str]:
    if id_type == 'nucleotide':
        id_info = subprocess_with_retries(['esearch', '-db', 'nuccore', '-query', f'{gb_id}'],
                                          text=True, input='', stdout=subprocess.PIPE)
    elif id_type == 'protein':
        id_info_protein = subprocess_with_retries(['esearch', '-db', 'protein', '-query', f'{gb_id}'],
                                                  text=True, input='', stdout=subprocess.PIPE)
        id_info = subprocess_with_retries(['elink', '-target', 'nuccore'], text=True,
                                          input=id_info_protein.stdout, stdout=subprocess.PIPE)
    else:
        raise AssertionError('Wrong --id-type !!!')

    id_info_docsum = subprocess_with_retries(['efetch', '-format', 'docsum'],
                                             text=True, input=id_info.stdout, stdout=subprocess.PIPE)
    id_title = subprocess_with_retries(['xtract', '-pattern', 'DocumentSummary', '-element', 'Title'],
                                       text=True, input=id_info_docsum.stdout, stdout=subprocess.PIPE)

    is_plasmid = ' plasmid ' in id_title.stdout

    assembly_info = subprocess_with_retries(['elink', '-target', 'assembly', '-name', 'nuccore_assembly'],
                                            text=True, input=id_info.stdout, stdout=subprocess.PIPE)
    return assembly_info.stdout, is_plasmid, id_title.stdout


def get_assembly_acc_and_species_name(assembly_info: str) -> Tuple[Sequence, Sequence]:
    assembly_info_fetched = subprocess_with_retries(['efetch', '-format', 'docsum'], text=True,
                                                    input=assembly_info, stdout=subprocess.PIPE)
    assembly_info_extracted = subprocess_with_retries(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'AssemblyAccession', 'SpeciesName'],
        text=True, input=assembly_info_fetched.stdout, stdout=subprocess.PIPE)

    assemblies = [line.strip().split('\t')[0] for line in assembly_info_extracted.stdout.strip().split('\n')]
    species_names = [line.strip().split('\t')[1] for line in assembly_info_extracted.stdout.strip().split('\n')]

    unique_ids = set([a.split('.')[0] for a in assemblies])
    # All assemblies are different
    if len(unique_ids) == len(assemblies):
        return assemblies, species_names
    # Not all assemblies are different
    else:
        # There might be old release assemblies as well, but we need the latest!
        if len(unique_ids) == 1:
            index_to_leave = assemblies.index(max(assemblies, key=lambda x: x[-1]))
            return [assemblies[index_to_leave]], [species_names[index_to_leave]]
        else:
            raise NotImplementedError


def get_and_filter_assembly_info(assembly_acc, is_plasmid, title) -> str:
    assembly = subprocess_with_retries(['esearch', '-db', 'assembly', '-query', assembly_acc],
                                       text=True, input='', stdout=subprocess.PIPE)
    assembly_info = subprocess_with_retries(['elink', '-target', 'nucleotide', '-name', 'assembly_nuccore_insdc'],
                                            text=True, input=assembly.stdout, stdout=subprocess.PIPE)
    assembly_info_fetched = subprocess_with_retries(['efetch', '-format', 'docsum'], text=True,
                                                    input=assembly_info.stdout, stdout=subprocess.PIPE)
    seqs_info_filtered = subprocess_with_retries(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'AccessionVersion', 'Title', 'Strain'],
        text=True, input=assembly_info_fetched.stdout, stdout=subprocess.PIPE)
    if is_plasmid:
        results = subprocess.run(['grep', f'{title.strip()}'], text=True, input=seqs_info_filtered.stdout,
                                 stdout=subprocess.PIPE)
    else:
        results = subprocess.run(['grep', '-v', 'plasmid'], text=True, input=seqs_info_filtered.stdout,
                                 stdout=subprocess.PIPE)
    return results.stdout


def get_strain_and_acc_ids(seqs_info):
    strings = seqs_info.strip().split(sep='\n')
    acc_ids = []
    try:
        acc_id, _, strain = strings[0].strip().split(sep='\t')
    except ValueError:
        try:
            acc_id, _ = strings[0].strip().split(sep='\t')
            strain = ''
        except ValueError:
            if strings == ['']:
                return '', ''
            else:
                raise NotImplementedError
    acc_ids.append(acc_id)
    acc_ids.extend(line.strip().split(sep='\t')[0] for line in strings[1:])
    return strain, acc_ids


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-file', type=click.Path(exists=True))
@click.argument('out-file', type=click.Path())
@click.option('--id-type', default='nucleotide', show_default=True,
              type=click.Choice(['nucleotide', 'protein']),
              help='type of GenBank ids in the list')
@click.option('--ncbi-api-key', default='')
def download_metadata_ncbi(in_file, out_file, id_type, ncbi_api_key):
    """
    IN_FILE -- contains a list of GenBank ids each on its own line.
    NOTE: requires "ncbi-entrez-direct" package!
    """
    os.environ['NCBI_API_KEY'] = ncbi_api_key

    gb_ids = read_gb_ids(in_file)

    column_names = ['GenBankID', 'AssemblyAccession', 'SpeciesName', 'Strain', 'IsPlasmid', 'GenBank']

    with open(out_file, 'w') as ouf:
        print('\t'.join(column_names), file=ouf)
        for num, gb_id in enumerate(gb_ids):
            print(f'{num + 1} Fetching metadata for {gb_id}...')
            try:
                assembly_info, is_plasmid, title = get_assembly_info(gb_id, id_type)
                assemblies, species_names = get_assembly_acc_and_species_name(assembly_info)

                for a, s in zip(assemblies, species_names):
                    try:
                        seqs_info = get_and_filter_assembly_info(a, is_plasmid, title)
                    except NoAssembly:
                        message = f'No assembly sequence data for the query {a} !!!'
                        print('\t'.join([gb_id, message]), file=ouf)
                        print(message)
                        continue

                    strain, acc_ids = get_strain_and_acc_ids(seqs_info)
                    print('\t'.join([gb_id, a, s, strain, str(is_plasmid), ','.join(acc_ids)]), file=ouf)
            except EmptyResult:
                message = f'NO RESULTS for the query {gb_id} !!!'
                print('\t'.join([gb_id, message]), file=ouf)
                print(message)
                continue

        # # add the info about LREC strains
        # print('\t'.join('chr_LREC237', 'Escherichia coli', 'LREC237', '-'), file=ouf)
        # for i in range(239, 263):
        #     print('\t'.join(f'chr_LREC{i}', 'Escherichia coli', f'LREC{i}', '-'))


if __name__ == '__main__':
    download_metadata_ncbi()
