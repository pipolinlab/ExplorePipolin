import click
import subprocess
from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.utilities.external_tools import subprocess_with_retries


def read_gb_ids(input_file):
    with open(input_file) as inf:
        return [line.strip() for line in inf if line != '']


def get_assembly_info(acc_id):
    acc_id_info = subprocess_with_retries(['esearch', '-db', 'nuccore', '-query', f'{acc_id}'],
                                          input='', stdout=subprocess.PIPE)

    assembly_info = subprocess_with_retries(['elink', '-target', 'assembly', '-name', 'nuccore_assembly'],
                                            input=acc_id_info.stdout, stdout=subprocess.PIPE)
    return assembly_info.stdout


def get_assembly_acc_and_species_name(assembly_info):
    assembly_info_fetched = subprocess_with_retries(['efetch', '-format', 'docsum'],
                                                    input=assembly_info, stdout=subprocess.PIPE)
    assembly_info_extracted = subprocess_with_retries(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'AssemblyAccession', 'SpeciesName'],
        input=assembly_info_fetched.stdout, stdout=subprocess.PIPE)
    info_strings = assembly_info_extracted.stdout.decode(encoding='UTF8').strip().split(sep='\n')
    # There might be old release assemblies as well, but we need the latest!
    if len(info_strings) == 1:
        return info_strings[0].strip().split(sep='\t')
    else:
        check_ids = []
        for line in info_strings:
            check_ids.append(line.strip().split('\t'))
        return max(check_ids, key=lambda x: x[0])


def get_and_filter_seqs_info(assembly_info):
    seqs_info = subprocess_with_retries(['elink', '-target', 'nucleotide', '-name', 'assembly_nuccore_insdc'],
                                        input=assembly_info, stdout=subprocess.PIPE)
    seqs_info_fetched = subprocess_with_retries(['efetch', '-format', 'docsum'],
                                                input=seqs_info.stdout, stdout=subprocess.PIPE)
    seqs_info_filtered = subprocess_with_retries(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'AccessionVersion', 'Title', 'Strain'],
        input=seqs_info_fetched.stdout, stdout=subprocess.PIPE)
    plasmids_excluded = subprocess.run(['grep', '-v', 'plasmid'],
                                       input=seqs_info_filtered.stdout, stdout=subprocess.PIPE)
    return plasmids_excluded.stdout.decode(encoding='UTF8')


def get_strain_and_acc_ids(seqs_info):
    strings = seqs_info.strip().split(sep='\n')
    acc_ids = []
    acc_id, _, strain = strings[0].strip().split(sep='\t')
    acc_ids.append(acc_id)
    acc_ids.extend(line.strip().split(sep='\t')[0] for line in strings[1:])
    return strain, acc_ids


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-file', type=click.Path(exists=True))
@click.argument('out-file', type=click.Path())
def download_metadata_ncbi(in_file, out_file):
    """
    IN_FILE -- contains a list of GenBank ids each on its own line
    NOTE: requires "ncbi-entrez-direct" package!
    """
    gb_ids = read_gb_ids(in_file)

    column_names = ['AssemblyAccession', 'SpeciesName', 'Strain', 'GenBank']

    with open(out_file, 'w') as ouf:
        print('\t'.join(column_names), file=ouf)
        for num, blast_id in enumerate(gb_ids):
            print(f'{num + 1} Fetching metadata for {blast_id}...')
            assembly_info = get_assembly_info(blast_id)
            assembly_acc, species_name = get_assembly_acc_and_species_name(assembly_info)
            seqs_info = get_and_filter_seqs_info(assembly_info)
            strain, acc_ids = get_strain_and_acc_ids(seqs_info)

            print('\t'.join([assembly_acc, species_name, strain, ','.join(acc_ids)]), file=ouf)

        # # add the info about LREC strains
        # print('\t'.join('chr_LREC237', 'Escherichia coli', 'LREC237', '-'), file=ouf)
        # for i in range(239, 263):
        #     print('\t'.join(f'chr_LREC{i}', 'Escherichia coli', f'LREC{i}', '-'))


if __name__ == '__main__':
    download_metadata_ncbi()
