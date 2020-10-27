from enum import Enum, auto

import click
import subprocess
from Bio import SearchIO
from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.utilities.external_tools import subprocess_with_retries

# TODO: write tests!


def get_gb_ids(acc_file, acc_type):
    if acc_type is AccType.LIST:
        with open(acc_file) as inf:
            gb_ids = [line.strip() for line in inf]
    else:
        blast_result = SearchIO.read(acc_file, 'blast-tab', comments=True)
        gb_ids = [i.id for i in blast_result]  # acc_ids are all unique here!
        gb_ids.append(blast_result.id)  # add also the query id
    return gb_ids


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


class AccType(Enum):
    LIST = auto()
    HITTABLE = auto()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('acc-file', type=click.Path(exists=True))
@click.argument('out-file', type=click.Path())
@click.option('--acc-type', default=AccType.LIST.name.lower(), show_default=True,
              type=click.Choice([acc.name.lower() for acc in AccType]))
def download_metadata_ncbi(acc_file: str, out_file: str, acc_type):
    """
    TODO: rewrite
    BLAST_TAB is a tabular file (hittable or outfmt=7) that was generated when searching E. coli genomes_ecoli
    against "reference" pi-polB from E. coli 3-373-03_S1_C2 (NZ_JNMI01000006.1). E. coli genome accessions
    will be extracted from the file and the required metadata will be fetch from NCBI into the OUT_FILE.
    NOTE: requires "ncbi-entrez-direct" package!
    """
    acc_type = AccType[acc_type.upper()]
    gb_ids = get_gb_ids(acc_file, acc_type)

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
    # metadata.append(Metadata('chr_LREC237', 'Escherichia coli', 'LREC237', '-'))
    # for i in range(239, 263):
    #     metadata.append(Metadata(f'chr_LREC{i}', 'Escherichia coli', f'LREC{i}', '-'))


if __name__ == '__main__':
    download_metadata_ncbi()
