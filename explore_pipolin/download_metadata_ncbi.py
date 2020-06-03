import click
import subprocess
from collections import namedtuple
from Bio import SearchIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def get_assembly_info(acc_id):
    acc_id_info = subprocess.run(['esearch', '-db', 'nuccore', '-query', f'{acc_id}'], stdout=subprocess.PIPE)
    assembly_info = subprocess.run(['elink', '-target', 'assembly', '-name', 'nuccore_assembly'],
                                   input=acc_id_info.stdout, stdout=subprocess.PIPE)
    return assembly_info.stdout


def get_assembly_acc_and_species_name(assembly_info):
    assembly_info_fetched = subprocess.run(['efetch', '-format', 'docsum'],
                                           input=assembly_info, stdout=subprocess.PIPE)
    assembly_info_extracted = subprocess.run(
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
    seqs_info = subprocess.run(['elink', '-target', 'nucleotide', '-name', 'assembly_nuccore_insdc'],
                                input=assembly_info, stdout=subprocess.PIPE)
    seqs_info_fetched = subprocess.run(['efetch', '-format', 'docsum'],
                                       input=seqs_info.stdout, stdout=subprocess.PIPE)
    seqs_info_filtered = subprocess.run(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'AccessionVersion', 'Title', 'Strain'],
        input=seqs_info_fetched.stdout, stdout=subprocess.PIPE)
    plasmids_excluded = subprocess.run(['grep', '-v', 'plasmid'],
                                       input=seqs_info_filtered.stdout, stdout=subprocess.PIPE)
    return plasmids_excluded.stdout.decode(encoding='UTF8')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('blast-tab', type=click.Path(exists=True))
@click.argument('out_file', type=click.Path())
def download_metadata_ncbi(blast_tab, out_file):
    """
    BLAST_TAB is a tabular file (hittable or outfmt=7) that was generated when searching E. coli genomes
    against "reference" pi-polB from E. coli 3-373-03_S1_C2 (NZ_JNMI01000006.1). E. coli genome accessions
    will be extracted from the file and the required metadata will be fetch from NCBI into the OUT_FILE.
    NOTE: requires "ncbi-entrez-direct" package!
    """
    blast_result = SearchIO.read(blast_tab, 'blast-tab', comments=True)
    blast_ids = [i.id for i in blast_result]   # acc_ids are all unique here!
    blast_ids.append(blast_result.id)   # query id == NZ_JNMI01000006.1

    Metadata = namedtuple('Metadata', ['assembly_id', 'species_name', 'strain_name', 'acc_ids'])
    metadata = []

    for num, blast_id in enumerate(blast_ids):
        print(f'{num + 1} Fetching metadata for {blast_id}...')
        assembly_info = get_assembly_info(blast_id)
        assembly_acc, species_name = get_assembly_acc_and_species_name(assembly_info)
        seqs_info = get_and_filter_seqs_info(assembly_info)

        strings = seqs_info.strip().split(sep='\n')
        acc_ids = []
        acc_id, _, strain = strings[0].strip().split(sep='\t')
        acc_ids.append(acc_id)
        acc_ids.extend(line.strip().split(sep='\t')[0] for line in strings[1:])

        metadata.append(Metadata(assembly_acc, species_name, strain, ','.join(acc_ids)))

    # add the info about LREC strains
    metadata.append(Metadata('chr_LREC237', 'Escherichia coli', 'LREC237', '-'))
    for i in range(239, 263):
        metadata.append(Metadata(f'chr_LREC{i}', 'Escherichia coli', f'LREC{i}', '-'))

    with open(out_file, 'w') as ouf:
        for line in metadata:
            print('\t'.join(list(line)), file=ouf)


if __name__ == '__main__':
    download_metadata_ncbi()
