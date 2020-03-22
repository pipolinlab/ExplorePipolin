#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import tempfile
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from utilities import CONTEXT_SETTINGS
from utilities import check_dir
from io import StringIO


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('blast-tab', type=click.Path(exists=True))
@click.argument('out-dir')
def download_metadata_ncbi(blast_tab, out_dir):
    """
    BLAST_TAB is a tabular file (hittable or outfmt=7) that was generated when searching E. coli genomes
    against "reference" pi-polB from E. coli 3-373-03_S1_C2 (NZ_JNMI01000006.1). E. coli genome accessions
    will be extracted from the file and the genome sequences in FASTA format will be downloaded from NCBI
    to the OUT_DIR.
    NOTE: requires "ncbi-entrez-direct" package!
    """
    blast_result = SearchIO.read(blast_tab, 'blast-tab', comments=True)
    acc_ids = [i.id for i in blast_result]   # acc_ids are all unique here!
    acc_ids.append(blast_result.id)   # query id == NZ_JNMI01000006.1

    with tempfile.TemporaryDirectory() as tmpdir:
        assembly_acc_ids = []
        species_names = []
        strain_names = []
        for acc_id in acc_ids:
            acc_id_info = subprocess.Popen(['esearch', '-db', 'nuccore', '-query', f'{acc_id}'],
                                       stdout=subprocess.PIPE)
            assembly_info = subprocess.Popen(['elink', '-target', 'assembly', '-name', 'nuccore_assembly'],
                                      stdin=acc_id_info.stdout, stdout=subprocess.PIPE)
            essembly_seqs_info = subprocess.Popen(['elink', '-target', 'nucleotide', '-name', 'assembly_nuccore_refseq'],
                                      stdin=assembly_info.stdout, stdout=subprocess.PIPE)
            assembly_info_fetched = subprocess.Popen(['efetch', '-format', 'docsum'],
                                                     stdin=assembly_info.stdout, stdout=subprocess.PIPE)
            assembly_info_extracted = subprocess.run(
                ['xtract', '-pattern', 'DocumentSummary', '-element', 'AssemblyAccession', 'SpeciesName'],
                                                       stdin=assembly_info_fetched.stdout, stdout=subprocess.PIPE)
            assembly_acc, species_name = assembly_info_extracted.stdout.decode(encoding='UTF8').split(sep='\t')
            assembly_acc_ids.append(assembly_acc)
            species_names.append(species_name)

            assembly_seqs_info_fetched = subprocess.Popen(['efetch', '-format', 'docsum'],
                                                          stdin=essembly_seqs_info.stdout, stdout=subprocess.PIPE)
            info_extracted = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'AccessionVersion',
                                       'Title', 'Strain'],
                                      stdin=assembly_seqs_info_fetched.stdout, stdout=subprocess.PIPE)
            acc_id_info.stdout.close()
            assembly_info.stdout.close()
            essembly_seqs_info.stdout.close()
            assembly_info_fetched.stdout.close()
            assembly_seqs_info_fetched.stdout.close()

            with open(os.path.join(tmpdir, f'{acc_id}' + '.tab'), 'w') as ouf:
                subprocess.run(['grep', '-v', 'plasmid'], stdin=info_extracted.stdout, stdout=ouf)
                info_extracted.stdout.close()

        check_dir(out_dir)
        for summary in os.listdir(tmpdir):
            with open(os.path.join(tmpdir, summary)) as inf:
                seqs = []
                for line in inf:
                    acc_id, _, strain = line.strip().split(sep='\t')
                    args = ['efetch', '-format', 'fasta', '-id', f'{acc_id}', '-db', 'sequences']
                    seqs.append(SeqIO.read(StringIO(subprocess.run(args, stdout=subprocess.PIPE).stdout.decode(encoding='UTF8')),
                                           format='fasta'))
                    strain_ = '_'.join(strain.split())
                strain_names.append(strain_)
                with open(os.path.join(out_dir, strain_ + '.fa'), 'w') as ouf:
                    SeqIO.write(seqs, ouf, format='fasta')


if __name__ == '__main__':
    download_metadata_ncbi()
