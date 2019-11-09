import subprocess

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def ncbi_acc_download(ids):
    for genome_id in ids:
        # TODO fix local path!
        subprocess.run(['/home/liubov/.local/bin/ncbi-acc-download', '-F', 'fasta', genome_id])


