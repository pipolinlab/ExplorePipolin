import os
from prefect import task
import subprocess


@task
def annotate_pipolins(gquery, pipolins_dir, proteins, root_dir):
    prokka_dir = os.path.join(root_dir, 'prokka')
    os.makedirs(prokka_dir, exist_ok=True)

    subprocess.run(['prokka', '--outdir', prokka_dir,
                    '--prefix', gquery.gquery_id,
                    '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--rfam', '--proteins', proteins, '--force',
                    '--locustag', gquery.gquery_id,
                    os.path.join(pipolins_dir, gquery.gquery_id + '.fa')])

    return prokka_dir
