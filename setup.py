import setuptools

setuptools.setup(
    name='PipolinFinder',
    version='0.0.1',
    packages=setuptools.find_packages(),
    install_requires=['click', 'Biopython', 'bcbio-gff', 'prefect'],
    python_requires='~=3.7',
    package_data={
        'pipolin_finder': ['data/attL.fa', 'data/HHpred_proteins.faa', 'data/pi-polB.faa']
    },
    entry_points={
        'console_scripts': ['pipolin_finder=pipolin_finder.main:explore_pipolins']
    }
)
