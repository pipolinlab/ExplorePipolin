import setuptools

setuptools.setup(
    name='ExplorePipolin',
    version='0.0.1',
    packages=setuptools.find_packages(),
    install_requires=['click', 'biopython', 'bcbio-gff', 'prefect'],
    python_requires='>=3.6',
    package_data={
        'explore_pipolin': ['data/attL.fa', 'data/HHpred_proteins.faa', 'data/pi-polB.faa']
    },
    entry_points={
        'console_scripts': ['explore_pipolin=explore_pipolin.main:explore_pipolin']
    }
)
