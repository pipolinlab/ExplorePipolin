import setuptools

setuptools.setup(
    name='ExplorePipolin',
    version='0.0.a1',
    packages=setuptools.find_packages(),
    install_requires=['click', 'biopython <= 1.76', 'bcbio-gff', 'prefect >= 0.11.5'],
    python_requires='>=3.6',
    package_data={
        'explore_pipolin': ['data/attL.fa', 'data/HHpred_proteins.faa', 'data/pi-polB.faa']
    },
    entry_points={
        'console_scripts': ['explore_pipolin=explore_pipolin.main:explore_pipolin',
                            'download_metadata_ncbi=explore_pipolin.download_metadata_ncbi:download_metadata_ncbi']
    },
    test_suite='tests',
    url='https://github.com/liubovch/ExplorePipolin',
    license=''
)
