import setuptools

setuptools.setup(
    name='ExplorePipolin',
    version='0.0.a1',
    packages=setuptools.find_packages(),
    install_requires=['click', 'numpy == 1.18.0', 'biopython <= 1.76', 'bcbio-gff', 'prefect >= 0.11.5'],
    # TODO: numpy version is only for MacOS: https://github.com/numpy/numpy/issues/15947
    python_requires='>=3.6',
    package_data={
        'explore_pipolin': ['data/*']
    },
    entry_points={
        'console_scripts': [
            'explore_pipolin=explore_pipolin.explore_pipolin:main',
            'download_taxon_accessions=explore_pipolin.download_taxon_accessions:main',
            'massive_screening=explore_pipolin.massive_screening:main',
            'collect_metadata=explore_pipolin.collect_metadata:main',
            'download_genomes=explore_pipolin.download_genomes:main',
        ]
    },
    test_suite='tests',
    url='https://github.com/liubovch/ExplorePipolin',
    license=''
)
