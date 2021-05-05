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
            'explore_pipolin=explore_pipolin.main:explore_pipolin',
            'massive_screening_ena_xml=explore_pipolin.massive_screening_ena_xml:massive_screening',
            'extract_metadata_ena_xml=explore_pipolin.extract_metadata_ena_xml:extract_metadata_all',
            'download_pipolin_genomes=explore_pipolin.download_pipolin_genomes:download_pipolin_genomes'
        ]
    },
    test_suite='tests',
    url='https://github.com/liubovch/ExplorePipolin',
    license=''
)
