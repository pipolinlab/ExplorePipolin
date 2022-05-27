import setuptools

setuptools.setup(
    name='ExplorePipolin',
    version='0.0.1',
    packages=setuptools.find_packages(),
    install_requires=['click', 'biopython <= 1.76', 'bcbio-gff', 'prefect >= 0.11.5'],
    python_requires='>=3.7',
    package_data={
        'explore_pipolin': ['data/*']
    },
    entry_points={
        'console_scripts': [
            'explore_pipolin=explore_pipolin.main:main',
        ]
    },
    test_suite='tests',
    url='https://github.com/liubovch/ExplorePipolin',
    license='GPLv3'
)
