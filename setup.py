import setuptools

setuptools.setup(
    name='explore-pipolin',
    version='0.0.1',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': ['pipolin_finder:main']
    }
)