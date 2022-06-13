#!/usr/bin/env python

from setuptools import setup
from gpas_covid_synthetic_reads import __version__

with open("README.md", "r") as f:
    README = f.read()

setup(
    name='gpas_covid_synthetic_reads',
    version=__version__,
    description='Create synthetic FASTQ files of different SARS-CoV-2 lineages for testing',
    author='Philip W Fowler',
    author_email='philip.fowler@ndm.ox.ac.uk',
    url='https://github.com/GenomePathogenAnalysisService/gpas_covid_synthetic_reads',
    scripts=['bin/gpas-covid-synreads-create.py',
             'bin/gpas-build-uploadcsv.py',
             'bin/gpas-covid-synreads-analyse.py',
             'bin/gpas-tb-synreads-create.py',
             'bin/gpas-test-fn4.py'],
    install_requires=[
        'pyfastaq >= 3.17.0',
        'pyyaml',
        'numpy',
        'tqdm',\
        'pandas',\
        'requests'
        ],
    python_requires='>=3.8',
    license="MIT",
    packages=['gpas_covid_synthetic_reads'],
    package_data={'': ['data/*']},
    include_package_data=True,
    zip_safe=False
    )
