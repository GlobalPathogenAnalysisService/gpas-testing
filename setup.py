#!/usr/bin/env python

from setuptools import setup
from gpas_covid_synthetic_reads import __version__

with open("README.md", "r") as f:
    README = f.read()

setup(
    name='gpas-covid-synthetic-reads',
    version=__version__,
    description='Create perfect FASTQ files for the SARS-CoV-2 WHO lineages for use in testing',
    author='Philip W Fowler',
    author_email='philip.fowler@ndm.ox.ac.uk',
    url='https://github.com/GenomePathogenAnalysisService/gpas-covid-synthetic-reads',
    scripts=['bin/gpas-covid-synthetic-reads.py'],
    install_requires=[
        'numpy >= 1.21.2',
        'pyfastaq >= 3.17.0'
        ],
    python_requires='>=3.8',
    license="MIT",
    package_data={'': ['data/*']},
    include_package_data=True,
    zip_safe=False
    )
