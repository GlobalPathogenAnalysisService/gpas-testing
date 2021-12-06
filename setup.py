#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as f:
    README = f.read()

setup(
    name='gpas-covid-perfect-reads',
    version='1.0',
    description='Create perfect FASTQ files for the SARS-CoV-2 WHO lineages for use in testing',
    author='Philip W Fowler',
    author_email='philip.fowler@ndm.ox.ac.uk',
    url='https://github.com/GenomePathogenAnalysisService/gpas-covid-perfect-reads',
    scripts=['bin/gpas-covid-perfect-reads.py'],
    install_requires=[
        'numpy >= 1.21.2',
        'pyfastaq >= 3.17.0'
        ],
    python_requires='>=3.8',
    license="MIT",
    package_data={'': ['config/MN908947.3.gbk','config/covid-artic-v3.json']},
    include_package_data=True,
    zip_safe=False
    )
