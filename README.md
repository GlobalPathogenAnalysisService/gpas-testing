# gpas-covid-perfect-reads
Create perfect FASTQ files for the SARS-CoV-2 WHO lineages for use in testing

## Installation

You first need to install [`gumpy`](https://github.com/oxfordmmm/gumpy) as

Once complete, 

```
$ git clone https://github.com/GenomePathogenAnalysisService/gpas-covid-perfect-reads.git
$ cd gpas-covid-perfect-reads
$ python3 setup.py install --user
```

Check it has installed correctly

```
$ gpas-covid-perfect-reads.py --help
usage: gpas-covid-perfect-reads.py [-h] [--reference REFERENCE] --variant_definitions VARIANT_DEFINITIONS
                                   [--variant_name VARIANT_NAME] [--primer_definition PRIMER_DEFINITION]
                                   [--output OUTPUT] [--read_length READ_LENGTH] [--depth DEPTH]

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        the GenBank file of the covid reference
  --variant_definitions VARIANT_DEFINITIONS
                        the path to the variant_definitions repository/folder from phe-genomics
  --variant_name VARIANT_NAME
                        a JSON file specifying the mutations to apply to the covid reference, if none supplied,
                        you'll get a wildtype sequence
  --primer_definition PRIMER_DEFINITION
                        the JSON file specifying the primer scheme used (default is covid-artic-v3.json)
  --output OUTPUT       the stem of the output file
  --read_length READ_LENGTH
                        the read length
  --depth DEPTH         the depth
```

## Usage

You will need to first install `variant_definitions` in an appropriate place in your filesystem

```
$ git clone https://github.com/phe-genomics/variant_definitions
```

First, we can simply create a set of perfect reads for the SARS-CoV-2 reference

```
$ gpas-covid-perfect-reads.py --variant_definitions ../variant_definitions/ --output reference
$ ls -lrt reference*
```
