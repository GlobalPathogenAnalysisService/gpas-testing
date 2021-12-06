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
$ usage: gpas-covid-perfect-reads.py [-h] --variant_definitions VARIANT_DEFINITIONS --output OUTPUT
                                   [--variant_name VARIANT_NAME] [--reference REFERENCE]
                                   [--primer_definition PRIMER_DEFINITION] [--read_length READ_LENGTH]
                                   [--depth DEPTH]

optional arguments:
  -h, --help            show this help message and exit
  --variant_definitions VARIANT_DEFINITIONS
                        the path to the variant_definitions repository/folder from phe-genomics
  --output OUTPUT       the stem of the output file
  --variant_name VARIANT_NAME
                        a JSON file specifying the mutations to apply to the covid reference (if not specified,
                        you'll get a wildtype sequence)
  --reference REFERENCE
                        the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will
                        be used)
  --primer_definition PRIMER_DEFINITION
                        the JSON file specifying the primer scheme used (if not specified, covid-artic-v3.json will
                        be used)
  --read_length READ_LENGTH
                        the read length (default value is 250)
  --depth DEPTH         the depth (default value is 500)
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
reference.fasta   reference_1.fastq reference_2.fastq
```

Next, let's create two Illumina paried Omicron FASTQ files using the default values for the read length (250) and depth (500).

```
$ gpas-covid-perfect-reads.py --variant_definitions ../variant_definitions/ --output omicron --variant_name Omicron --output omicron_r250_d500
$ ls omicron_r250_d500*
omicron_r250_d500.fasta   omicron_r250_d500_1.fastq omicron_r250_d500_2.fastq
```

Now we can increase the average depth 

```
$ gpas-covid-perfect-reads.py --variant_definitions ../variant_definitions/ --output omicron --variant_name Omicron --output omicron_r250_d1000 --depth 1000
$ ls omicron_r250_d1000*
omicron_r250_d1000.fasta   omicron_r250_d1000_1.fastq omicron_r250_d1000_2.fastq
```

These pairs of `fastq` files (after being compressed using `gzip`) can be used for 
* automated end-to-end testing of bioinformatics workflows, such as GPAS, since you can knows what the expected consensus sequence (the `fasta` file) should be.
* preliminary testing of new variants 

Philip W Fowler, 6 Dec 2021

