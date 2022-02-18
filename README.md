# gpas-covid-synthetic-reads
Create perfect FASTQ files for the SARS-CoV-2 WHO lineages for use in testing

## Installation

You first need to install [`gumpy`](https://github.com/oxfordmmm/gumpy) as the code uses gumpy to build the genome of each sample.

Once complete, 

```
$ git clone https://github.com/GenomePathogenAnalysisService/gpas-covid-synthetic-reads.git
$ cd gpas-covid-synthetic-reads
$ python3 setup.py install --user
```

Check it has installed correctly

```
$ usage: gpas-covid-synthetic-reads.py [-h]
                                     [--variant_definitions VARIANT_DEFINITIONS]
                                     [--pango_definitions PANGO_DEFINITIONS]
                                     [--output OUTPUT]
                                     [--variant_name VARIANT_NAME]
                                     [--reference REFERENCE] --tech TECH
                                     [--primers PRIMERS [PRIMERS ...]]
                                     [--read_length READ_LENGTH]
                                     [--read_stddev READ_STDDEV]
                                     [--depth DEPTH [DEPTH ...]]
                                     [--snps SNPS [SNPS ...]]
                                     [--repeats REPEATS]
                                     [--error_rate ERROR_RATE [ERROR_RATE ...]]
                                     [--drop_amplicons DROP_AMPLICONS [DROP_AMPLICONS ...]]
                                     [--write_fasta]
                                     [--bias_amplicons BIAS_AMPLICONS [BIAS_AMPLICONS ...]]
                                     [--bias_primers BIAS_PRIMERS [BIAS_PRIMERS ...]]
                                     [--drop_forward_amplicons DROP_FORWARD_AMPLICONS [DROP_FORWARD_AMPLICONS ...]]

optional arguments:
  -h, --help            show this help message and exit
  --variant_definitions VARIANT_DEFINITIONS
                        the path to the variant_definitions repository/folder
                        from phe-genomics
  --pango_definitions PANGO_DEFINITIONS
                        the path to the constellations repository/folder from
                        cov-lineages
  --output OUTPUT       the stem of the output file
  --variant_name VARIANT_NAME
                        the name of the variant, default is Reference
  --reference REFERENCE
                        the GenBank file of the covid reference (if not
                        specified, the MN908947.3.gbk reference will be used)
  --tech TECH           whether to generate illumina (paired) or nanopore
                        (unpaired) reads
  --primers PRIMERS [PRIMERS ...]
                        the name of the primer schema, must be on of articv3,
                        articv4, midnight1200, ampliseq
  --read_length READ_LENGTH
                        if specified, the read length in bases, otherwise
                        defaults to the whole amplicon
  --read_stddev READ_STDDEV
                        the standard deviation in the read lengths (default
                        value is 0)
  --depth DEPTH [DEPTH ...]
                        the depth (default value is 500)
  --snps SNPS [SNPS ...]
                        the number of snps to randomly introduce into the
                        sequence
  --repeats REPEATS     how many repeats to create
  --error_rate ERROR_RATE [ERROR_RATE ...]
                        the percentage base error rate (default value is 0.0)
  --drop_amplicons DROP_AMPLICONS [DROP_AMPLICONS ...]
                        the number (int) of one or more amplicons to drop i.e.
                        have no reads.
  --write_fasta         whether to write out the FASTA file for the variant
  --bias_amplicons BIAS_AMPLICONS [BIAS_AMPLICONS ...]
                        whether to introduce an incorrect SNP in one or more
                        specified amplicons
  --bias_primers BIAS_PRIMERS [BIAS_PRIMERS ...]
                        whether to introduce an incorrect SNP in both primers
                        of an amplicon
  --drop_forward_amplicons DROP_FORWARD_AMPLICONS [DROP_FORWARD_AMPLICONS ...]
                        the names of one or more amplicons where there will be
                        no reads mapping to the forward strand.
```

## Usage

You need to install a set of SARS-CoV-2 variant definitions. Two flavours are supported: `constellation` from `cov-lineages` which is the team behind `pangolin`, or `variant_definitions` from `phe-genomics` which support `aln2type`. We originally used the latter as the description of the genetic variation found in each lineage is easier to use, but we then found that e.g. an Omicron as defined by `variant-definitions` was not classified as such by `pangolin` and `constellations` appears to be more rapidly updated. Hence whilst both are offered, we currently recommend `constellations`.

To install [`variant_definitions`](https://github.com/phe-genomics/variant_definitions) curated by `phe-genomics` in an appropriate place in your filesystem

```
$ git clone https://github.com/phe-genomics/variant_definitions
```
Whilst to install [`constellations`](https://github.com/cov-lineages/constellations) curated by `cov-lineages` issue

```
$ git clone https://github.com/cov-lineages/constellations.git
```

First, we can simply create a set of perfect Illumina reads for the SARS-CoV-2 reference

```
$ gpas-covid-synthetic-reads.py --pango_definitions ../constellations/ --output reference --tech illumina --variant_name reference --write_fasta
$ ls -lrt reference*
reference.fasta   reference_1.fastq reference_2.fastq
$ gzip *fastq
```

Next, let's create two Illumina paried Omicron FASTQ files using the default values for the read length (250) and depth (500).

```
$ gpas-covid-synthetic-reads.py --pango_definitions ../constellations/ --output omicron --tech illumina --variant_name omicron --write_fasta
$ ls omicron*
omicron.fasta   omicron_1.fastq omicron_2.fastq
$ gzip omicron*fastq
```

Now we can increase the average depth from the default of 500 to 1000

```
$ gpas-covid-synthetic-reads.py --pango_definitions ../constellations/ --output omicron_d1000 --tech illumina --variant_name omicron --write_fasta --depth 1000
$ ls omicron_d1000*
omicron_d1000.fasta   omicron_d1000_1.fastq omicron_d1000_2.fastq
$ gzip omicron*fastq
```

You can create a larger set of different, synthetic samples for testing using the `--snps` or `--error_rate` flags. Both are stochastic so repeats will produce different samples. Note that the synthetic samples are about 10-100x smaller than real samples so cannot be used for true benchmarking.

These pairs of `fastq` files (after being compressed using `gzip`) can be used for 
* automated end-to-end testing of bioinformatics workflows, such as GPAS, since you can knows what the expected consensus sequence (the `fasta` file) should be.
* preliminary testing of new variants 

Philip W Fowler, 18 Feb 2022

