# gpas_testing
Create perfect FASTQ files for the SARS-CoV-2 WHO lineages for use in testing

## Installation

First, let's get the repository and work in a virtual environment

```
git clone git@github.com:GenomePathogenAnalysisService/gpas-testing.git
cd gpas-testing
python -m venv env
source env/bin/activate
```

Now let's manually download the [`constellations`](https://github.com/cov-lineages/constellations) SARS-CoV-2 lineage definitions

```
git clone https://github.com/cov-lineages/constellations.git
```

Now we can automatically install the rest of dependencies

```
pip install -e .
```

I've specified the `-e` or `--editable` flag so that if you make changes to the `gpas-testing` they will flow through automatically into the installed version so you don't need to reinstall each time. This isn't required if you are using a static version.

Once complete, let's check the modest number of unit tests work

```
py.test tests/
```

This repository contains three scripts that help in creating and analysing batches of SARS-CoV-2 samples using GPAS.

### `gpas-synreads-covid-create.py` -- creating synthetic SARS-CoV-2 reads

```
usage: gpas-synreads-covid-create.py [-h] [--variant_definitions VARIANT_DEFINITIONS] [--pango_definitions PANGO_DEFINITIONS]
                                     [--output OUTPUT] [--variant_name VARIANT_NAME] [--reference REFERENCE] --tech TECH
                                     [--primers PRIMERS [PRIMERS ...]] [--read_length READ_LENGTH]
                                     [--read_stddev READ_STDDEV] [--depth DEPTH [DEPTH ...]] [--depth_stddev DEPTH_STDDEV]
                                     [--snps SNPS [SNPS ...]] [--repeats REPEATS] [--error_rate ERROR_RATE [ERROR_RATE ...]]
                                     [--drop_amplicons DROP_AMPLICONS [DROP_AMPLICONS ...]] [--write_fasta]
                                     [--bias_amplicons BIAS_AMPLICONS [BIAS_AMPLICONS ...]]
                                     [--bias_primers BIAS_PRIMERS [BIAS_PRIMERS ...]]
                                     [--drop_forward_amplicons DROP_FORWARD_AMPLICONS [DROP_FORWARD_AMPLICONS ...]]

optional arguments:
  -h, --help            show this help message and exit
  --variant_definitions VARIANT_DEFINITIONS
                        the path to the variant_definitions repository/folder from phe-genomics
  --pango_definitions PANGO_DEFINITIONS
                        the path to the constellations repository/folder from cov-lineages
  --output OUTPUT       the stem of the output file
  --variant_name VARIANT_NAME
                        the name of the variant, default is Reference
  --reference REFERENCE
                        the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will be used)
  --tech TECH           whether to generate illumina (paired) or nanopore (unpaired) reads
  --primers PRIMERS [PRIMERS ...]
                        the name of the primer schema, must be on of articv3, articv4, midnight1200, ampliseq
  --read_length READ_LENGTH
                        if specified, the read length in bases, otherwise defaults to the whole amplicon
  --read_stddev READ_STDDEV
                        the standard deviation in the read lengths (default value is 0)
  --depth DEPTH [DEPTH ...]
                        the depth (default value is 500)
  --depth_stddev DEPTH_STDDEV
                        the standard deviation of the depth distribution (default value is 0)
  --snps SNPS [SNPS ...]
                        the number of snps to randomly introduce into the sequence
  --repeats REPEATS     how many repeats to create
  --error_rate ERROR_RATE [ERROR_RATE ...]
                        the percentage base error rate (default value is 0.0)
  --drop_amplicons DROP_AMPLICONS [DROP_AMPLICONS ...]
                        the number (int) of one or more amplicons to drop i.e. have no reads.
  --write_fasta         whether to write out the FASTA file for the variant
  --bias_amplicons BIAS_AMPLICONS [BIAS_AMPLICONS ...]
                        whether to introduce an incorrect SNP in one or more specified amplicons
  --bias_primers BIAS_PRIMERS [BIAS_PRIMERS ...]
                        whether to introduce an incorrect SNP in both primers of an amplicon
  --drop_forward_amplicons DROP_FORWARD_AMPLICONS [DROP_FORWARD_AMPLICONS ...]
                        the names of one or more amplicons where there will be no reads mapping to the forward strand.
```                        

As described in Installation, you will have already installed the `constellation` definitions of (pangolin) SARS-CoV-2 lineages. The code also copes with the phe-genomics set, however at the time of writing these were not being updated as frequently. In addition, the definitions are not interchangeable e.g. if you want to create a genome that `pangolin` will classify as `BA.2` you need to use the `pangolin` definitions, which are `constellations`. Hence whilst both are offered, we currently recommend `constellations`.

First, we can simply create a set of perfect `cBA.1` Illumina reads 

```
gpas-synreads-covid-create.py --pango_definitions constellations/ --tech illumina --variant_name cBA.1 --write_fasta
ls illumina-*
illumina-articv3-cBA.1-cov-0snps-500d-0.0e-0.fasta   illumina-articv3-cBA.1-cov-0snps-500d-0.0e-0_2.fastq
illumina-articv3-cBA.1-cov-0snps-500d-0.0e-0_1.fastq
gzip illumina*fastq
```
Note that the `fastq` files were written uncompressed, so the first thing we need to do is `gzip` them. Also note that since we didn't tell the code what to call the files (using `--output`) it has given the files a hierarchial descriptive name. As you can see lot's of options have taken their default values, e.g. the default amplicon scheme is Artic v3, the default depth is 500 reads, there are no errors and no SNPs.

Let's make somthing a bit more realistic; a `cBA.3` sample that has 'been' sequenced with Nanopore using the Midnight primers. This sample has a mean depth of 500 and a standard deviation of 100 so there is a reasonable probability one or more amplicons only have a read depth of 250. One of the amplicons (2) has no reads associated with it and all the reads have a 2% error rate. This sample is 2 SNPs away from the `cBA.3` reference definition.

```
gpas-synreads-covid-create.py --pango_definitions constellations/ --tech nanopore --variant_name cBA.3 --write_fasta --depth 400 --read_stddev 100 --primers midnight1200 --snps 2 --error_rate 2 --drop_amplicons 2 --write_fasta
ls nanopore-*
nanopore-midnight1200-cBA.3-cov-2snps-400d-0.02e-a2da-0.fasta nanopore-midnight1200-cBA.3-cov-2snps-400d-0.02e-a2da-0.fastq
$ gzip nanopore*fastq
```

This degree of customisability let's you build a batch of samples to test "genetic space" e.g. one could create a batch containing all the common SARS-CoV-2 lineages in circulation, or a batch containing a range of error rates to explore at what point GPAS can no longer produce a consensus genome etc.

For demonstration purposes, let's build a simple batch of Illumina cBA.2 samples.

```
mkdir batch-0
cd batch-0
gpas-synreads-covid-create.py --pango_definitions constellations/ --tech illumina --variant_name cBA.2 --write_fasta --depth 300 --read_stddev 100 --primers articv4 --snps 2 --error_rate 1 --write_fasta --repeats 3
gzip *fastq
```

In this folder we have three pairs of gzipped FASTQ files we can upload and, crucially, the genome (FASTA file) that was used to create them, hence we have an *expectation* that we can compare back to later on and assert equality.

Now we can use the next script

## `gpas-build-uploadcsv.py` - automatically creating an upload CSV

One can always write manually an upload CSV but is time-consuming, boring and can introduce mistakes, hence this helper script. If you run it in a folder that contains gzipped fastq files (like `batch-0`) and give it some arguments it will create a valid upload CSV for you that you can then use either with the Electron Client / App or the CLI to upload the batch to GPAS for processing.

```
gpas-build-uploadcsv.py --help
usage: gpas-build-uploadcsv.py [-h] [--country COUNTRY] [--tech TECH] [--file_type FILE_TYPE] [--tag_file TAG_FILE]
                               [--uuid_length UUID_LENGTH] [--number_of_tags NUMBER_OF_TAGS] [--old_format]
                               [--organisation ORGANISATION]

optional arguments:
  -h, --help            show this help message and exit
  --country COUNTRY     the name of the country where the samples were collected
  --tech TECH           whether to generate illumina (paired) or nanopore (unpaired) reads
  --file_type FILE_TYPE
                        whether to look for FASTQ or BAM files
  --tag_file TAG_FILE   a plaintext file with one row per tag
  --uuid_length UUID_LENGTH
                        whether to use a long or short UUID4
  --number_of_tags NUMBER_OF_TAGS
                        how many tags to give each sample. Can be zero, or up to the number of rows in <tag_file>. Default is 2 so as to test the delimiter
  --old_format          whether to use the original upload CSV format and headers
  --organisation ORGANISATION
                        the name of the organisation (the user must belong to it otherwise validation will fail)
```

We will only focus here on the new (May 2022 / Electron Client v1.0.7 / gpas-uploader v1.1.1) format for the upload CSV so the last two arguments will not concern us. First we need to make a plaintext file containing each tag we wish the code to pick from on its own line. I've made one called `tags.txt` which has a single line containing `test`.

Given the sequencing technology used, the code will detect the samples where it is run and then autogenerate e.g.  a local `batch` and `run_number` for each sample as well as randomly create a `collection_date` from the recent past. It does not add any `control`, `region` or `district`. You can always edit the file afterwards to customise it (e.g. by adding deliberate mistakes if you wish to test the validation rules in the uploader). It also renames all the fastq files by appending either a long or short `uuid` string. This ensures that each upload CSV is different, even if it contains the "same" samples. The result is printed to STDOUT so you need to redirect it to a file if you wish to keep it.

```
gpas-build-uploadcsv.py --country USA --tag_file tags.txt --uuid_length short --number_of_tags 1 --tech illumina > upload.csv
```
You can now use this upload CSV either with the Electron Client or the gpas cli.

## `gpas-synreads-covid-create.py` - checking if the GPAS consensus genome matches what the reads were built from

This script compares the output consensus genomes to the genome used to build the FASTQ files and will declare success if the input genome contains the output genome and the latter is longer than 29k bases. It also compares the input lineage to the detected lineage. Since it requires the sample to be run through GPAS I won't describe it more here but it does assume you have the file mapping the local to gpas identifiers, which here I will call `sample_names.csv`.

## Cultured "truth" samples

These can be obtained from these two ENA projects: [PRJEB50520](https://www.ebi.ac.uk/ena/browser/view/PRJEB50520) and [PRJEB51850](https://www.ebi.ac.uk/ena/browser/view/PRJEB51850) with `vcf` files describing the variation w.r.t the Wuhan reference found [here](https://github.com/iqbal-lab-org/covid-truth-datasets).

## `tb-synreads` - Generate TB synthetic reads from a specified file
Use a file which specifies variants to generate synthetic paired FASTQ files. 
```
usage: tb-synreads [-h] --reference REFERENCE --depth DEPTH --read_length READ_LENGTH --variant_file VARIANT_FILE [--error_rate ERROR_RATE] --output OUTPUT

options:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        Path to the reference genome
  --depth DEPTH         Depth of the reads
  --read_length READ_LENGTH
                        Read length
  --variant_file VARIANT_FILE
                        Path to a file detailing mutations
  --error_rate ERROR_RATE
                        The percentage base error rate as a decimal i.e in range [0,1]
  --output OUTPUT       Path to the stem of the VCF
```
File specification
```
#This is a comment

#Specify a SNP
<pos> <ref> <alt>

#Specify an insertion
<pos> ins <bases>

#Specify a deletion
<pos> del <number of bases>
```
Has been developed so that the output can be piped to timing utilities: `tb-synreads ... | ts '[%H:%M:%.S]` for timestamped output

## Docker Use

There is a Dockerfile provided to make your own container or an image available in Dockerhub (oxfordmmm/gpas_testing). Both versions put the `constellations` and `variant_definitions` directories at the root of the container.

(PWF: I've messed with this so it needs updating..)

The container can be run with a command such as:

```
docker run -v /path/to/output:/output oxfordmmm/gpas_testing python3 /gpas_testing/bin/gpas_testing.py  --pango_definitions /constellations/ --output /output/reference --tech illumina --variant_name reference --write_fasta 

docker run -v .:/output oxfordmmm/gpas_testing python3 /gpas_testing/bin/gpas-synreads-covid-create.py --pango_definitions /constellations/ --tech illumina --variant_name cBA.1 --write_fasta --output /output/BA.1
```

Philip W Fowler, 7 July 2022

