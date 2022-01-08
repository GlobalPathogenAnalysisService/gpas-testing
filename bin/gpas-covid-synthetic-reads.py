#!/usr/bin/env python3

import copy, glob, argparse, random, pkg_resources

import numpy, yaml, pyfastaq, pandas
from tqdm import tqdm
import gumpy

import gpas_covid_synthetic_reads as gcsr

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_definitions",required=True,help="the path to the variant_definitions repository/folder from phe-genomics ")
    parser.add_argument("--output",required=True,help="the stem of the output file")
    parser.add_argument("--variant_name",default='Reference',help="the name of the variant, default is Reference")
    parser.add_argument("--reference",required=False,default=pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/MN908947.3.gbk'),help="the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will be used)")
    parser.add_argument("--primers",default='artic-v3',help="the name of the primer schema, must be on of artic-v3, artic-v4, midnight-1200")
    parser.add_argument("--tech",default='illumina',help="whether to generate illumina (paired) or nanopore (unpaired) reads")
    parser.add_argument("--read_length",default=250,type=int,help="the read length in bases (default value is 250)")
    parser.add_argument("--read_stddev",default=0,type=int,help="the standard deviation in the read lengths (default value is 0)")
    parser.add_argument("--depth",default=500,type=int,help="the depth (default value is 500)")
    parser.add_argument("--snps",nargs='+',default=0,type=int,help="the number of snps to randomly introduce into the sequence")
    parser.add_argument("--repeats",default=1,type=int,help="whether to repeat building the FASTQ files")
    parser.add_argument("--error_rate",default=0.0,type=float,help="the percentage base error rate (default value is 0.0)")
    parser.add_argument("--write_fasta", dest="write_fasta",action="store_true", help="whether to write out the FASTA file for the variant")
    options = parser.parse_args()

    # load in the covid reference using gumpy
    covid_reference=gumpy.Genome(options.reference)

    assert 100>options.error_rate>=0

    error_rate=options.error_rate/100.

    bases={'A','T','C','G'}

    assert options.primers in ['articv3','articv4','midnight1200']

    if options.primers=='articv3':
        primer_scheme_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/artic-v3.qcovid.tsv')
    elif options.primers=='articv4':
        primer_scheme_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/artic-v4.qcovid.tsv')
    elif options.primers=='midnight1200':
        primer_scheme_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/midnight-1200.qcovid.tsv')

    # load the definitions of the amplicons
    primers=pandas.read_csv(primer_scheme_file,\
                        sep='\t',
                        names=['pool','name','seq','bool1','bool2','start'])

    def assign_end(row):
        return(row['start']+len(row.seq)-1)

    primers['end']=primers.apply(assign_end,axis=1)

    def find_handedness(row):
        return(row['name'].split('_')[-1])

    primers['hand']=primers.apply(find_handedness,axis=1)

    left=primers.loc[primers.hand=='LEFT']
    right=primers.loc[primers.hand=='RIGHT']

    amplicons=left.merge(right,left_on='pool',right_on='pool',suffixes=['_left','_right'])
    amplicons=amplicons[['pool','start_left','end_right']]
    amplicons.rename(columns={'pool':'name','start_left':'start','end_right':'end'},inplace=True)

    # only if a variant has been specified, otherwise output reference
    if options.variant_name == 'Reference':

        variant_definitions=False

        variant=copy.deepcopy(covid_reference)

        description = "Reference"

        index_lookup=variant.nucleotide_index

    else:

        variant_definitions=gcsr.load_variant_definitions(options.variant_definitions)

        # check that the variant that has been specified has a YAML file!
        assert options.variant_name in variant_definitions.keys(), "specified variant not defined here "+options.variant_definitions

        variant=gcsr.VariantGenome(covid_reference, variant_definitions[options.variant_name])

        description = options.variant_name
        #+"_"+options.primer_definition.split('.')[0]+"_readlength_"+str(options.read_length)+"_depth_"+str(options.depth)

        index_lookup=variant.index_lookup

        variant=variant.variant

    if options.snps>0:

        snp_indices=numpy.random.choice(variant.nucleotide_index,size=options.snps,replace=False)

        for i in snp_indices:
            mask=variant.nucleotide_index==i
            current_base=variant.nucleotide_sequence[mask]
            bases={'a','t','c','g'}
            new_bases = bases ^ set(current_base)
            new_base=random.choice(list(new_bases))
            variant.nucleotide_sequence[mask]=new_base

    genome=variant.build_genome_string()

    for repeat in tqdm(range(options.repeats)):

        if options.repeats==1:
            outputstem=options.output
        else:
            outputstem=options.output+"-"+str(repeat)


        if options.write_fasta:
            variant.save_fasta(outputstem+".fasta",\
                                fixed_length=False,\
                                overwrite_existing=True,\
                                description=description)

        ref = pyfastaq.sequences.Fasta(id_in=description,seq_in=genome.upper())
        ref = ref.to_Fastq([40] * len(ref))

        if options.tech=='illumina':

            with open(outputstem+"_1.fastq", "w") as f1, open(outputstem+"_2.fastq", "w") as f2:

                for idx,row in amplicons.iterrows():

                    amplicon_name = row['name']
                    # We'll make a read pair that looks like this:
                    #
                    #  |------ read1 ------------->
                    #                      <--------- read2 -----|
                    #  |-------------- amplicon -----------------|

                    assert options.read_length < row["end"] - row["start"] + 1 < 2 * options.read_length

                    start = numpy.where(index_lookup==row["start"])[0][0]
                    end = numpy.where(index_lookup==row["end"])[0][0]

                    for i in range(0, int(options.depth / 2)):

                        length = int(numpy.random.normal(options.read_length,options.read_stddev))

                        read1 = ref.subseq(start, start + length)
                        read2 = ref.subseq(end - length, end)
                        read2.revcomp()

                        if error_rate>0:
                            read1.seq=gcsr.mutate_read(read1.seq,error_rate=error_rate)
                            read2.seq=gcsr.mutate_read(read2.seq,error_rate=error_rate)

                        read1.id = f"{amplicon_name}.{i} /1"
                        read2.id = f"{amplicon_name}.{i} /2"
                        print(read1, file=f1)
                        print(read2, file=f2)

                        read1.id = f"{amplicon_name}.{i}.2 /2"
                        read2.id = f"{amplicon_name}.{i}.2 /1"
                        print(read1, file=f2)
                        print(read2, file=f1)

        elif options.tech=='nanopore':

            with open(outputstem+".fastq", "w") as f1:

                for idx,row in amplicons.iterrows():

                    amplicon_name = row['name']

                    start = numpy.where(index_lookup==row["start"])[0][0]
                    end = numpy.where(index_lookup==row["end"])[0][0]

                    for i in range(0, int(options.depth / 2)):

                        length = int(numpy.random.normal(options.read_length,options.read_stddev))
                        read1 = ref.subseq(start, start + length)

                        if error_rate>0:
                            read1.seq=gcsr.mutate_read(read1.seq,error_rate=error_rate)

                        read1.id = f"{amplicon_name}.{i} forward"
                        print(read1, file=f1)

                    for i in range(0, int(options.depth / 2)):

                        length = int(numpy.random.normal(options.read_length,options.read_stddev))
                        read2 = ref.subseq(end - length, end)
                        read2.revcomp()

                        if error_rate>0:
                            read2.seq=gcsr.mutate_read(read2.seq,error_rate=error_rate)

                        read2.id = f"{amplicon_name}.{i}.2 reverse"
                        print(read2, file=f1)

        else:
            raise ValueError("--tech must be one of illumina or nanopore!")
