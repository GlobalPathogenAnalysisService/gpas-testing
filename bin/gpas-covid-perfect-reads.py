#!/usr/bin/env python3

import copy, glob, json, yaml, argparse, random

import numpy
import pyfastaq
import gumpy

import gpas_covid_perfect_reads as gcpr




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_definitions",required=True,help="the path to the variant_definitions repository/folder from phe-genomics ")
    parser.add_argument("--output",required=True,help="the stem of the output file")
    parser.add_argument("--variant_name",required=False,help="a JSON file specifying the mutations to apply to the covid reference (if not specified, you'll get a wildtype sequence)")
    parser.add_argument("--reference",required=False,default='config/MN908947.3.gbk',help="the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will be used)")
    parser.add_argument("--primer_definition",default='config/covid-artic-v3.json',help="the JSON file specifying the primer scheme used (if not specified, covid-artic-v3.json will be used)")
    parser.add_argument("--tech",default='illumina',help="whether to generate illumina (paired) or nanopore (unpaired) reads")
    parser.add_argument("--read_length",default=250,type=int,help="the read length in bases (default value is 250)")
    parser.add_argument("--read_stddev",default=0,type=int,help="the standard deviation in the read lengths (default value is 0)")
    parser.add_argument("--depth",default=500,type=int,help="the depth (default value is 500)")
    parser.add_argument("--error_rate",default=0.0,type=float,help="the percentage base error rate (default value is 0.0)")
    options = parser.parse_args()

    # load in the covid reference using gumpy
    covid_reference=gumpy.Genome(options.reference)

    assert 100>options.error_rate>=0

    error_rate=options.error_rate/100.

    bases={'A','T','C','G'}

    # load the definitions of the amplicons
    with open(options.primer_definition) as f:
        amplicons = json.load(f)

    # only if a variant has been specified, otherwise output reference
    if not options.variant_name:

        variant_definitions=False

        description = "Reference"

        genome=covid_reference.build_genome_string()

        index_lookup=covid_reference.nucleotide_index

    else:

        variant_definitions=gcpr.load_variant_definitions('../variant_definitions')

        # check that the variant that has been specified has a YAML file!
        assert options.variant_name in variant_definitions.keys(), "specified variant not defined here "+options.variant_definitions

        variant=gcpr.VariantGenome(covid_reference, variant_definitions[options.variant_name])

        description = options.variant_name
        #+"_"+options.primer_definition.split('.')[0]+"_readlength_"+str(options.read_length)+"_depth_"+str(options.depth)

        genome=variant.genome

        index_lookup=variant.index_lookup

    # # now read the FASTA file back in
    # refs = {}
    # pyfastaq.tasks.file_to_dict(options.output+".fasta", refs)
    # assert len(refs) == 1
    # ref = list(refs.values())[0]
    ref = pyfastaq.sequences.Fasta(id_in=description,seq_in=genome.upper())
    ref = ref.to_Fastq([40] * len(ref))

    if options.tech=='illumina':

        with open(options.output+"_1.fastq", "w") as f1, open(options.output+"_2.fastq", "w") as f2:

            for amplicon_name, amplicon_d in amplicons["amplicons"].items():

                # We'll make a read pair that looks like this:
                #
                #  |------ read1 ------------->
                #                      <--------- read2 -----|
                #  |-------------- amplicon -----------------|

                assert options.read_length < amplicon_d["end"] - amplicon_d["start"] + 1 < 2 * options.read_length

                start = numpy.where(index_lookup==amplicon_d["start"])[0][0]
                end = numpy.where(index_lookup==amplicon_d["end"])[0][0]

                for i in range(0, int(options.depth / 2)):

                    length = int(numpy.random.normal(options.read_length,options.read_stddev))

                    read1 = ref.subseq(start, start + length)
                    read2 = ref.subseq(end - length, end)
                    read2.revcomp()

                    read1.seq=mutate_read(read1.seq,error_rate)
                    read2.seq=mutate_read(read2.seq,error_rate)

                    read1.id = f"{amplicon_name}.{i} /1"
                    read2.id = f"{amplicon_name}.{i} /2"
                    print(read1, file=f1)
                    print(read2, file=f2)

                    read1.id = f"{amplicon_name}.{i}.2 /2"
                    read2.id = f"{amplicon_name}.{i}.2 /1"
                    print(read1, file=f2)
                    print(read2, file=f1)

    elif options.tech=='nanopore':

        with open(options.output+".fastq", "w") as f1:
            for amplicon_name, amplicon_d in amplicons["amplicons"].items():

                start = numpy.where(index_lookup==amplicon_d["start"])[0][0]
                end = numpy.where(index_lookup==amplicon_d["end"])[0][0]

                for i in range(0, int(options.depth / 2)):

                    length = int(numpy.random.normal(options.read_length,options.read_stddev))
                    read1 = ref.subseq(start, start + length)

                    read1.seq=mutate_read(read1.seq,error_rate)
                    read1.id = f"{amplicon_name}.{i} forward"
                    print(read1, file=f1)

                for i in range(0, int(options.depth / 2)):

                    length = int(numpy.random.normal(options.read_length,options.read_stddev))
                    read2 = ref.subseq(end - length, end)
                    read2.revcomp()

                    read2.seq=mutate_read(read2.seq,error_rate)
                    read2.id = f"{amplicon_name}.{i}.2 reverse"
                    print(read2, file=f1)

    else:
        raise ValueError("--tech must be one of illumina or nanopore!")
