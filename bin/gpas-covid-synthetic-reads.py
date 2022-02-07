#!/usr/bin/env python3

import copy, glob, argparse, random, pkg_resources

import numpy, yaml, pyfastaq, pandas
import gumpy

import gpas_covid_synthetic_reads as gcsr

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_definitions",default=False,help="the path to the variant_definitions repository/folder from phe-genomics ")
    parser.add_argument("--pango_definitions",default=False,help="the path to the constellations repository/folder from cov-lineages ")
    parser.add_argument("--output",required=False,help="the stem of the output file")
    parser.add_argument("--variant_name",default='Reference',help="the name of the variant, default is Reference")
    parser.add_argument("--reference",required=False,default=pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/MN908947.3.gbk'),help="the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will be used)")
    parser.add_argument("--tech",required=True,help="whether to generate illumina (paired) or nanopore (unpaired) reads")
    parser.add_argument("--primers",nargs='+',default=['articv3'],help="the name of the primer schema, must be on of articv3, articv4, midnight1200, ampliseq")
    parser.add_argument("--read_length",default=None,type=int,help="if specified, the read length in bases, otherwise defaults to the whole amplicon")
    parser.add_argument("--read_stddev",default=0,type=int,help="the standard deviation in the read lengths (default value is 0)")
    parser.add_argument("--depth",nargs='+',default=[500],type=int,help="the depth (default value is 500)")
    parser.add_argument("--snps",nargs='+',default=[0],type=int,help="the number of snps to randomly introduce into the sequence")
    parser.add_argument("--repeats",default=1,type=int,help="how many repeats to create")
    parser.add_argument("--error_rate",nargs='+',default=[0.0],type=float,help="the percentage base error rate (default value is 0.0)")
    parser.add_argument("--drop_amplicons", nargs='+', type=int, required=False,help="the number (int) of one or more amplicons to drop i.e. have no reads.")
    parser.add_argument("--write_fasta", action="store_true", help="whether to write out the FASTA file for the variant")
    parser.add_argument("--bias_amplicons", nargs='+', type=int, default=None, help="whether to introduce an incorrect SNP in one or more specified amplicons")
    parser.add_argument("--bias_primers", nargs='+', type=int, default=None, help="whether to introduce an incorrect SNP in both primers of an amplicon")
    parser.add_argument("--drop_forward_amplicons", nargs='+', type=int, required=False,help="the names of one or more amplicons where there will be no reads mapping to the forward strand.")
    options = parser.parse_args()

    # load in the covid reference using gumpy
    covid_reference=gumpy.Genome(options.reference)

    error_rates=numpy.array(options.error_rate)/100.

    bases = {'a','t','c','g'}

    assert options.tech in ['illumina', 'nanopore']

    assert not (options.variant_definitions and options.pango_definitions), 'cannot specify both variant_definitions and pango_definitions'

    if options.drop_amplicons is not None:
        assert len(options.primers)==1, 'can only specify dropped amplicons for a single primer scheme!'

    for primer_name in options.primers:
        assert primer_name in ['articv3','articv4','midnight1200', 'ampliseq']

        if primer_name=='articv3':
            amplicon_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/covid-artic-v3.amplicons.csv')
        elif primer_name=='articv4':
            amplicon_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/covid-artic-v4.amplicons.csv')
        elif primer_name=='midnight1200':
            amplicon_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/covid-midnight-1200.amplicons.csv')
        elif primer_name=='ampliseq':
            amplicon_file = pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/covid-ampliseq-v1.amplicons.csv')

        # load the definitions of the amplicons
        amplicons = pandas.read_csv(amplicon_file,\
                                    dtype = {'start': int, 'end': int, 'length': int, 'number': int,
                                             'start_left': int, 'end_left': int,
                                             'start_right': int, 'end_right': int,
                                             'start_amplicon': int, 'end_amplicon': int})

        if options.drop_amplicons is not None:
            for i in options.drop_amplicons:
                assert i in list(amplicons['number']), 'amplicon '+i+' not found in '+primer_name+' scheme'

        if options.drop_forward_amplicons is not None and options.drop_amplicons is not None:
            for i in options.drop_forward_amplicons:
                assert i not in options.drop_amplicons, 'cannot both drop an amplicon andd drop its forward strand'
                assert i in list(amplicons['number']), 'amplicon '+i+' not found in '+primer_name+' scheme'

        if options.bias_amplicons is not None and options.drop_amplicons is not None:
            for i in options.bias_amplicons:
                assert i not in options.drop_amplicons, 'cannot both bias and drop an amplicon'
                assert i in list(amplicons['number']), 'amplicon '+i+' not found in '+primer_name+' scheme'

        if options.bias_primers is not None and options.drop_amplicons is not None:
            for i in options.bias_primers:
                assert i not in options.drop_amplicons, 'cannot both bias and drop an amplicon'
                assert i in list(amplicons['number']), 'amplicon '+i+' not found in '+primer_name+' scheme'

        if options.variant_definitions:

            variant_definitions = gcsr.load_variant_definitions(options.variant_definitions)

            # check that the variant that has been specified has a YAML file!
            if options.variant_name != 'reference':
                assert options.variant_name in variant_definitions.keys(), "specified variant not defined here "+options.variant_definitions

            variant = gcsr.VariantGenome(covid_reference,
                                         variant_definitions,
                                         options.variant_name)

            variant_source = 'phe'

        elif options.pango_definitions:

            variant_definitions = gcsr.load_pango_definitions(options.pango_definitions)

            # check that the variant that has been specified has a JSON file!
            if options.variant_name != 'reference':
                assert options.variant_name in variant_definitions.keys(), "specified variant not defined here "+options.variant_definitions

            variant = gcsr.PangoGenome(covid_reference,
                                       variant_definitions,
                                       options.variant_name)

            variant_source = 'cov'

        else:
            raise ValueError('must specify either --variant_definitions or --pango_definitions')


        description = options.variant_name
        #+"_"+options.primer_definition.split('.')[0]+"_readlength_"+str(options.read_length)+"_depth_"+str(options.depth)

        index_lookup=variant.index_lookup

        # drop any amplicons specified
        if options.drop_amplicons is not None:

            for chosen_amplicon in options.drop_amplicons:

                row = amplicons[amplicons.number == chosen_amplicon]

                mask = (covid_reference.nucleotide_index >= int(row['start_amplicon'])) & (covid_reference.nucleotide_index <= int(row['end_amplicon']))

                variant.expected.nucleotide_sequence[mask] = 'n'

        if options.drop_forward_amplicons is not None:

            for chosen_amplicon in options.drop_forward_amplicons:

                row = amplicons[amplicons.number == chosen_amplicon]

                mask = (covid_reference.nucleotide_index > int(row['start_left'])) & (covid_reference.nucleotide_index < int(row['end_right']))

                variant.expected.nucleotide_sequence[mask] = 'n'

        if options.bias_amplicons is not None:

            for chosen_amplicon in options.bias_amplicons:

                row = amplicons[amplicons.number == chosen_amplicon]

                idx = numpy.random.choice(numpy.arange(int(row['start_amplicon']),int(row['end_amplicon'])))
                mask = covid_reference.nucleotide_index == idx
                ref_base = covid_reference.nucleotide_sequence[mask]
                new_bases = bases ^ set(ref_base)
                new_base = random.choice(list(new_bases))

                variant.input1.nucleotide_sequence[mask] = new_base
                variant.input2.nucleotide_sequence[mask] = ref_base
                variant.expected.nucleotide_sequence[mask] = 'n'

        if options.bias_primers is not None:

            for chosen_amplicon in options.bias_primers:

                row = amplicons[amplicons.number == chosen_amplicon]

                idx = numpy.random.choice(numpy.arange(int(row['start_left']),int(row['end_left'])))
                mask = covid_reference.nucleotide_index == idx
                ref_base = covid_reference.nucleotide_sequence[mask]
                new_bases = bases ^ set(ref_base)
                new_base = random.choice(list(new_bases))
                variant.input1.nucleotide_sequence[mask] = new_base
                variant.input2.nucleotide_sequence[mask] = ref_base
                variant.expected.nucleotide_sequence[mask] = 'n'

                idx = numpy.random.choice(numpy.arange(int(row['start_right']),int(row['end_right'])))
                mask = covid_reference.nucleotide_index == idx
                ref_base = covid_reference.nucleotide_sequence[mask]
                new_bases = bases ^ set(ref_base)
                new_base = random.choice(list(new_bases))
                variant.input1.nucleotide_sequence[mask] = new_base
                variant.input2.nucleotide_sequence[mask] = ref_base
                variant.expected.nucleotide_sequence[mask] = 'n'

        for snps in options.snps:

            current_variant = copy.deepcopy(variant)

            if snps > 0:

                snp_indices = numpy.random.choice(current_variant.expected.nucleotide_index,size=snps,replace=False)

                for i in snp_indices:
                    mask = current_variant.expected.nucleotide_index == i
                    current_base = current_variant.expected.nucleotide_sequence[mask]

                    # don't mutate an N
                    if current_base == 'n':
                        continue

                    new_bases = bases ^ set(current_base)
                    new_base = random.choice(list(new_bases))
                    current_variant.expected.nucleotide_sequence[mask] = new_base
                    current_variant.input1.nucleotide_sequence[mask] = new_base
                    current_variant.input2.nucleotide_sequence[mask] = new_base

            input1_genome = current_variant.input1.build_genome_string()
            input1_genome = pyfastaq.sequences.Fasta(id_in=current_variant.name,
                                                     seq_in=input1_genome.upper())
            input1_genome = input1_genome.to_Fastq([40] * len(input1_genome))

            input2_genome = current_variant.input2.build_genome_string()
            input2_genome = pyfastaq.sequences.Fasta(id_in=current_variant.name,
                                                     seq_in=input2_genome.upper())
            input2_genome = input2_genome.to_Fastq([40] * len(input2_genome))

            expected_genome = current_variant.expected.build_genome_string()
            expected_genome = pyfastaq.sequences.Fasta(id_in=current_variant.name,
                                                    seq_in=expected_genome.upper())
            expected_genome = expected_genome.to_Fastq([40] * len(expected_genome))

            for depth in options.depth:

                for error_rate in error_rates:

                    assert 100>error_rate>=0

                    for repeat in range(options.repeats):

                        if options.output is None:
                            outputstem = options.tech + '-' + primer_name.lower() + '-' + options.variant_name + '-' +\
                                         variant_source + '-' + str(snps) + 'snps-' + str(depth) + 'd-' +\
                                         str(error_rate) + 'e-'
                            if options.drop_amplicons is not None:
                                outputstem += ''.join('a' + str(i) for i in options.drop_amplicons)
                                outputstem += 'da-'
                            if options.drop_forward_amplicons is not None:
                                outputstem += ''.join('a' + str(i) for i in options.drop_forward_amplicons)
                                outputstem += 'df-'
                            if options.bias_amplicons is not None:
                                outputstem += ''.join('a' + str(i) for i in options.bias_amplicons)
                                outputstem += 'ba-'
                            if options.bias_primers is not None:
                                outputstem += ''.join('a' + str(i) for i in options.bias_primers)
                                outputstem += 'bp-'
                            outputstem += str(repeat)
                        else:
                            outputstem = options.output

                        if options.write_fasta:

                            current_variant.expected.save_fasta(outputstem+".fasta",\
                                                        fixed_length=False,\
                                                        overwrite_existing=True,\
                                                        description=variant.name)

                        if options.tech=='illumina':

                            with open(outputstem+"_1.fastq", "w") as f1, open(outputstem+"_2.fastq", "w") as f2:

                                for idx,row in amplicons.iterrows():

                                    amplicon_number = row['number']

                                    if options.drop_amplicons is not None and amplicon_number in options.drop_amplicons:
                                        continue
                                    # We'll make a read pair that looks like this:
                                    #
                                    #  |------ read1 ------------->
                                    #                      <--------- read2 -----|
                                    #  |-------------- amplicon -----------------|
                                    if options.read_length is not None:
                                        assert options.read_length < row["end"] - row["start"] + 1 < 2 * options.read_length

                                    start = numpy.where(index_lookup == row["start"])[0][0]
                                    end = numpy.where(index_lookup == row["end"])[0][0]

                                    for i in range(0, int(depth / 2)):

                                        if options.read_length is None:
                                            length = end - start
                                        else:
                                            length = int(numpy.random.normal(options.read_length, options.read_stddev))

                                        if options.bias_amplicons is not None or options.bias_primers is not None:
                                            if (amplicon_number in options.bias_amplicons) or (amplicon_number in options.bias_primers):
                                                read1 = input2_genome.subseq(start, start + length)
                                                read2 = input2_genome.subseq(end - length, end)
                                            else:
                                                read1 = input1_genome.subseq(start, start + length)
                                                read2 = input1_genome.subseq(end - length, end)
                                        elif options.drop_forward_amplicons is not None:
                                            read1 = input1_genome.subseq(start, start + length)
                                            read2 = input1_genome.subseq(end - length, end)
                                        else:
                                            read1 = expected_genome.subseq(start, start + length)
                                            read2 = expected_genome.subseq(end - length, end)

                                        # reverse complement the second read
                                        read2.revcomp()

                                        if error_rate > 0:
                                            read1.seq = gcsr.mutate_read(read1.seq, error_rate=error_rate)
                                            read2.seq = gcsr.mutate_read(read2.seq, error_rate=error_rate)

                                        read1.id = f"{amplicon_number}.{i} /1"
                                        read2.id = f"{amplicon_number}.{i} /2"

                                        if options.drop_forward_amplicons is not None:
                                             if amplicon_number not in options.drop_forward_amplicons:
                                                 print(read1, file=f1)
                                        else:
                                            print(read1, file=f1)
                                        print(read2, file=f2)

                                        read1.id = f"{amplicon_number}.{i}.2 /2"
                                        read2.id = f"{amplicon_number}.{i}.2 /1"
                                        print(read1, file=f2)
                                        if options.drop_forward_amplicons is not None:
                                            if amplicon_number not in options.drop_forward_amplicons:
                                                print(read2, file=f1)
                                        else:
                                                print(read2, file=f1)


                        elif options.tech == 'nanopore':

                            with open(outputstem+".fastq", "w") as f1:

                                for idx, row in amplicons.iterrows():

                                    amplicon_number = row['number']

                                    if options.drop_amplicons is not None and amplicon_number in options.drop_amplicons:
                                        continue

                                    start = numpy.where(index_lookup == row["start"])[0][0]
                                    end = numpy.where(index_lookup == row["end"])[0][0]

                                    for i in range(0, int(depth / 2)):

                                        if options.drop_forward_amplicons is not None and amplicon_number in options.drop_forward_amplicons:
                                            continue

                                        if options.read_length is None:
                                            length = end - start
                                        else:
                                            length = int(numpy.random.normal(options.read_length, options.read_stddev))

                                        if options.bias_amplicons is not None or options.bias_primers is not None:
                                            if (amplicon_number in options.bias_amplicons) or (amplicon_number in options.bias_primers):
                                                read1 = input2_genome.subseq(start, start + length)
                                            else:
                                                read1 = input1_genome.subseq(start, start + length)
                                        else:
                                            read1 = expected_genome.subseq(start, start + length)

                                        if error_rate > 0:
                                            read1.seq = gcsr.mutate_read(read1.seq, error_rate=error_rate)

                                        read1.id = f"{amplicon_number}.{i} forward"
                                        print(read1, file=f1)

                                    for i in range(0, int(depth / 2)):

                                        if options.read_length is None:
                                            length = end - start
                                        else:
                                            length = int(numpy.random.normal(options.read_length, options.read_stddev))

                                        if options.bias_amplicons is not None or options.bias_primers is not None:
                                            if (amplicon_number in options.bias_amplicons) or (amplicon_number in options.bias_primers):
                                                read2 = input2_genome.subseq(end - length, end)
                                            else:
                                                read2 = input1_genome.subseq(end - length, end)
                                        elif options.drop_forward_amplicons is not None:
                                            read2 = input1_genome.subseq(end - length, end)
                                        else:
                                            read2 = expected_genome.subseq(end - length, end)

                                        read2.revcomp()

                                        if error_rate > 0:
                                            read2.seq = gcsr.mutate_read(read2.seq, error_rate=error_rate)

                                        read2.id = f"{amplicon_number}.{i}.2 reverse"
                                        print(read2, file=f1)

                        else:
                            raise ValueError("--tech must be one of illumina or nanopore!")
