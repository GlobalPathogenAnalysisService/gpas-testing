#!/usr/bin/env python3

import copy, glob, argparse, random, pkg_resources, pathlib

import numpy, pyfastaq, pandas
import gumpy

from tqdm import tqdm

import gpas_testing

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--output",required=False,help="the stem of the output file")
    parser.add_argument("--reference",required=False,default=pkg_resources.resource_filename("gpas_testing", 'data/NC_000962.3.gbk.gz'),help="the GenBank file of the covid reference (if not specified, the NC_000962.3.gbk reference will be used)")
    parser.add_argument("--tech",required=True,help="whether to generate illumina (paired) or nanopore (unpaired) reads")
    parser.add_argument("--read_length",default=None,type=int,help="if specified, the read length in bases, otherwise defaults to the whole amplicon")
    parser.add_argument("--read_stddev",default=0,type=int,help="the standard deviation in the read lengths (default value is 0)")
    parser.add_argument("--depth",default=50,type=int,help="the depth (default value is 50)")
    parser.add_argument("--snps",nargs='+',default=[0],type=int,help="the number of snps to randomly introduce into the sequence")
    parser.add_argument("--variant_file",required=False,type=str,help="an optional file containing one or more genetic variants to add to the sample")
    parser.add_argument("--repeats",default=1,type=int,help="how many repeats to create")
    parser.add_argument("--error_rate",nargs='+',default=[0.0],type=float,help="the percentage base error rate (default value is 0.0)")
    parser.add_argument("--write_fasta", action="store_true", help="whether to write out the FASTA file for the variant")
    parser.add_argument("--debug", action="store_true", help="whether to write out helpful statements to STDOUT")
    options = parser.parse_args()

    # load in the covid reference using gumpy
    reference=gumpy.Genome(options.reference, show_progress_bar=True)

    if options.debug:
        print('read reference!')

    error_rates=numpy.array(options.error_rate)/100.

    bases = {'a','t','c','g'}

    assert options.tech in ['illumina', 'nanopore']

    for snps in options.snps:

        if options.debug:
            print("copying reference...")

        current_variant = copy.deepcopy(reference)

        if options.debug:
            print("done!")

        if snps > 0:

            if options.debug:
                print("introducing SNPs...")

            snp_indices = numpy.random.choice(current_variant.nucleotide_index,size=snps,replace=False)

            for i in snp_indices:
                mask = current_variant.nucleotide_index == i
                current_base = current_variant.nucleotide_sequence[mask]

                # don't mutate an N
                if current_base == 'n':
                    continue

                new_bases = bases ^ set(current_base)
                new_base = random.choice(list(new_bases))
                current_variant.nucleotide_sequence[mask] = new_base

                print(i, current_base[0], new_base)
            
            if options.variant_file:

                assert pathlib.Path(options.variant_file).is_file(), "file does not exist!"

                with open(options.variant_file) as handle:
                    offset = 0 #Offset to account for indels
                    for line in handle:

                        cols = line.rstrip().split(' ')
                        pos = int(cols[0]) + offset
                        if cols[1] == 'ins':
                            #Ins
                            bases = list(cols[2])
                            #Insert the bases
                            current_variant.nucleotide_sequence = numpy.insert(current_variant.nucleotide_sequence, pos+1, bases)
                            #Add to the offset to allow adjusted pos
                            offset += len(bases)
                        if cols[1] == 'del':
                            #Del
                            #Del specified by `<pos> del <n bases>`
                            bases = int(cols[2])
                            current_variant.nucleotide_sequence = numpy.delete(current_variant.nucleotide_sequence, range(pos+1, pos+bases+1))
                            #Subtract from offset to allow adjusted pos
                            offset -= bases
                        else:
                            #SNPs
                            mask = current_variant.nucleotide_index == pos
                            current_base = current_variant.nucleotide_sequence[mask]
                            assert current_base == cols[1], cols[1]

                            current_variant.nucleotide_sequence[mask] = cols[2]

                        print(cols[0], cols[1], cols[2])

            if options.debug:
                print("done!")

        if options.debug:
            print("making input genome")

        input_genome = ''.join(i for i in current_variant.nucleotide_sequence)
        input_genome = pyfastaq.sequences.Fasta(id_in=current_variant.name,
                                                 seq_in=input_genome.upper())
        input_genome = input_genome.to_Fastq([40] * len(input_genome))

        if options.debug:
            print("done!")

        depth = options.depth

        for error_rate in error_rates:

            assert 100>error_rate>=0

            for repeat in range(options.repeats):

                # build an automatic output string if no output is specified or it is a pure directory
                if options.output is None or pathlib.Path(options.output).is_dir():
                    outputstem = reference.name + '-' + options.tech + '-' + str(snps) + 'snps-' + str(depth) + 'd-' +\
                                    str(error_rate) + 'e-'
                    outputstem += str(repeat)
                else:
                    outputstem = options.output

                # if options.write_fasta:
                #
                #     current_variant.save_fasta(outputstem+".fasta",\
                #                                 fixed_length=False,\
                #                                 overwrite_existing=True,\
                #                                 description=current_variant.name)

                coverage = numpy.zeros(reference.length)
                read_number = 0
                prev_coverage = 0

                if options.tech=='illumina':

                    with open(outputstem+"_1.fastq", "w") as f1, open(outputstem+"_2.fastq", "w") as f2:

                        while numpy.mean(coverage)<depth:

                            current_coverage = int(numpy.mean(coverage) / 1) * 1

                            if  current_coverage > prev_coverage:
                                prev_coverage = current_coverage

                            read_length =  int(numpy.random.normal(options.read_length, options.read_stddev))

                            read_start = numpy.random.choice(reference.nucleotide_index)

                            if read_start + read_length > reference.length:
                                continue

                            for i in range(read_length):
                                coverage[read_start+i]+=2

                            read1 = input_genome.subseq(read_start, read_start + read_length)
                            read2 = input_genome.subseq(read_start, read_start + read_length)
                            read2.revcomp()

                            if error_rate > 0:
                                read1.seq = gpas_testing.mutate_read(read1.seq, error_rate=error_rate)
                                read2.seq = gpas_testing.mutate_read(read2.seq, error_rate=error_rate)

                            read1.id = f"{read_number}.{i} /1"
                            read2.id = f"{read_number}.{i} /2"

                            print(read1, file=f1)
                            print(read2, file=f2)

                            read1.id = f"{read_number}.{i}.2 /2"
                            read2.id = f"{read_number}.{i}.2 /1"
                            print(read1, file=f2)
                            print(read2, file=f1)

                            read_number+=1

                elif options.tech == 'nanopore':

                    with open(outputstem+".fastq", "w") as f1:

                        while numpy.mean(coverage)<depth:

                            current_coverage = int(numpy.mean(coverage) / 10) * 10

                            if  current_coverage > prev_coverage:
                                prev_coverage = current_coverage

                            read_length =  int(numpy.random.normal(options.read_length, options.read_stddev))

                            read_start = numpy.random.choice(reference.nucleotide_index)

                            if read_start + read_length > reference.length:
                                continue

                            for i in range(read_length):
                                coverage[read_start+i]+=1

                            read1 = input_genome.subseq(read_start, read_start + read_length)

                            if error_rate > 0:
                                read1.seq = gpas_testing.mutate_read(read1.seq, error_rate=error_rate)

                            if numpy.random.choice([True, False]):
                                read1.id = f"{read_number}.{i} forward"
                                print(read1, file=f1)
                            else:
                                read1.id = f"{read_number}.{i}.2 reverse"
                                read1.revcomp()
                                print(read1, file=f1)

                            read_number+=1

                else:
                    raise ValueError("--tech must be one of illumina or nanopore!")

                print(read_number, numpy.min(coverage), numpy.max(coverage), numpy.sum(coverage>=10))
