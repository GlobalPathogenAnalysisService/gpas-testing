#!/usr/bin/env python3

import copy, glob, json, yaml, argparse

import numpy
import pyfastaq
import gumpy

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_definitions",required=True,help="the path to the variant_definitions repository/folder from phe-genomics ")
    parser.add_argument("--output",required=True,help="the stem of the output file")
    parser.add_argument("--variant_name",required=False,help="a JSON file specifying the mutations to apply to the covid reference (if not specified, you'll get a wildtype sequence)")
    parser.add_argument("--reference",required=False,default='config/MN908947.3.gbk',help="the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will be used)")
    parser.add_argument("--primer_definition",default='config/covid-artic-v3.json',help="the JSON file specifying the primer scheme used (if not specified, covid-artic-v3.json will be used)")
    parser.add_argument("--read_length",default=250,type=int,help="the read length (default value is 250)")
    parser.add_argument("--depth",default=500,type=int,help="the depth (default value is 500)")
    options = parser.parse_args()

    # load in the covid reference using gumpy
    covid_reference=gumpy.Genome(options.reference)

    variant_reference=copy.deepcopy(covid_reference)

    # load the definitions of the amplicons
    with open(options.primer_definition) as f:
        amplicons = json.load(f)

    # only if a variant has been specified, otherwise output reference
    if options.variant_name:

        variant_defintion_files=glob.glob(options.variant_definitions+'/variant_yaml/*.yml')

        variant_definitions={}

        for i in variant_defintion_files:

            with open(i) as INPUT:
                a=yaml.safe_load(INPUT)
                if 'who-label' in a.keys() and a['who-label'] is not None:
                    variant_definitions[a['who-label']]=a

        # check that the variant that has been specified has a YAML file!
        assert options.variant_name in variant_definitions.keys(), "specified variant not defined here "+options.variant_definitions

        for mutation in variant_definitions[options.variant_name]['variants']:

            mutation_type = mutation['type']

            if mutation_type=='SNP':

                # build the mask to isolate the nucleotide
                mask=variant_reference.nucleotide_index==mutation['one-based-reference-position']

                # check the definition agrees with the reference sequence
                assert variant_reference.nucleotide_sequence[mask]==mutation['reference-base'].lower()

                # only now mutate the base
                variant_reference.nucleotide_sequence[mask]=mutation['variant-base'].lower()

            elif mutation_type=='MNP':

                length=len(mutation['reference-base'])

                for i in range(length):

                    mask=variant_reference.nucleotide_index==(mutation['one-based-reference-position']+i)

                    # check the definition agrees with the reference sequence
                    assert variant_reference.nucleotide_sequence[mask][0]==mutation['reference-base'].lower()[i]

                    # only now mutate the base
                    variant_reference.nucleotide_sequence[mask]=mutation['variant-base']

            elif mutation_type=='insertion':

                mask=variant_reference.nucleotide_index==mutation['one-based-reference-position']+1

                variant_reference.is_indel[mask]=True
                variant_reference.indel_length[mask]=len(mutation['variant-base'])-1
                variant_reference.indel_nucleotides[mask]=mutation['variant-base'][1:].lower()

            elif mutation_type=='deletion':

                mask=variant_reference.nucleotide_index==mutation['one-based-reference-position']+1
                variant_reference.is_indel[mask]=True
                variant_reference.indel_length[mask]=-1*(len(mutation['reference-base'])-1)
                variant_reference.indel_nucleotides[mask]=mutation['reference-base'][1:].lower()

        assert variant_reference!=covid_reference

    else:

        variant_definitions=False

    if options.variant_name:
        description = options.variant_name #+"_"+options.primer_definition.split('.')[0]+"_readlength_"+str(options.read_length)+"_depth_"+str(options.depth)
    else:
        description = "Reference"

    variant_reference.save_fasta(options.output+".fasta",\
                            fixed_length=False,\
                            overwrite_existing=True,\
                            description=description)


    # now read the FASTA file back in
    refs = {}
    pyfastaq.tasks.file_to_dict(options.output+".fasta", refs)
    assert len(refs) == 1
    ref = list(refs.values())[0]
    ref = ref.to_Fastq([40] * len(ref))


    with open(options.output+"_1.fastq", "w") as f1, open(options.output+"_2.fastq", "w") as f2:
        for amplicon_name, amplicon_d in amplicons["amplicons"].items():
            # We'll make a read pair that looks like this:
            #
            #  |------ read1 ------------->
            #                      <--------- read2 -----|
            #  |-------------- amplicon -----------------|
            assert options.read_length < amplicon_d["end"] - amplicon_d["start"] + 1 < 2 * options.read_length

            read1 = ref.subseq(amplicon_d["start"], amplicon_d["start"] + options.read_length)
            read2 = ref.subseq(amplicon_d["end"] - options.read_length, amplicon_d["end"])
            read2.revcomp()

            for i in range(0, int(options.depth / 2)):
                read1.id = f"{amplicon_name}.{i} /1"
                read2.id = f"{amplicon_name}.{i} /2"
                print(read1, file=f1)
                print(read2, file=f2)
                read1.id = f"{amplicon_name}.{i}.2 /2"
                read2.id = f"{amplicon_name}.{i}.2 /1"
                print(read2, file=f1)
                print(read1, file=f2)
