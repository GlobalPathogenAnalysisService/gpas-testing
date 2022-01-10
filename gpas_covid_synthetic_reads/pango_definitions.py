import pathlib
import copy
import random
import json
from collections import defaultdict
import pkg_resources

import numpy
import pandas


def load_lineages_dataframe():

    cov_lineages = pkg_resources.resource_filename('gpas_covid_synthetic_reads','data/cov-lineages.csv')

    lineages_reference = pandas.read_csv(cov_lineages)

    return lineages_reference

def load_pango_definitions(path, lineages_reference):

    pango_definitions = {}

    constellations_path = pathlib.Path(path) / "constellations/definitions"

    for i in lineages_reference['pango_lineage']:

        with open(constellations_path / (i + '.json') ) as INPUT:

            pango_definitions[i]=json.load(INPUT)

    return pango_definitions

def create_amino_acid_to_codon(genome):

    spike = genome.build_gene('S')

    amino_acid_to_codon=defaultdict(list)

    for codon in spike.codon_to_amino_acid:

        if 'x' not in codon and 'z' not in codon and 'o' not in codon:

            amino_acid=spike.codon_to_amino_acid[codon]

            amino_acid_to_codon[amino_acid].append(codon)

    return amino_acid_to_codon

def determine_closet_codon(current_codon, possible_codons):

    max_bases=4

    for putative_codon in possible_codons:

        bases_to_change=sum([i!=j for i,j in zip(current_codon,putative_codon)])

        if bases_to_change<max_bases:

            new_codon = putative_codon

            max_bases=bases_to_change

    return(new_codon)


class PangoGenome(object):

    def __init__(self, sample_genome, pango_definition, lineage):

        amino_acid_to_codon=create_amino_acid_to_codon(sample_genome)

        constellation_to_gumpy_lookup = {'spike': 'S',
                                 'S': 'S',
                                 's' : 'S',
                                 'm': 'M',
                                 'M': 'M',
                                 'e': 'E',
                                 'E': 'E',
                                 'n': 'N',
                                 'N': 'N',
                                 'orf1ab': 'ORF1ab',
                                 '1ab': 'ORF1ab',
                                 'ORF1ab' : 'ORF1ab',
                                 'ORF3a': 'ORF3a',
                                 'orf3a': 'ORF3a',
                                 'ORF7a': 'ORF7a',
                                 'orf6': 'ORF6',
                                 'ORF8': 'ORF8',
                                 '8': 'ORF8'}

        sample = {}

        all_mutations_identified = True

        for i in pango_definition['sites']:

            cols = i.split(':')

            if cols[0] == 'nuc':

                # deal with nasty insertions
                if '+' in cols[1]:

                    ref = cols[1][0]
                    pos = cols[1].find('+')
                    idx = int(cols[1][1:pos])
                    alt = cols[1][pos+1:].lower()

                    mask = sample_genome.nucleotide_index == idx

                    sample_genome.is_indel[mask] = True
                    sample_genome.indel_length[mask] = 1 * len(alt)
                    sample_genome.indel_nucleotides[mask] = alt

                else:

                    # split into the ref, alt bases and the genome index
                    ref = cols[1][0]
                    idx = int(cols[1][1:-1])
                    alt = cols[1][-1]

                    mask = sample_genome.nucleotide_index == idx

                    # have to bypass since this appears twice in the BA.2 definition...
                    if i !='nuc:C15714T':
                        # insist that the specified ref/before base matches what is the reference genome
                        assert sample_genome.nucleotide_sequence[mask][0] == ref.lower(), (sample_genome.nucleotide_sequence[mask][0], cols)

                    # now mutate the base
                    sample_genome.nucleotide_sequence[mask] = alt.lower()

            elif cols[0] == 'del':

                idx = int(cols[1])

                number_bases_deleted = int(cols[2])

                mask=sample_genome.nucleotide_index == idx

                sample_genome.is_indel[mask] = True
                sample_genome.indel_nucleotides[mask] = -1 * number_bases_deleted

            elif cols[0] in constellation_to_gumpy_lookup.keys() and constellation_to_gumpy_lookup[cols[0]] in sample_genome.genes.keys():

                # translate to what the gene is called in gumpy (and thence the GenBank file)
                gene_name = constellation_to_gumpy_lookup[cols[0]]

                assert sample_genome.contains_gene(gene_name)

                # build the Gene object
                gene = sample_genome.build_gene(gene_name)

                number_of_amino_acids=sum([i.isalpha() for i in cols[1]])

                if cols[1][-1]=='*':
                    number_of_amino_acids+=1

                if cols[1][-1] == '-':

                    aa_num = int(cols[1][number_of_amino_acids:-1])

                    for i in range(number_of_amino_acids):

                        ref = cols[1][i]

                        mask = gene.amino_acid_number == aa_num + i

                        assert gene.amino_acid_sequence[mask] == ref

                    # find out the genome nucleotide indices corresponding to this codon
                    idx=gene.nucleotide_index[gene.gene_position==aa_num]

                    # create a genome mask
                    mask=numpy.isin(sample_genome.nucleotide_index,idx)

                    sample_genome.is_indel[mask] = True
                    sample_genome.indel_nucleotides[mask] = -1 * (3 * number_of_amino_acids)

                else:

                    # this must be an even number
                    assert number_of_amino_acids % 2 == 0, (cols, number_of_amino_acids)

                    # split into the ref, alt amino acids and the amino acid number
                    if number_of_amino_acids == 2:

                        refs = [cols[1][0]]
                        aa_nums = [int(cols[1][1:-1])]
                        alts = [cols[1][-1]]
                        if cols[1][-1] == '*':
                            alts = ['!']

                    elif number_of_amino_acids == 4:

                        refs = [cols[1][0], cols[1][1]]
                        aa_nums = [int(cols[1][2:-2]), int(cols[1][2:-2])+1]
                        alts = [cols[1][-2], cols[1][-1]]

                    for ref,aa_num,alt in zip(refs,aa_nums,alts):

                        # insist that the provided REF amino acid matches that in the genome
                        mask = gene.amino_acid_number == aa_num

                        assert gene.amino_acid_sequence[mask] == ref

                        # find out the current codon
                        current_codon = gene.codons[mask][0]

                        # find out what codons would encode the ALT amino acid
                        possible_codons = amino_acid_to_codon[alt]

                        # work out which of these requires the fewest number of base changes
                        new_codon = determine_closet_codon(current_codon, possible_codons)

                        # find out the genome nucleotide indices corresponding to this codon
                        idx=gene.nucleotide_index[gene.gene_position==aa_num]

                        # create a genome mask
                        mask=numpy.isin(sample_genome.nucleotide_index,idx)

                        # mutate the genome to the new codon which corresponds to the mutated amino acid
                        sample_genome.nucleotide_sequence[mask] = numpy.array([i for i in new_codon])

        #     elif cols[0].upper() in constellation_genome['genes'].keys():
        #         print('mutation in gene, harder')
            # elif cols[0] in name_to_gene_lookup.keys():
            #     print(lineage,'mutation in gene, but called after name, harder', cols)
            # elif cols[0].upper()=='ORF1AB':
            #     print(lineage,'mutation in ORF1ab', cols)
            else:
                all_mutations_identified = False
                # print(lineage,i)

        if all_mutations_identified:
            print(lineage)
        sample_genome.save_fasta(lineage+".fasta")
