import copy
import numpy

import gpas_covid_synthetic_reads as gcsr

class PangoGenome(object):

    def __init__(self, reference_genome, pango_definition, name, indels=None):

        self.name = name
        if self.name!='reference':
            self.definition = pango_definition[self.name]
        self.reference = reference_genome
        self.expected = copy.deepcopy(reference_genome)
        self.index_lookup=copy.deepcopy(reference_genome.nucleotide_index)
        if indels is None:
            self.indels = []
        else:
            self.indels = indels

        self.amino_acid_to_codon = gcsr.create_amino_acid_to_codon(self.expected)

        if self.name != 'reference':
            self._create_variant()
            assert self.expected != self.reference

        self.input1 = copy.deepcopy(self.expected)
        self.input2 = copy.deepcopy(self.expected)


    def _create_variant(self):

        constellation_to_gumpy_lookup = {'spike': 'S',
                                 'S': 'S',
                                 's' : 'S',
                                 'm': 'M',
                                 'M': 'M',
                                 'e': 'E',
                                 'E': 'E',
                                 'n': 'N',
                                 'N': 'N',
                                 'orf1ab': 'orf1ab',
                                 '1ab': 'orf1ab',
                                 'ORF1ab' : 'orf1ab',
                                 'ORF3a': 'ORF3a',
                                 'orf3a': 'ORF3a',
                                 'ORF7a': 'ORF7a',
                                 'orf6': 'ORF6',
                                 'ORF8': 'ORF8',
                                 '8': 'ORF8'}

        sample = {}

        all_mutations_identified = True

        for i in self.definition['sites']:

            cols = i.split(':')

            if cols[0] == 'nuc':

                # deal with nasty insertions
                if '+' in cols[1]:

                    ref = cols[1][0]
                    pos = cols[1].find('+')
                    idx = int(cols[1][0:pos])
                    alt = cols[1][pos+1:].lower()

                    mask = self.expected.nucleotide_index == idx

                    self.expected.is_indel[mask] = True
                    self.expected.indel_length[mask] = 1 * len(alt)
                    self.expected.indel_nucleotides[mask] = alt
                    self.indels.append((idx, 1*len(alt)))

                else:

                    # split into the ref, alt bases and the genome index
                    ref = cols[1][0]
                    idx = int(cols[1][1:-1])
                    alt = cols[1][-1]

                    mask = self.expected.nucleotide_index == idx

                    # have to bypass since this appears twice in the BA.2 definition...
                    if i !='nuc:C15714T':
                        # insist that the specified ref/before base matches what is the reference genome
                        assert self.expected.nucleotide_sequence[mask][0] == ref.lower(), (self.expected.nucleotide_sequence[mask][0], cols)

                    # now mutate the base
                    self.expected.nucleotide_sequence[mask] = alt.lower()

            elif cols[0] == 'del':

                idx = int(cols[1])

                number_bases_deleted = int(cols[2])

                mask=self.expected.nucleotide_index == idx

                mask2 = numpy.isin(self.expected.nucleotide_index, [i for i in range(idx, idx + number_bases_deleted)])
                deleted_nucs = self.expected.nucleotide_sequence[mask2]

                bases_deleted = ''.join(i for i in deleted_nucs)

                self.expected.is_indel[mask] = True
                self.expected.indel_nucleotides[mask] = bases_deleted
                self.expected.indel_length[mask] = -1 * number_bases_deleted
                self.indels.append((idx, -1 * number_bases_deleted))

            elif cols[0] in constellation_to_gumpy_lookup.keys() and constellation_to_gumpy_lookup[cols[0]] in self.expected.genes.keys():

                # translate to what the gene is called in gumpy (and thence the GenBank file)
                gene_name = constellation_to_gumpy_lookup[cols[0]]

                assert self.expected.contains_gene(gene_name)

                # build the Gene object
                gene = self.expected.build_gene(gene_name)

                number_of_amino_acids=sum([i.isalpha() for i in cols[1]])

                if cols[1][-1]=='*':
                    number_of_amino_acids+=1

                if cols[1][-1] == '-':

                    aa_num = int(cols[1][number_of_amino_acids:-1])

                    for i in range(number_of_amino_acids):

                        ref = cols[1][i]

                        mask = gene.amino_acid_number == aa_num + i

                        assert gene.amino_acid_sequence[mask] == ref

                    # find out the genome nucleotide indices corresponding to the start of this codon
                    idx=gene.nucleotide_index[gene.gene_position==aa_num][0]

                    # create a genome mask
                    mask=numpy.isin(self.expected.nucleotide_index,idx)
                    self.expected.is_indel[mask] = True
                    self.expected.indel_length[mask] = -1 * (3 * number_of_amino_acids)
                    self.indels.append((idx, -1 * (3 * number_of_amino_acids)))

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
                        possible_codons = self.amino_acid_to_codon[alt]

                        # work out which of these requires the fewest number of base changes
                        new_codon = gcsr.determine_closet_codon(current_codon, possible_codons)

                        # find out the genome nucleotide indices corresponding to this codon
                        idx=gene.nucleotide_index[gene.gene_position==aa_num]

                        # create a genome mask
                        mask=numpy.isin(self.expected.nucleotide_index,idx)

                        # mutate the genome to the new codon which corresponds to the mutated amino acid
                        self.expected.nucleotide_sequence[mask] = numpy.array([i for i in new_codon])

        #     elif cols[0].upper() in constellation_genome['genes'].keys():
        #         print('mutation in gene, harder')
            # elif cols[0] in name_to_gene_lookup.keys():
            #     print(lineage,'mutation in gene, but called after name, harder', cols)
            # elif cols[0].upper()=='ORF1AB':
            #     print(lineage,'mutation in ORF1ab', cols)
            else:
                all_mutations_identified = False
                # print(lineage,i)

        # if all_mutations_identified:
        #     print(self.name)
        # self.expected.save_fasta(self.name+".fasta")
