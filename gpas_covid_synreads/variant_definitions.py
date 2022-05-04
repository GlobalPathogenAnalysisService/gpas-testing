import copy
import numpy


class VariantGenome(object):

    def __init__(self, reference_genome, variant_definition, name):

        self.name = name
        assert self.name in variant_definition.keys(), 'specified lineage not found in definitions! '
        if self.name != 'reference':
            self.definition=variant_definition[self.name]
        self.reference=reference_genome
        self.expected = copy.deepcopy(reference_genome)
        self.index_lookup=copy.deepcopy(reference_genome.nucleotide_index)

        if self.name != 'reference':
            self._create_variant()
            assert self.expected != self.reference

        self.input1 = copy.deepcopy(self.expected)
        self.input2 = copy.deepcopy(self.expected)

    def _create_variant(self):

        for mutation in self.definition['variants']:

            mutation_type = mutation['type']

            if mutation_type=='SNP':

                # build the mask to isolate the nucleotide
                mask=self.expected.nucleotide_index==mutation['one-based-reference-position']

                # check the definition agrees with the reference sequence
                assert self.expected.nucleotide_sequence[mask]==mutation['reference-base'].lower()

                # only now mutate the base
                self.expected.nucleotide_sequence[mask]=mutation['variant-base'].lower()

            elif mutation_type=='MNP':

                length=len(mutation['reference-base'])

                for i in range(length):

                    mask=self.expected.nucleotide_index==(mutation['one-based-reference-position']+i)

                    # check the definition agrees with the reference sequence
                    assert self.expected.nucleotide_sequence[mask][0]==mutation['reference-base'].lower()[i]

                    # only now mutate the base
                    self.expected.nucleotide_sequence[mask]=mutation['variant-base'].lower()[i]

            elif mutation_type=='insertion':

                indel_length=(len(mutation['variant-base'])-1)
                indel_index=mutation['one-based-reference-position']+1

                mask=self.expected.nucleotide_index==indel_index

                self.expected.is_indel[mask]=True
                self.expected.indel_length[mask]=indel_length
                self.expected.indel_nucleotides[mask]=mutation['variant-base'][1:].lower()

                mask1=self.index_lookup[(self.index_lookup<=indel_index)]
                mask2=numpy.repeat(indel_index,indel_length)
                mask3=self.index_lookup[self.index_lookup>indel_index]
                self.index_lookup=numpy.concatenate([mask1,mask2,mask3])

            elif mutation_type=='deletion':

                indel_length=-1*(len(mutation['reference-base'])-1)
                indel_index=mutation['one-based-reference-position']+1

                mask=self.expected.nucleotide_index==indel_index
                self.expected.is_indel[mask]=True
                self.expected.indel_length[mask]=indel_length
                self.expected.indel_nucleotides[mask]=mutation['reference-base'][1:].lower()

                mask2=(self.index_lookup<indel_index) | (self.index_lookup>=(indel_index-indel_length))
                self.index_lookup=self.index_lookup[mask2]
