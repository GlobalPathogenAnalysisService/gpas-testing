import pathlib, copy, numpy, random

import yaml
import pandas

def define_amplicon(tmp, amplicons, reference_genome):

    chosen_amplicon = tmp['name']
    row = amplicons[amplicons.name == chosen_amplicon]
    mask = (reference_genome.nucleotide_index >= int(row['end_left'])) & (reference_genome.nucleotide_index < int(row['start_right']))

    for idx, row in amplicons.iterrows():

        if row['name'] == chosen_amplicon:
            continue

        current_amplicon_mask = (reference_genome.nucleotide_index >= int(row['start'])) & (reference_genome.nucleotide_index < int(row['end']))

        overlap_region = current_amplicon_mask & mask

        overlap_region = numpy.logical_not(overlap_region)

        mask = overlap_region & mask

    return(pandas.Series( [reference_genome.nucleotide_index[mask].min(), reference_genome.nucleotide_index[mask].max()]))


def load_variant_definitions(path):

    variant_definitions_path=pathlib.Path(path) / "variant_yaml/"
    variant_definitions_files=variant_definitions_path.glob('**/*.yml')

    variant_definitions={}
    for i in variant_definitions_files:
        with open(i) as INPUT:
            a=yaml.safe_load(INPUT)
            if 'who-label' in a.keys() and a['who-label'] is not None:
                variant_definitions[a['who-label']]=a
    return(variant_definitions)

def mutate_read(read,error_rate=0,snps=0,debug_mutations=None):

    assert error_rate<=1

    if debug_mutations is not None:
        mask=numpy.isin(numpy.arange(0,len(read)),list(debug_mutations.keys()))
    elif snps>0:
        positions=random.sample(range(len(read)),snps)
        mask=numpy.isin(numpy.arange(0,len(read)),positions)
    elif error_rate>0:
        mask=numpy.random.uniform(size=len(read))<error_rate
    else:
        raise ValueError('Read will not be mutated!')

    # only if there are more than zero mutations
    if numpy.sum(mask)==0:

        return(read)

    else:

        bases={'A','T','C','G'}

        # create a list of new mutations
        new_bases=[]

        # convert the read into an array of chars
        r=numpy.array(list(read))

        if debug_mutations is None:

            for i in r[mask]:
                new_bases.append(random.sample(bases ^ set(i),1)[0])

        else:

            new_bases=list(debug_mutations.values())
            print(new_bases,list(debug_mutations.values()))

        # set the new values
        r[mask]=numpy.array(new_bases)

        # recreate the sequence
        sequence=''.join(i for i in r)

        return(sequence)

class VariantGenome(object):

    def __init__(self, reference_genome, variant_definition):

        self.definition=variant_definition
        self.name = variant_definition['who-label']
        self.reference=reference_genome
        self.variant=copy.deepcopy(reference_genome)
        self.index_lookup=copy.deepcopy(reference_genome.nucleotide_index)

        self._create_variant()

        assert self.variant != self.reference

        self.genome=self.variant.build_genome_string()

        # self.variant.save_fasta(self.name+".fasta",\
        #                         fixed_length=False,\
        #                         overwrite_existing=True,\
        #                         description=self.name)
        #
    def _create_variant(self):

        for mutation in self.definition['variants']:

            mutation_type = mutation['type']

            if mutation_type=='SNP':

                # build the mask to isolate the nucleotide
                mask=self.variant.nucleotide_index==mutation['one-based-reference-position']

                # check the definition agrees with the reference sequence
                assert self.variant.nucleotide_sequence[mask]==mutation['reference-base'].lower()

                # only now mutate the base
                self.variant.nucleotide_sequence[mask]=mutation['variant-base'].lower()

            elif mutation_type=='MNP':

                length=len(mutation['reference-base'])

                for i in range(length):

                    mask=self.variant.nucleotide_index==(mutation['one-based-reference-position']+i)

                    # check the definition agrees with the reference sequence
                    assert self.variant.nucleotide_sequence[mask][0]==mutation['reference-base'].lower()[i]

                    # only now mutate the base
                    self.variant.nucleotide_sequence[mask]=mutation['variant-base'].lower()[i]

            elif mutation_type=='insertion':

                indel_length=(len(mutation['variant-base'])-1)
                indel_index=mutation['one-based-reference-position']+1

                mask=self.variant.nucleotide_index==indel_index

                self.variant.is_indel[mask]=True
                self.variant.indel_length[mask]=indel_length
                self.variant.indel_nucleotides[mask]=mutation['variant-base'][1:].lower()

                mask1=self.index_lookup[(self.index_lookup<=indel_index)]
                mask2=numpy.repeat(indel_index,indel_length)
                mask3=self.index_lookup[self.index_lookup>indel_index]
                self.index_lookup=numpy.concatenate([mask1,mask2,mask3])

            elif mutation_type=='deletion':

                indel_length=-1*(len(mutation['reference-base'])-1)
                indel_index=mutation['one-based-reference-position']+1

                mask=self.variant.nucleotide_index==indel_index
                self.variant.is_indel[mask]=True
                self.variant.indel_length[mask]=indel_length
                self.variant.indel_nucleotides[mask]=mutation['reference-base'][1:].lower()

                mask2=(self.index_lookup<indel_index) | (self.index_lookup>=(indel_index-indel_length))
                self.index_lookup=self.index_lookup[mask2]
