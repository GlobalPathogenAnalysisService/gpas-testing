import pkg_resources
import pathlib
import random
import numpy
import pandas
import json
import yaml
from collections import defaultdict

def define_amplicon(tmp, amplicons, reference_genome):

    chosen_amplicon = tmp['name']
    row = amplicons[amplicons.name == chosen_amplicon]
    # PWF: this used to be >= but end_left is the index of the last base in the left primer
    mask = (reference_genome.nucleotide_index > int(row['end_left'])) & (reference_genome.nucleotide_index < int(row['start_right']))

    for idx, row in amplicons.iterrows():

        if row['name'] == chosen_amplicon:
            continue

        # PWF: this used to be >= but end_left is the index of the last base in the left primer
        current_amplicon_mask = (reference_genome.nucleotide_index > int(row['end_left'])) & (reference_genome.nucleotide_index < int(row['start_right']))

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
                variant_definitions[a['who-label'].lower()]=a
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


def load_lineages_dataframe():

    cov_lineages = pkg_resources.resource_filename('gpas_covid_synthetic_reads','data/cov-lineages.csv')

    lineages_reference = pandas.read_csv(cov_lineages)

    return lineages_reference


def load_pango_definitions(path):

    lineages_reference = load_lineages_dataframe()

    pango_definitions = {}

    constellations_path = pathlib.Path(path) / "constellations/definitions"

    for i in lineages_reference['pango_lineage']:
        who_label = lineages_reference[lineages_reference['pango_lineage']==i]['who_label'].values[0].lower()

        if who_label == 'epsilon' and i == 'cB.1.427':
            continue
        elif who_label == 'omicron':
            if i == 'cB.1.1.529':
                who_label = 'omicronB.1.1.529'
            elif i == 'cBA.1':
                who_label = 'omicronBA.1'
            elif i == 'cBA.2':
                who_label = 'omicronBA.2'
            else:
                continue

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
