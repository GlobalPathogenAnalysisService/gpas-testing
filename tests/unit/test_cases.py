import gumpy, copy, numpy, pkg_resources

import gpas_covid_synthetic_reads as gcsr

covid_reference=gumpy.Genome(pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/MN908947.3.gbk'))

def test_variant_definitions():

    variant_definitions=gcsr.load_variant_definitions('../variant_definitions')

    assert 'alpha' in variant_definitions.keys()
    assert 'delta' in variant_definitions.keys()
    assert 'omicron' in variant_definitions.keys()
    assert 'B.1.1.7' not in variant_definitions.keys()
    assert 'Alpha' not in variant_definitions.keys()

def test_mutate_read():

    # assert gcsr.mutate_read('AAAAA',0)=='AAAAA'

    assert gcsr.mutate_read('AAAAA',debug_mutations={0:'T'})=='TAAAA'

    assert gcsr.mutate_read('AAAAA',debug_mutations={0:'T',4:'C'})=='TAAAC'

    assert gcsr.mutate_read('AAAAA',debug_mutations={1:'A'})=='AAAAA'

def test_gamma():

    # Chose Gamma because it has one insertion and one deletion

    variant_definitions=gcsr.load_variant_definitions('variant_definitions')

    gamma=gcsr.VariantGenome(covid_reference, variant_definitions, 'gamma')

    start=11280
    stop=11305

    mask=(gamma.reference.nucleotide_index>start) & (gamma.reference.nucleotide_index<stop)
    test_string=''.join(i for i in gamma.reference.nucleotide_sequence[mask])
    assert test_string=='tagtttgtctggttttaagctaaa'

    test_indices=gamma.index_lookup[(gamma.index_lookup>start) & (gamma.index_lookup<stop)]
    mask=numpy.isin(gamma.reference.nucleotide_index,test_indices)
    test_string=''.join(i for i in gamma.reference.nucleotide_sequence[mask])
    assert test_string=='tagtttgaagctaaa'

    start=28258
    stop=28268
    mask=(gamma.reference.nucleotide_index>start) & (gamma.reference.nucleotide_index<stop)
    test_string=''.join(i for i in gamma.reference.nucleotide_sequence[mask])
    assert test_string=='aacgaacaa'

    gamma_genome = gamma.expected.build_genome_string()

    assert len(gamma_genome)==len(gamma.index_lookup)
    i0=numpy.where(gamma.index_lookup==start)[0][0]+1
    i1=numpy.where(gamma.index_lookup==stop)[0][0]
    assert gamma_genome[i0:i1]=='aacgaacaaacaa'


# for i in c*csv; do awk -F ',' '(NR>1) {print FILENAME, $5}' $i; done
