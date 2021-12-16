import gumpy, copy, numpy

import gpas_covid_perfect_reads as gcpr

covid_reference=gumpy.Genome('config/MN908947.3.gbk')

def test_variant_definitions():

    variant_definitions=gcpr.load_variant_definitions('../variant_definitions')

    assert 'Alpha' in variant_definitions.keys()
    assert 'Delta' in variant_definitions.keys()
    assert 'Omicron' in variant_definitions.keys()
    assert 'B.1.1.7' not in variant_definitions.keys()
    assert 'alpha' not in variant_definitions.keys()

def test_mutate_read():

    assert gcpr.mutate_read('AAAAA',0)=='AAAAA'

    assert gcpr.mutate_read('AAAAA',debug_mutations={0:'T'})=='TAAAA'

    assert gcpr.mutate_read('AAAAA',debug_mutations={0:'T',4:'C'})=='TAAAC'

    assert gcpr.mutate_read('AAAAA',debug_mutations={1:'A'})=='AAAAA'

def test_gamma():

    # Chose Gamma because it has one insertion and one deletion

    variant_definitions=gcpr.load_variant_definitions('../variant_definitions')

    gamma=gcpr.VariantGenome(covid_reference, variant_definitions['Gamma'])

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

    assert len(gamma.genome)==len(gamma.index_lookup)
    i0=numpy.where(gamma.index_lookup==start)[0][0]+1
    i1=numpy.where(gamma.index_lookup==stop)[0][0]
    assert gamma.genome[i0:i1]=='aacgaacaaacaa'
