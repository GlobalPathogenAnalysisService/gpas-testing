import gumpy, copy, numpy, pkg_resources

import gpas_covid_synthetic_reads as gcsr

covid_reference=gumpy.Genome(pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/MN908947.3.gbk'))

def test_pango_definitions():

    pango_definitions=gcsr.load_pango_definitions('constellations')

    assert 'cBA.1' in pango_definitions.keys()
    assert 'B.1.1.7' not in pango_definitions.keys()

def test_mutate_read():

    # assert gcsr.mutate_read('AAAAA',0)=='AAAAA'

    assert gcsr.mutate_read('AAAAA',debug_mutations={0:'T'})=='TAAAA'

    assert gcsr.mutate_read('AAAAA',debug_mutations={0:'T',4:'C'})=='TAAAC'

    assert gcsr.mutate_read('AAAAA',debug_mutations={1:'A'})=='AAAAA'

def test_gamma():

    # Chose Gamma because it has one insertion and one deletion
    pango_definitions=gcsr.load_pango_definitions('constellations')

    gamma=gcsr.PangoGenome(covid_reference, pango_definitions, 'cP.1')

    test_string = gamma.reference.build_genome_string(nucleotide_index_range=(11280,11305))
    assert test_string == 'ctagtttgtctggttttaagctaaa'

    test_string = gamma.expected.build_genome_string(nucleotide_index_range=(11280,11305))
    assert test_string == 'ctagttttaagctaaa'

    gamma_genome = gamma.expected.build_genome_string()

    assert len(gamma_genome) == 29903-9+4
