import gumpy, copy, numpy, pkg_resources

import gpas_testing

covid_reference=gumpy.Genome(pkg_resources.resource_filename("gpas_testing", 'data/MN908947.3.gbk'))

def test_pango_definitions():

    pango_definitions=gpas_testing.load_pango_definitions('constellations')

    assert 'cBA.1' in pango_definitions.keys()
    assert 'B.1.1.7' not in pango_definitions.keys()

def test_mutate_read():

    # assert gpas_testing.mutate_read('AAAAA',0)=='AAAAA'

    assert gpas_testing.mutate_read('AAAAA',debug_mutations={0:'T'})=='TAAAA'

    assert gpas_testing.mutate_read('AAAAA',debug_mutations={0:'T',4:'C'})=='TAAAC'

    assert gpas_testing.mutate_read('AAAAA',debug_mutations={1:'A'})=='AAAAA'

def test_gamma():

    # Chose Gamma because it has one insertion and one deletion
    pango_definitions=gpas_testing.load_pango_definitions('constellations')

    gamma=gpas_testing.PangoGenome(covid_reference, pango_definitions, 'cP.1')

    test_string = gamma.reference.build_genome_string(nucleotide_index_range=(11280,11305))
    assert test_string == 'ctagtttgtctggttttaagctaaa'

    test_string = gamma.expected.build_genome_string(nucleotide_index_range=(11280,11305))
    assert test_string == 'ctagttttaagctaaa'

    gamma_genome = gamma.expected.build_genome_string()

    assert len(gamma_genome) == 29903-9+4
