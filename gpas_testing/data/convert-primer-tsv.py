#!/usr/bin/env python3

import pandas

import gumpy
import gpas_testing

covid_reference = gumpy.Genome('MN908947.3.gbk')

for primer_name in ['covid-artic-v3', 'covid-artic-v4', 'covid-midnight-1200', 'covid-ampliseq-v1']:

    primer_scheme_file = primer_name + '.vwf.tsv'

    # load the definitions of the amplicons
    primers = pandas.read_csv(primer_scheme_file,
                              sep='\t',
                              dtype={'start': int},
                              skiprows=1,
                              names=['pool', 'name', 'left_or_right', 'seq', 'start'])

    def assign_end(row):
        return(row['start']+len(row.seq)-1)

    primers['end'] = primers.apply(assign_end, axis=1)

    def find_handedness(row):
        return(row['name'].split('_')[-1])

    primers['hand'] = primers.apply(find_handedness, axis=1)

    def number_primer(row):
        if 'pool' in row['pool']:
            return int(row['pool'].split('_pool')[0].split('_')[-1])
        else:
            return int(row['pool'].split('_')[-1])

    # primers = primers[['pool', 'number']]

    primers['number'] = primers.apply(number_primer, axis=1)

    primers.to_csv(primer_name + '.primers.csv', index=False)

    left = primers.loc[primers.hand == 'LEFT']
    right = primers.loc[primers.hand == 'RIGHT']

    amplicons = left.merge(right, left_on='pool', right_on='pool', suffixes=['_left', '_right'])
    amplicons = amplicons[['pool', 'number_left', 'start_left', 'end_left', 'start_right', 'end_right']]
    amplicons.rename(columns={'pool': 'name', 'number_left': 'number'}, inplace=True)
    amplicons['start'] = amplicons['start_left']
    amplicons['end'] = amplicons['end_right']
    amplicons['length'] = amplicons['end'] - amplicons['start']

    amplicons[['start_amplicon', 'end_amplicon']] = amplicons.apply(gpas_testing.define_amplicon,
                                                                    amplicons=amplicons,
                                                                    reference_genome=covid_reference,
                                                                    axis=1)

    amplicons.to_csv(primer_name + '.amplicons.csv', index=False)
