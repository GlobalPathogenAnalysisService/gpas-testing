#! /usr/bin/env python

import pandas
import glob
import os
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", required=True, help="the batch uuid that identifies the folder on sp3:/work/output where the results are stored")
    parser.add_argument("--tech", required=True, help="either illumina or nanopore")
    parser.add_argument("--outfile", default='results.csv.gz', help="the name of the output CSV file")
    options = parser.parse_args()

    # this CSV let's us map from local UID to GPAS UID
    sample_names = pandas.read_csv('sample_names.csv', names=['in_fasta','out_fasta'])

    # as long as the filenames only contain a single underscore, we can find the stem which contains metatdata when testing
    def find_stem(row):
        return(row['in_fasta'].split('_')[0])

    sample_names['filestem'] = sample_names.apply(find_stem, axis=1)

    def find_characteristics(row):

        cols = row['filestem'].split('-')

        snps = int(cols[4].split('snps')[0])
        depth = int(cols[5].split('d')[0])
        error = float(cols[6].split('e')[0])

        if cols[7][0] == 'a':
            dropped_amplicons = cols[7]
            repeat = int(cols[8].split('.fasta')[0])
        else:
            dropped_amplicons = None
            repeat=int(cols[7].split('.fasta')[0])

        return pandas.Series([cols[0], cols[1], cols[2], cols[3], snps, depth, error, dropped_amplicons, repeat])

    sample_names[['tech', 'primer_scheme', 'lineage', 'lineage_definition', 'snps', 'depth', 'error', 'dropped_amplicons', 'repeat']] = sample_names.apply(find_characteristics,axis=1)

    def read_fasta(filename):

        if os.path.exists(filename):
            INPUT = open(filename, 'r')
            INPUT.readline()
            genome=''
            for i in INPUT:
                genome += i.rstrip()
        else:
            genome = ''

        return(genome)

    def compare_genomes(row):

        in_genome = read_fasta(row['filestem'] + '.fasta')
        out_genome = read_fasta(options.batch + '/consensus_seqs/' + row['out_fasta'] + '.fasta')

        # only mark as correctly assembled if the length of the output string is >29000 and it is found in the input string
        match = out_genome in in_genome
        length = True if len(out_genome) > 29000 else False
        result = match and length

        return(pandas.Series([result, len(in_genome), len(out_genome), in_genome, out_genome]))

    sample_names[['assembled_correct', 'len_in', 'len_out', 'in_genome', 'out_genome']]=sample_names.apply(compare_genomes,axis=1)



    print("Assembled correctly:\n", sample_names.assembled_correct.value_counts())

    def read_pango(row):

        pangofile = options.batch + '/analysis/pango/' + options.tech + '/' + row['out_fasta'] + '_lineage_report.csv'

        if os.path.exists(pangofile):
            df = pandas.read_csv(pangofile)
            return pandas.Series([df['lineage'].values[0], df['conflict'].values[0], df['ambiguity_score'].values[0],
                                  df['scorpio_call'].values[0], df['scorpio_support'].values[0],
                                  df['scorpio_conflict'].values[0]])
        else:

            return pandas.Series([None, None, None, None, None, None])

    sample_names[['pango_lineage', 'pango_conflict', 'pango_ambiguity_score', 'scorpio_call', 'scorpio_support', 'scorpio_conflict']] = sample_names.apply(read_pango, axis=1)

    sample_names = sample_names[['in_fasta', 'out_fasta', 'filestem', 'tech', 'primer_scheme', 'lineage',
                                'lineage_definition', 'snps', 'depth', 'error', 'dropped_amplicons',
                                'repeat', 'assembled_correct', 'len_in', 'len_out',
                                'pango_lineage', 'pango_conflict',
                                'pango_ambiguity_score', 'scorpio_call', 'scorpio_support',
                                'scorpio_conflict', 'in_genome', 'out_genome']]
    
    sample_names.to_csv(options.outfile, index=False)
            
