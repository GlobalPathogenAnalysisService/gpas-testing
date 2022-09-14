#!/usr/bin/env python3

import glob
import uuid
import os
import argparse
import random
import pkg_resources
from datetime import date, timedelta

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--country",default="GBR",help="the name of the country where the samples were collected")
    parser.add_argument("--tech",default='illumina',help="whether to generate illumina (paired) or nanopore (unpaired) reads")
    parser.add_argument("--file_type",default='fastq',help="whether to look for FASTQ or BAM files")
    parser.add_argument("--tag_file",default=pkg_resources.resource_filename("gpas_testing", 'data/tags.txt'),help="a plaintext file with one row per tag")
    parser.add_argument("--uuid_length",default='long',help="whether to use a long or short UUID4")
    parser.add_argument("--number_of_tags",default=1,type=int,help="how many tags to give each sample. Can be zero, or up to the number of rows in <tag_file>. Default is 1.")
    parser.add_argument("--no_rename", action="store_true", help="whether to rename the files or not.")
    parser.add_argument("--old_format", action="store_true", help="whether to use the original upload CSV format and headers")
    parser.add_argument("--organisation",default="University of Oxford",help="the name of the organisation (the user must belong to it otherwise validation will fail)")
    options = parser.parse_args()

    assert options.tech in ['illumina','nanopore']
    assert options.file_type in ['fastq', 'bam']
    assert options.uuid_length in ['long', 'short']

    tags=[]
    with open(options.tag_file,'r') as INPUT:
        for line in INPUT:
            tags.append(line.rstrip())

    assert options.number_of_tags > 0
    assert options.number_of_tags <= len(tags)
    
    if options.file_type == 'fastq':

        if options.tech == 'illumina':

            if options.old_format:
                header = 'name,fastq1,fastq2,organisation,tags,specimenOrganism,host,collectionDate,country,submissionTitle,submissionDescription,instrument_platform,instrument_model,flowcell'
            else:
                header = 'batch,run_number,sample_name,fastq1,fastq2,control,collection_date,tags,country,region,district,specimen_organism,host,instrument_platform,primer_scheme'

            if glob.glob('*_1.fastq.gz'):
                illumina_includes_r = False
                file_list = glob.glob('*_1.fastq.gz')
                file_extensions = ['_1.fastq.gz','_2.fastq.gz']
            else:
                illumina_includes_r = True
                file_list = glob.glob('*_R1.fastq.gz')
                file_extensions = ['_R1.fastq.gz','_R2.fastq.gz']
                
        elif options.tech == 'nanopore':
            if options.old_format:
                header = 'name,fastq,organisation,tags,specimenOrganism,host,collectionDate,country,submissionTitle,submissionDescription,instrument_platform,instrument_model,flowcell'
            else:
                header = 'batch,run_number,sample_name,fastq,control,collection_date,tags,country,region,district,specimen_organism,host,instrument_platform,primer_scheme'
            file_list = glob.glob('*.fastq.gz')
            file_extensions = ['.fastq.gz']
    elif options.file_type == 'bam':
        if options.old_format:
            header = 'name,bam,organisation,tags,specimenOrganism,host,collectionDate,country,submissionTitle,submissionDescription,instrument_platform,instrument_model,flowcell'
        else:
            header = 'batch,run_number,sample_name,bam,control,collection_date,tags,country,region,district,specimen_organism,host,instrument_platform,primer_scheme'
        file_list = glob.glob('*.bam')
        file_extensions = ['.bam']

    print(header)

    def build_rest_of_line_old():
        rest_of_line = options.organisation+','
        if options.number_of_tags==0:
            rest_of_line+=','
        elif options.number_of_tags==len(tags):
            for i in tags:
                rest_of_line+=i+':'
            rest_of_line=rest_of_line[:-1]+','
        else:
            for i in random.sample(tags,options.number_of_tags):
                rest_of_line+=i+':'
            rest_of_line=rest_of_line[:-1]+','
        rest_of_line+='SARS-CoV-2,human,'

        date_collected=str(date.today()-timedelta(days=random.choice(range(7))))
        rest_of_line+=date_collected+','
        rest_of_line+=options.country+','
        rest_of_line+='covid study,study of covid,'
        if options.tech=='illumina':
            rest_of_line+="Illumina,HiSeq,96"
        elif options.tech=='nanopore':
            rest_of_line+="Nanopore,GridION,96"

        return(rest_of_line)

    def build_rest_of_line():
        # no control specified at present
        rest_of_line = ','
        date_collected=str(date.today()-timedelta(days=random.choice(range(7))))
        rest_of_line+=date_collected+','
        if options.number_of_tags==0:
            rest_of_line+=','
        elif options.number_of_tags==len(tags):
            for i in tags:
                rest_of_line+=i+':'
            rest_of_line=rest_of_line[:-1]+','
        else:
            for i in random.sample(tags,options.number_of_tags):
                rest_of_line+=i+':'
            rest_of_line=rest_of_line[:-1]+','
        rest_of_line+=options.country+','
        # adminstrative domain not specified 
        rest_of_line += ','
        # district not specified 
        rest_of_line += ','
        rest_of_line+='SARS-CoV-2,human,'
        if options.tech=='illumina':
            rest_of_line+="Illumina,auto"
        elif options.tech=='nanopore':
            rest_of_line+="Nanopore,auto"

        return(rest_of_line)

    if not options.old_format:
        batch = 'B-' + str(uuid.uuid4())[-6:]
        run = str(uuid.uuid4())[-5:]

    for i in file_list:

        if options.file_type == 'fastq':
            if options.tech=='illumina':
                if illumina_includes_r:
                    filename=i.split('_R1.fastq.gz')[0]
                else:
                    filename=i.split('_1.fastq.gz')[0]
            elif options.tech=='nanopore':
                filename=i.split('.fastq.gz')[0]
        else:
            filename=i.split('.bam')[0]

        if '_' in filename and not options.no_rename:
            lineage=i.split('_')[0]
        else:
            lineage=filename

        if options.uuid_length == 'long':
            uid=str(uuid.uuid4())
        else:
            uid=str(uuid.uuid4())[-4:]

        if options.old_format:
            rest_of_line=build_rest_of_line_old()
        else:
            rest_of_line=build_rest_of_line()

        for file_extension in file_extensions:
            if not options.no_rename:
                os.rename(filename+file_extension,lineage+"_"+uid+file_extension)

        if options.file_type == 'fastq':
            if options.tech=='illumina':
                if options.old_format:
                    if illumina_includes_r:
                        line=lineage+"_"+uid+','+lineage+"_"+uid+'_R1.fastq.gz,'+lineage+"_"+uid+'_R2.fastq.gz,'+rest_of_line                        
                    else:
                        line=lineage+"_"+uid+','+lineage+"_"+uid+'_1.fastq.gz,'+lineage+"_"+uid+'_2.fastq.gz,'+rest_of_line
                elif options.no_rename:
                    if illumina_includes_r:
                        line=batch+','+run+','+lineage+','+lineage+'_R1.fastq.gz,'+lineage+'_R2.fastq.gz,'+rest_of_line
                    else:
                        line=batch+','+run+','+lineage+','+lineage+'_1.fastq.gz,'+lineage+'_2.fastq.gz,'+rest_of_line
                else:
                    if illumina_includes_r:
                        line=batch+','+run+','+lineage+'_'+uid+','+lineage+"_"+uid+'_R1.fastq.gz,'+lineage+"_"+uid+'_R2.fastq.gz,'+rest_of_line
                    else:
                        line=batch+','+run+','+lineage+'_'+uid+','+lineage+"_"+uid+'_1.fastq.gz,'+lineage+"_"+uid+'_2.fastq.gz,'+rest_of_line                        
            elif options.tech=='nanopore':
                if options.old_format:
                    line=lineage+"_"+uid+','+lineage+"_"+uid+'.fastq.gz,'+rest_of_line
                elif options.no_rename:
                    line=batch+','+run+','+lineage+','+lineage+'.fastq.gz,'+rest_of_line
                else:
                    line=batch+','+run+','+lineage+'_'+uid+','+lineage+"_"+uid+'.fastq.gz,'+rest_of_line
        else:
            if options.old_format:
                line=lineage+"_"+uid+','+lineage+"_"+uid+'.bam,'+rest_of_line
            elif options.no_rename:
                line=batch+','+run+','+lineage+','+lineage+'.bam,'+rest_of_line
            else:
                line=batch+','+run+','+lineage+"_"+uid+','+lineage+"_"+uid+'.bam,'+rest_of_line
        print(line)
