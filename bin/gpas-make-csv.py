import glob, uuid, os, argparse, random, pkg_resources

from datetime import date, timedelta

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--organisation",default="University of Oxford",help="the name of the organisation (the user must belong to it otherwise validation will fail)")
    parser.add_argument("--country",default="United Kingdom",help="the name of the country where the samples were collected")
    parser.add_argument("--tech",default='illumina',help="whether to generate illumina (paired) or nanopore (unpaired) reads")
    parser.add_argument("--file_type",default='fastq',help="whether to look for FASTQ or BAM files")
    parser.add_argument("--tag_file",default=pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/tags.txt'),help="a plaintext file with one row per tag")
    parser.add_argument("--uuid_length",default='long',help="whether to use a long or short UUID4")
    parser.add_argument("--number_of_tags",default=2,type=int,help="how many tags to give each sample. Can be zero, or up to the number of rows in <tag_file>. Default is 2 so as to test the delimiter")
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
            header='name,fastq1,fastq2,organisation,tags,specimenOrganism,host,collectionDate,country,submissionTitle,submissionDescription,instrument_platform,instrument_model,flowcell'
            file_list = glob.glob('*_1.fastq.gz')
            file_extensions = ['_1.fastq.gz','_2.fastq.gz']
        elif options.tech == 'nanopore':
            header='name,fastq,organisation,tags,specimenOrganism,host,collectionDate,country,submissionTitle,submissionDescription,instrument_platform,instrument_model,flowcell'
            file_list = glob.glob('*.fastq.gz')
            file_extensions = ['.fastq.gz']
            file_extensions = ['_1.fastq.gz','_2.fastq.gz']
    elif options.file_type == 'bam':
        header='name,bam,organisation,tags,specimenOrganism,host,collectionDate,country,submissionTitle,submissionDescription,instrument_platform,instrument_model,flowcell'
        file_list = glob.glob('*.bam')
        file_extensions = ['.bam']

    print(header)

    def build_rest_of_line():
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

    for i in file_list:

        if options.file_type == 'fastq':
            if options.tech=='illumina':
                filename=i.split('_1.fastq.gz')[0]
            elif options.tech=='nanopore':
                filename=i.split('.fastq.gz')[0]
        else:
            filename=i.split('.bam')[0]

        if '_' in filename:
            lineage=i.split('_')[0]
        else:
            lineage=filename

        if options.uuid_length == 'long':
            uid=str(uuid.uuid4())
        else:
            uid=str(uuid.uuid4())[-4:]

        rest_of_line=build_rest_of_line()

        for file_extension in file_extensions:
            os.rename(filename+file_extension,lineage+"_"+uid+file_extension)

        if options.file_type == 'fastq':
            if options.tech=='illumina':
                line=lineage+"_"+uid+','+lineage+"_"+uid+'_1.fastq.gz,'+lineage+"_"+uid+'_2.fastq.gz,'+rest_of_line
            elif options.tech=='nanopore':
                line=lineage+"_"+uid+','+lineage+"_"+uid+'.fastq.gz,'+rest_of_line
        else:
            line=lineage+"_"+uid+','+lineage+"_"+uid+'.bam,'+rest_of_line

        print(line)
