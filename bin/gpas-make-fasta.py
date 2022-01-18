

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--output",required=False,help="the stem of the output file")
    parser.add_argument("--variant_name",default='Reference',help="the name of the variant, default is Reference")
    parser.add_argument("--reference",required=False,default=pkg_resources.resource_filename("gpas_covid_synthetic_reads", 'data/MN908947.3.gbk'),help="the GenBank file of the covid reference (if not specified, the MN908947.3.gbk reference will be used)")
    parser.add_argument("--write_fasta", action="store_true", help="whether to write out the FASTA file for the variant")
    options = parser.parse_args()
