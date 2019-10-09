#!/usr/bin/env python
import os.path
import sys
import tempfile
import argparse
import oncogemini.version

def examples(parser, args):

    print( "[stats] - report basic statistics about your variants:")
    print( "   oncogemini stats --tstv my.db")
    print( "   oncogemini stats --tstv-coding my.db")
    print( "   oncogemini stats --sfs my.db")
    print( "   oncogemini stats --snp-counts my.db")
    print("")

    print( "[query] - explore the database with ad hoc queries:")
    print( "   oncogemini query -q \"select * from variants where is_lof = 1 and aaf <= 0.01\" my.db")
    print( "   oncogemini query -q \"select chrom, pos, gt_bases.NA12878 from variants\" my.db")
    print( "   oncogemini query -q \"select chrom, pos, in_omim, clin_sigs from variants\" my.db")
    print("")

    print( "[dump] - convenient \"data dumps\":")
    print( "   oncogemini dump --variants my.db")
    print( "   oncogemini dump --genotypes my.db")
    print( "   oncogemini dump --samples my.db")
    print("")

    print( "[region] - access variants in specific genomic regions:")
    print( "   oncogemini region --reg chr1:100-200 my.db")
    print( "   oncogemini region --gene TP53 my.db")
    print("")

    print( "[tools] - there are also many specific tools available")
    print( "   1. Find truncal variants.")
    print( "     oncogemini truncal my.db")
    print("")

    exit()

def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='oncogemini', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed oncogemini version",
                        action="version",
                        version="%(prog)s " + str(oncogemini.version.__version__))
    parser.add_argument('--annotation-dir', dest='annotation_dir',
                             help='Path to the annotation database.\n'
                                'This argument is optional and if given will take precedence over the default location stored in the oncogemini config file.')
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    #########################################
    # $ gemini examples
    #########################################
    parser_examples = subparsers.add_parser('examples',
                                            help='show usage examples')
    parser_examples.set_defaults(func=examples)

    #########################################
    # $ gemini amend
    #########################################
    parser_amend = subparsers.add_parser('amend',
                                         help="Amend an already loaded OncoGEMINI database.")
    parser_amend.add_argument('db',
                              metavar='db',
                              help='The name of the database to be amended.')
    parser_amend.add_argument('--sample',
                              metavar='sample',
                              default=None,
                              help='New sample information file to load')
    parser_amend.add_argument('--clear',
                              default=False,
                              action="store_true",
                              help='Set all values in this column to NULL before loading.')

    def amend_fn(parser, args):
        from oncogemini import gemini_amend
        gemini_amend.amend(parser, args)
    parser_amend.set_defaults(func=amend_fn)

    #########################################
    # $ gemini query
    #########################################
    parser_query = subparsers.add_parser('query',
            help='issue ad hoc SQL queries to the DB')
    parser_query.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_query.add_argument('-q',
            dest='query',
            metavar='QUERY_STR',
            help='The query to be issued to the database')
    parser_query.add_argument('--gt-filter',
            dest='gt_filter',
            metavar='STRING',
            help='Restrictions to apply to genotype values')
    parser_query.add_argument('--show-samples',
                              dest='show_variant_samples',
                              action='store_true',
                              default=False,
                              help=('Add a column of all sample names with a variant to each '
                                    'variant.'))
    parser_query.add_argument('--show-families',
                              dest='show_families',
                              action='store_true',
                              default=False,
                              help=('Add a column listing all of the families '
                                    'with a variant to each variant.'))
    parser_query.add_argument('--family-wise',
                              dest='family_wise',
                              default=False,
                              action='store_true',
                              help=('Perform the sample-filter on a family-wise '
                                    'basis.'))
    parser_query.add_argument('--min-kindreds',
                              dest='min_kindreds',
                              default=1,
                              type=int,
                              help=('Minimum number of families for a variant passing '
                                    'a family-wise filter to be in.'))
    parser_query.add_argument('--sample-delim',
            dest='sample_delim',
            metavar='STRING',
            help='The delimiter to be used with the --show-samples option.',
            default=',')

    parser_query.add_argument('--header',
            dest='use_header',
            action='store_true',
            help='Add a header of column names to the output.',
            default=False)
    parser_query.add_argument('--sample-filter',
                              dest='sample_filter',
                              help='SQL filter to use to filter the sample table',
                              default=None)
    parser_query.add_argument('--in',
                              dest='in_subject',
                              nargs='*',
                              help=('A variant must be in either all, none or any '
                                    'samples passing the --sample-query filter.'),
                              choices=['all', 'none', 'any', 'only', 'not'],
                              default=['any'])
    parser_query.add_argument('--format',
                              dest='format',
                              default='default',
                              help='Format of output (JSON, TPED or default)')
    parser_query.add_argument('--region',
                              dest='region',
                              default=None,
                              help=('Restrict query to this region, '
                                    'e.g. chr1:10-20.'))
    parser_query.add_argument('--carrier-summary-by-phenotype',
                              dest='carrier_summary',
                              default=None,
                              help=('Output columns of counts of carriers and '
                                    'non-carriers stratified by the given '
                                    'sample phenotype column'))

    def query_fn(parser, args):
        from oncogemini import gemini_query
        gemini_query.query(parser, args)

    parser_query.set_defaults(func=query_fn)

    #########################################
    # $ gemini dump
    #########################################
    parser_dump = subparsers.add_parser('dump',
            help='shortcuts for extracting data from the DB')
    parser_dump.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_dump.add_argument('--variants',
            dest='variants',
            action='store_true',
            help='Report all rows/columns from the variants table.',
            default=False)
    parser_dump.add_argument('--genotypes',
            dest='genotypes',
            action='store_true',
            help='Report all rows/columns from the variants table \nwith one line per sample/genotype.',
            default=False)
    parser_dump.add_argument('--samples',
            dest='samples',
            action='store_true',
            help='Report all rows/columns from the samples table.',
            default=False)
    parser_dump.add_argument('--header',
            dest='use_header',
            action='store_true',
            help='Add a header of column names to the output.',
            default=False)
    parser_dump.add_argument('--sep',
            dest='separator',
            metavar='STRING',
            help='Output column separator',
            default="\t")
    parser_dump.add_argument('--tfam',
                             dest='tfam',
                             action='store_true',
                             default=False,
                             help='Output sample information to TFAM format.')
    def dump_fn(parser, args):
        from oncogemini import gemini_dump
        gemini_dump.dump(parser, args)
    parser_dump.set_defaults(func=dump_fn)

    #########################################
    # $ gemini region
    #########################################
    parser_region = subparsers.add_parser('region',
            help='extract variants from specific genomic loci')
    parser_region.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_region.add_argument('--reg',
            dest='region',
            metavar='STRING',
            help='Specify a chromosomal region chr:start-end')
    parser_region.add_argument('--gene',
            dest='gene',
            metavar='STRING',
            help='Specify a gene of interest')
    parser_region.add_argument('--header',
            dest='use_header',
            action='store_true',
            help='Add a header of column names to the output.',
            default=False)
    parser_region.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_region.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_region.add_argument('--show-samples',
                               dest='show_variant_samples',
                               action='store_true',
                               default=False,
                                help=('Add a column of all sample names with a variant to each '
                                      'variant.'))
    parser_region.add_argument('--format',
                              dest='format',
                              default='default',
                              help='Format of output (JSON, TPED or default)')
    def region_fn(parser, args):
        from oncogemini import gemini_region
        gemini_region.region(parser, args)
    parser_region.set_defaults(func=region_fn)

    #########################################
    # $ gemini stats
    #########################################
    parser_stats = subparsers.add_parser('stats',
            help='compute useful variant stastics')
    parser_stats.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_stats.add_argument('--tstv',
            dest='tstv',
            action='store_true',
            help='Report the overall ts/tv ratio.',
            default=False)
    parser_stats.add_argument('--tstv-coding',
            dest='tstv_coding',
            action='store_true',
            help='Report the ts/tv ratio in coding regions.',
            default=False)
    parser_stats.add_argument('--tstv-noncoding',
            dest='tstv_noncoding',
            action='store_true',
            help='Report the ts/tv ratio in non-coding regions.',
            default=False)
    parser_stats.add_argument('--snp-counts',
            dest='snp_counts',
            action='store_true',
            help='Report the count of each type of SNP (A->G, G->T, etc.).',
            default=False)
    parser_stats.add_argument('--sfs',
            dest='sfs',
            action='store_true',
            help='Report the site frequency spectrum of the variants.',
            default=False)
    parser_stats.add_argument('--mds',
            dest='mds',
            action='store_true',
            help='Report the pairwise genetic distance between the samples.',
            default=False)
    parser_stats.add_argument('--vars-by-sample',
            dest='variants_by_sample',
            action='store_true',
            help='Report the number of variants observed in each sample.',
            default=False)
    parser_stats.add_argument('--gts-by-sample',
            dest='genotypes_by_sample',
            action='store_true',
            help='Report the count of each genotype class obs. per sample.',
            default=False)
    parser_stats.add_argument('--summarize',
            dest='query',
            metavar='QUERY_STR',
            default=None,
            help='The query to be issued to the database to summarize')
    parser_stats.add_argument('--gt-filter',
            dest='gt_filter',
            metavar='STRING',
            help='Restrictions to apply to genotype values')
    def stats_fn(parser, args):
        from oncogemini import gemini_stats
        gemini_stats.stats(parser, args)
    parser_stats.set_defaults(func=stats_fn)

    #########################################
    # gemini annotate
    #########################################
    parser_get = subparsers.add_parser('annotate',
            help='Add new columns for custom annotations')
    parser_get.add_argument('db',
            metavar='db',
            help='The name of the database to be updated.')
    parser_get.add_argument('-f',
            dest='anno_file',
            required=True,
            help='The TABIX\'ed BED file containing the annotations')
    parser_get.add_argument('-c',
            dest='col_names',
            help='The name(s) of the BED column(s) to be added to the variant table.'
            'If the input file is a VCF, then this is the name of the info field to pull.')
    parser_get.add_argument('-a',
            dest='anno_type',
            help='How should the annotation file be used? (def. extract)',
            default="extract",
            choices=['boolean', 'count', 'extract'])
    parser_get.add_argument('-e',
            dest='col_extracts',
            help='Column(s) to extract information from for list annotations.'
            'If the input is VCF, then this defaults to the fields specified in `-c`.')
    parser_get.add_argument('-t',
            dest='col_types',
            help='What data type(s) should be used to represent the new values '
                 'in the database? '
                 'Any of {integer, float, text}')
    parser_get.add_argument('-o',
            dest='col_operations',
            help='Operation(s) to apply to the extract column values '
                  'in the event that a variant overlaps multiple annotations '
                  'in your annotation file (-f).'
                  'Any of {sum, mean, median, min, max, mode, list, uniq_list, first, last}')
    parser_get.add_argument('--region-only',
            dest='region_only',
            action='store_true',
            default=False,
            help='If set, only region coordinates will be considered when annotating variants.'
                 'The default is to annotate using region coordinates as well as REF and ALT variant values'
                 'This option is only valid if annotation is a VCF file')
    def annotate_fn(parser, args):
        from oncogemini import gemini_annotate
        gemini_annotate.annotate(parser, args)
    parser_get.set_defaults(func=annotate_fn)

    #########################################
    # gemini windower
    #########################################
    parser_get = subparsers.add_parser('windower',
            help='Compute statistics across genome \"windows\"')
    parser_get.add_argument('db',
            metavar='db',
            help='The name of the database to be updated.')
    parser_get.add_argument('-w',
            dest='window_size',
            type=int,
            default=1000000,
            help='The name of the column to be added to the variant table.')
    parser_get.add_argument('-s',
            dest='step_size',
            type=int,
            default=0,
            help="The step size for the windows in bp.\n")
    parser_get.add_argument('-t',
            dest='analysis_type',
            help='The type of windowed analysis requested.',
            choices=['nucl_div', 'hwe'],
            default='hwe')
    parser_get.add_argument('-o',
            dest='op_type',
            help='The operation that should be applied to the -t values.',
            choices=['mean', 'median', 'min', 'max', 'collapse'],
            default='mean')
    def windower_fn(parser, args):
        from oncogemini import gemini_windower
        gemini_windower.windower(parser, args)
    parser_get.set_defaults(func=windower_fn)

    #########################################
    # gemini db_info
    #########################################
    parser_get = subparsers.add_parser('db_info',
            help='Get the names and types of cols. database tables')
    parser_get.add_argument('db',
            metavar='db',
            help='The name of the database to be updated.')
    def dbinfo_fn(parser, args):
        from oncogemini import gemini_dbinfo
        gemini_dbinfo.db_info(parser, args)
    parser_get.set_defaults(func=dbinfo_fn)

    #########################################
    # $ gemini set_somatic
    #########################################
    parser_set_somatic = subparsers.add_parser("set_somatic",
                          help="Tag somatic mutations (is_somatic) by comparint tumor/normal pairs.")
    parser_set_somatic.add_argument('db', metavar='db',
            help='The name of the database to be updated.')

    parser_set_somatic.add_argument('--minDP',
            dest='minDP',
            type=int,
            default=0,
            help='Minimum depth required in all samples (default is 0)')

    parser_set_somatic.add_argument('--minGQ',
            dest='minGQ',
            type=float,
            default=0,
            help='Minimum genotype quality required in all samples (default is 0)')

#    parser_set_somatic.add_argument('--min-somatic-score',
#            dest='min_somatic_score',
#            type=float,
#            default=0,
#            help='The min somatic score (SSC) (def: %(default)s).')

    parser_set_somatic.add_argument('--normAF',
            dest='normAF',
            type=float,
            default=0,
            help='The max freq. of the alt. allele in the normal sample (def: %(default)s).')

    parser_set_somatic.add_argument('--normCount',
            dest='normCount',
            type=int,
            default=0,
            help='The max count. of the alt. allele in the normal sample (def: %(default)s).')

    parser_set_somatic.add_argument('--normDP',
            dest='normDP',
            type=int,
            default=0,
            help='The minimum depth allowed in the normal sample to believe somatic (def: %(default)s).')

    parser_set_somatic.add_argument('--tumAF',
            dest='tumAF',
            type=float,
            default=0,
            help='The min freq. of the alt. allele in the tumor sample (def: %(default)s).')

    parser_set_somatic.add_argument('--tumCount',
            dest='tumCount',
            type=int,
            default=0,
            help='The min count. of the alt. allele in the tumor sample (def: %(default)s).')

    parser_set_somatic.add_argument('--tumDP',
            dest='tumDP',
            type=int,
            default=0,
            help='The minimum depth allowed in the tumor sample to believe somatic (def: %(default)s).')

#    parser_set_somatic.add_argument('--chrom',
#            dest='chrom',
#            metavar='STRING',
#            help='A specific chromosome on which to tag somatic mutations. (def: %(default)s).',
#            default=None,
#            )

    parser_set_somatic.add_argument('--dry-run',
            dest='dry_run',
            action='store_true',
            help='Don\'t set the is_somatic flag, just report what _would_ be set. For testing parameters.',
            default=False)

    def set_somatic_fn(parser, args):
        from oncogemini import gemini_set_somatic
        gemini_set_somatic.set_somatic(parser, args)
    parser_set_somatic.set_defaults(func=set_somatic_fn)

    #########################################
    # $ gemini update
    #########################################
    parser_update = subparsers.add_parser("update", help="Update oncogemini software and data files.")
    parser_update.add_argument("--devel", help="Get the latest development version instead of the release",
                               action="store_true", default=False)
    parser_update.add_argument("--dataonly", help="Only update data, not the underlying libraries.",
                               action="store_true", default=False)
    parser_update.add_argument("--nodata", help="Do not install data dependencies",
                               dest="install_data", action="store_false", default=True)
    parser_update.add_argument("--extra", help="Add additional non-standard genome annotations to include",
                               action="append", default=[], choices=["gerp_bp","cadd_score"])
    parser_update.add_argument("--tooldir", help="Directory for third party tools (ie /usr/local) update")
    def update_fn(parser, args):
        from oncogemini import gemini_update
        gemini_update.release(parser, args)
    parser_update.set_defaults(func=update_fn)


    #########################################
    # $ gemini roh
    #########################################
    parser_hom_run = subparsers.add_parser('roh',
            help='Identify runs of homozygosity')
    parser_hom_run.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_hom_run.add_argument('--min-snps',
            dest='min_snps',
            metavar="INTEGER",
            type=int,
            default=25,
            help='Minimum number of homozygous snps expected in a run (def. 25)')
    parser_hom_run.add_argument('--min-total-depth',
            dest='min_total_depth',
            metavar="INTEGER",
            type=int,
            default=20,
            help="""The minimum overall sequencing depth required"""
                 """for a SNP to be considered (def = 20).""")
    parser_hom_run.add_argument('--min-gt-depth',
            dest='min_genotype_depth',
            metavar="INTEGER",
            type=int,
            default=0,
            help="""The minimum required sequencing depth underlying a given sample's genotype"""
                 """for a SNP to be considered (def = 0).""")
    parser_hom_run.add_argument('--min-size',
            metavar="INTEGER",
            dest='min_size',
            type=int,
            default=100000,
            help='Minimum run size in base pairs (def. 100000)')
    parser_hom_run.add_argument('--max-hets',
            metavar="INTEGER",
            dest='max_hets',
            type=int,
            default=1,
            help='Maximum number of allowed hets in the run (def. 1)')
    parser_hom_run.add_argument('--max-unknowns',
            metavar="INTEGER",
            type=int,
            dest='max_unknowns',
            default=3,
            help='Maximum number of allowed unknowns in the run (def. 3)')
    parser_hom_run.add_argument('-s',
            dest='samples',
            default=None,
            help='Comma separated list of samples to screen for ROHs. e.g S120,S450')
    def homozygosity_runs_fn(parser, args):
        from gemini.tool_homozygosity_runs import run
        run(parser, args)
    parser_hom_run.set_defaults(func=homozygosity_runs_fn)


    #########################################
    # $ gemini fusions
    #########################################
    parser_fusions = subparsers.add_parser('fusions',
                                         help="Identify somatic fusion genes from a OncoGEMINI database.")
    parser_fusions.add_argument('db',
                              metavar='db',
                              help='The name of the database to be queried.')
    parser_fusions.add_argument('--in_cosmic_census',
                                action='store_true',
                                help='One or both genes in fusion is in COSMIC cancer census')
    parser_fusions.add_argument('--min_qual',
                                dest='min_qual',
                                metavar='FLOAT',
                                type=float,
                                default=None,
                                help='The min variant quality (VCF QUAL) (def: %(default)s).')
    parser_fusions.add_argument('--evidence_type',
                                metavar='STR',
                                dest='evidence_type',
                                type=str,
                                default=None,
                                help='The supporting evidence types for the variant ("PE", "SR", or "PE,SR").')

    def fusions_fn(parser, args):
        from oncogemini.tool_fusions import run
        run(parser, args)
    parser_fusions.set_defaults(func=fusions_fn)

    #########################################
    # $ gemini truncal
    #########################################
    parser_truncal = subparsers.add_parser('truncal',
            help='Report truncal mutations')
    parser_truncal.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_truncal.add_argument('--minDP',
            dest='minDP',
            metavar='INTEGER',
            help='Minimum depth required in all samples default is 0)')
    parser_truncal.add_argument('--minGQ',
            dest='minGQ',
            metavar='INTEGER',
            help='Minimum genotype quality required in all samples (default is 0)')
    parser_truncal.add_argument('--maxNorm',
            dest='maxNorm',
            metavar='FLOAT',
            help='Optional: specify a maximum normal sample AF to allow (default is 0)')
    parser_truncal.add_argument('--patient',
            dest='patient',
            metavar='STRING',
            help='Specify a patient to filter (should correspond to a patient_id in ped file)')
    parser_truncal.add_argument('--samples',
            dest='samples',
            metavar='STRING',
            help='Optional: rather than including all samples, a string of comma-separated specified samples to use (default is "All")')
    parser_truncal.add_argument('--increase',
            dest='increase',
            metavar='FLOAT',
            help='Optional: add amount to increase truncal AF filter between normal and tumor samples (default is 0)')
    parser_truncal.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned (default is "*", which returns every column)')
    parser_truncal.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_truncal.add_argument('--purity',
            action="store_true",
            help='Using purity estimates in ped file, make corrections to AF to be used')
    parser_truncal.add_argument('--somatic_only',
            action="store_true",                       
            help='Only include variants that have been marked as somatic via the set_somatic command')
    parser_truncal.add_argument('--cancers',
            dest='cancers',
            metavar='STRING',
            help='Restrict results to variants/genes associated with specific cancer types by entering a comma-separated string of cancer type abbreviations (see documents for abbreviations) REQUIRES that db include civic_gene_abbrevations and/or cgi_gene_abbreviations')

    def truncal_fn(parser, args):
        from oncogemini.gemini_truncal import truncal
        truncal(parser, args)
    parser_truncal.set_defaults(func=truncal_fn)

    #########################################
    # $ gemini loh
    #########################################
    parser_loh = subparsers.add_parser('loh',
            help='Report loss of heterozygosity (LOH) mutations between normal tissue and all tumor samples')
    parser_loh.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_loh.add_argument('--minDP',
            dest='minDP',
            metavar='INTEGER',
            help='Minimum depth required in all samples default is 0)')
    parser_loh.add_argument('--minGQ',
            dest='minGQ',
            metavar='INTEGER',
            help='Minimum genotype quality required in all samples (default is 0)')
    parser_loh.add_argument('--maxNorm',
            dest='maxNorm',
            metavar='FLOAT',
            help='Specify a maximum normal sample AF to allow (default is 0.7)')
    parser_loh.add_argument('--minNorm',
            dest='minNorm',
            metavar='FLOAT',
            help='Specify a minimum normal sample AF to allow (default is 0.3)')
    parser_loh.add_argument('--minTumor',
            dest='minTumor',
            metavar='FLOAT',
            help='Specify a minimum AF for tumor samples to require (default is 0.8)')
    parser_loh.add_argument('--patient',
            dest='patient',
            metavar='STRING',
            help='Specify a patient to filter (should correspond to a patient_id in ped file)')
    parser_loh.add_argument('--samples',
            dest='samples',
            metavar='STRING',
            help='Rather than including all samples, enter a string of comma-separated specified samples to use (default is "All")')
    parser_loh.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A comma-separated list of columns that you would like returned (default is "*", which returns every column)')
    parser_loh.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_loh.add_argument('--purity',
            action="store_true",
            help='Using purity estimates in cancer manidest, make corrections to AF to be used')
    parser_loh.add_argument('--specific',
            dest='specific',
            metavar='STRING',
            help='Search for LOH variants in a single sample compared to the sample(s) that precede it (must specify single sample included among --samples, also --minNorm, --maxNorm will now apply to the preceding sample)')
    parser_loh.add_argument('--cancers',
            dest='cancers',
            metavar='STRING',
            help='Restrict results to variants/genes associated with specific cancer types by entering a comma-separated string of cancer type abbreviations (see documents for abbreviations) REQUIRES that db include civic_gene_abbrevations and/or cgi_gene_abbreviations')

    def loh_fn(parser, args):
        from oncogemini.gemini_loh import loh
        loh(parser, args)
    parser_loh.set_defaults(func=loh_fn)

    #########################################
    # $ gemini bottleneck
    #########################################
    parser_bottleneck = subparsers.add_parser('bottleneck',
            help='Report bottleneck mutations')
    parser_bottleneck.add_argument('db',
            metavar='db',
            help='The name of the database to be queried')
    parser_bottleneck.add_argument('--minDP',
            dest='minDP',
            metavar='INTEGER',                       
            help='Minimum depth required in all samples default is 0)')                       
    parser_bottleneck.add_argument('--minGQ',
            dest='minGQ',
            metavar='INTEGER',
            help='Minimum genotype quality required in all samples (default is 0)')
    parser_bottleneck.add_argument('--maxNorm',
            dest='maxNorm',
            metavar='FLOAT',
            help='Specify a maximum normal sample AF to allow (default is 0)')
    parser_bottleneck.add_argument('--minSlope',
            dest='minSlope',
            metavar='FLOAT',
            help='Minimum slope required for the AFs across samples (default is 0.05)')
    parser_bottleneck.add_argument('--minR',
            dest='minR',
            metavar='FLOAT',
            help='Minimum r correlation coefficient required for AFs (default is 0.5)')
    parser_bottleneck.add_argument('--samples',
            dest='samples',
            metavar='STRING',
            help='Rather than including all samples, a string of comma-separated specified samples to use (default is "All")')
    parser_bottleneck.add_argument('--minEnd',
            dest='minEnd',
            metavar='FLOAT',
            help='Minimum AF required of the sample representing the final timepoint (default is 0)')                        
    parser_bottleneck.add_argument('--endDiff',
            dest='endDiff',
            metavar='FLOAT',
            help='Minimum required AF difference between the samples representing the first and final timepoints (default is 0)')
    parser_bottleneck.add_argument('--patient',
            dest='patient',
            metavar='STRING',
            help='Specify a patient to filter (should correspond to a patient_id in ped file)')
    parser_bottleneck.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned (default is "*", which returns every column)')
    parser_bottleneck.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_bottleneck.add_argument('--purity',
            action="store_true",                       
            help='Using purity estimates in ped file, make corrections to AF to be used')
    parser_bottleneck.add_argument('--somatic_only',
            action="store_true",                       
            help='Only include variants that have been marked as somatic via the set_somatic command')
    parser_bottleneck.add_argument('--cancers',
            dest='cancers',
            metavar='STRING',
            help='Restrict results to variants/genes associated with specific cancer types by entering a comma-separated string of cancer type abbreviations (see documents for abbreviations) REQUIRES that db include civic_gene_abbrevations and/or cgi_gene_abbreviations')

    def bottleneck_fn(parser, args):
        from oncogemini.gemini_bottleneck import bottleneck
        bottleneck(parser, args)
    parser_bottleneck.set_defaults(func=bottleneck_fn)

    #########################################
    # $ gemini unique
    #########################################
    parser_unique = subparsers.add_parser('unique',
            help='Report unique mutations for specified sample(s)')
    parser_unique.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_unique.add_argument('--minDP',
            dest='minDP',
            metavar='INTEGER',
            help='Minimum depth required in all samples default is 0)')
    parser_unique.add_argument('--minGQ',
            dest='minGQ',
            metavar='INTEGER',
            help='Minimum genotype quality required in all samples (default is 0)')
    parser_unique.add_argument('--specific',
            dest='specific',
            metavar='STRING',
            help='Identify unique variants that exist only in samples from this comma-separated list')
    parser_unique.add_argument('--maxOthers',
            dest='maxOthers',
            metavar='FLOAT',
            help='Specify a maximum sample AF to allow in other samples (default is 0)')
    parser_unique.add_argument('--patient',
            dest='patient',
            metavar='STRING',
            help='Specify a patient to filter (should correspond to a patient_id in ped file)')
    parser_unique.add_argument('--samples',
            dest='samples',
            metavar='STRING',
            help='Rather than including all samples in filters, a string of comma-separated specified samples to use (default is "All")')
    parser_unique.add_argument('--increase',
            dest='increase',
            metavar='FLOAT',
            help='Add amount to increase AF filter between unique and other samples (default is 0)')
    parser_unique.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned (default is "*", which returns every column)')
    parser_unique.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_unique.add_argument('--purity',
            action="store_true",
            help='Using purity estimates in ped file, make corrections to AF to be used')
    parser_unique.add_argument('--somatic_only',
            action="store_true",                       
            help='Only include variants that have been marked as somatic via the set_somatic command')
    parser_unique.add_argument('--cancers',
            dest='cancers',
            metavar='STRING',
            help='Restrict results to variants/genes associated with specific cancer types by entering a comma-separated string of cancer type abbreviations (see documents for abbreviations) RE\
QUIRES that db include civic_gene_abbrevations and/or cgi_gene_abbreviations')

    def unique_fn(parser, args):
        from oncogemini.gemini_unique import unique
        unique(parser, args)
    parser_unique.set_defaults(func=unique_fn)


    #######################################################
    # parse the args and call the selected function
    #######################################################
    import operator
    subparsers._choices_actions.sort(key=operator.attrgetter('dest'))
    for k in sorted(subparsers.choices):
        subparsers.choices[k] = subparsers.choices.pop(k)

    args = parser.parse_args()

    # make sure database is found if provided
    try:
        args.func(parser, args)
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

def xor(arg1, arg2):
    return bool(arg1) ^ bool(arg2)


if __name__ == "__main__":
    main()
