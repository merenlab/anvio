#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

import sys

import anvio
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.mcgclassifier import MetagenomeCentricGeneClassifier


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ShaiberAlon']
__description__ = "A program to classify genes according to coverage across multiple metagenomes"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()

    mcg = MetagenomeCentricGeneClassifier(args)

    if args.output_file_prefix:
        mcg.write_output_to_files = True

    args.init_gene_coverages = False if args.get_samples_stats_only else True
    args.init_split_coverage_values_per_nt = True
    args.populate_nt_level_coverage = True
    summary = summarizer.ProfileSummarizer(args)
    summary.init()
    if not args.init_gene_coverages:
        summary.init_split_coverage_values_per_nt_dict()

    bin_names_in_collection = summary.bin_ids
    if args.bin_ids_file:
        filesnpaths.is_file_exists(args.bin_ids_file)
        bin_names_of_interest = [line.strip() for line in open(args.bin_ids_file).readlines()]

        missing_bins = [b for b in bin_names_of_interest if b not in bin_names_in_collection]
        if len(missing_bins):
            raise ConfigError("Some bin names you declared do not appear to be in the collection %s. "
                               "These are the bins that are missing: %s, these are the bins that are "
                               "actually in your collection: %s" % (args.collection_name, missing_bins, bin_names_in_collection))
    elif args.bin_id:
        if args.bin_id not in bin_names_in_collection:
            raise ConfigError("The bin you declared, %s, does not appear to be in the collection %s." \
                              % (args.bin_id, args.collection_name))
        bin_names_of_interest = [args.bin_id]
    else:
        bin_names_of_interest = bin_names_in_collection

    for bin_name in bin_names_of_interest:
        # This is where things are happening
        _bin = summarizer.Bin(summary, bin_name)
        mcg.init(gene_level_coverage_stats_dict = _bin.gene_level_coverage_stats_dict, \
                    split_coverage_values_per_nt_dict = _bin.split_coverage_values_per_nt_dict, \
                    additional_description = bin_name)
        if not args.get_samples_stats_only:
            mcg.classify()
        else:
            mcg.init_samples_coverage_stats_dict()
            mcg.save_samples_information()


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)
    group0 = parser.add_argument_group('ESSENTIAL INPUTS', "You must supply a merged profile db (along with a matching contigs db)")
    group0.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    group0.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    group1 = parser.add_argument_group('ESSENTIAL OUTPUTS', 'The outputs of the algorithm are: an anvio additional layers format file with the classification \
                                                            information for genes. An anvio samples information file with detectino information per sample. \
                                                            In addition, when a profile database is given then a gene-coverages, and gene-detection tables \
                                                            would also be saved. All files are created with the prefix that is provided by the user.')
    group1.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    group1.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB = parser.add_argument_group('ADDITIONAL STUFF',"Parameters to provide pre-existing additional layers, samples-information files, so that the \
                                                            outputs would be added to these files")
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id', {'required': False}))
    groupB.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file', {'required': False}))
    groupB.add_argument('--exclude-samples', metavar='FILE', default=None, help='List of samples to exclude for the analysis.')
    groupB.add_argument('--include-samples', metavar='FILE', default=None, help='List of samples to include for the analysis.')
    groupB.add_argument(*anvio.A('gen-figures'), **anvio.K('gen-figures'))
    groupB.add_argument(*anvio.A('get-samples-stats-only'), **anvio.K('get-samples-stats-only'))
    groupB.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))
    groupC = parser.add_argument_group("PARAMETERS","Parameters to determine cut-offs for the gene-classifier")
    groupC.add_argument('--alpha', '--genome-detection-uncertainty', metavar='NUM', type=float, default=0.25,
                        help='Determines the range of sample detection values that are considered negative, ambiguous or positive. \
                        Min of 0 and smaller than 0.5, default of 0.25. For exmaple for the default samples with detection \
                        below  0.5-0.25 = 0.25 will be considered negative (i.e. donot contain the genome), samples with detection \
                        between 0.25 and 0.75 would be ambiguous (and hence would not be used for the classification), and samples with \
                        detection above 0.75 would be considered positive (i.e. contain the genome).')
    groupC.add_argument(*anvio.A('outliers-threshold'), **anvio.K('outliers-threshold'))
    groupC.add_argument(*anvio.A('zeros-are-outliers'), **anvio.K('zeros-are-outliers'))
    return parser.get_args(parser)


if __name__ == '__main__':
    main()
