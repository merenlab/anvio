#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.parsers as parsers
import anvio.terminal as terminal

from anvio.constants import levels_of_taxonomy
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import TableForLayerAdditionalData


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['single-profile-db', 'layer-taxonomy-txt',]
__provides__ = ['layer-taxonomy',]
__description__ = ("Import layers-level taxonomy into an anvi'o additional layer data table in an "
                   "anvi'o single-profile database")


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
    run = terminal.Run()
    progress = terminal.Progress()

    if not args.input_files:
        raise ConfigError("You need to use '--input-files' parameter to list file(s) that is/are required the "
                           "parser you chose. Please see the documentation for details.")

    if utils.is_profile_db_merged(args.profile_db) and not utils.is_blank_profile(args.profile_db):
        raise ConfigError("This will not work with merged profile databases :(")

    if args.min_abundance < 0 or args.min_abundance > 5:
        raise ConfigError("The minimum abundance value for taxon names to be included should be equal or "
                          "greater than 0%, and lower than 5% :/")

    parser_obj = parsers.get_parser_obj('taxonomy_layers', args.parser)

    parser = parser_obj(args.input_files, t.splits_taxonomy_table_structure, run=run, progress=progress)
    taxonomy_dict = parser.process()

    if not len(taxonomy_dict):
        raise ConfigError("Your parser (%s) returned an empty result for your input file. Something must have gone wrong." % (args.parser))

    # build input for layer additional data tables
    sample_name = dbops.ProfileDatabase(args.profile_db).meta['sample_id']
    for level in levels_of_taxonomy:
        taxonomy_dict_level = {}

        total_count = sum(taxonomy_dict[level].values())
        taxa_excluded = set([])
        taxa_excluded_count = 0

        for taxon in taxonomy_dict[level]:
            if taxonomy_dict[level][taxon] * 100 / total_count < args.min_abundance:
                taxa_excluded.add(taxon)
                taxa_excluded_count += taxonomy_dict[level][taxon]
            else:
                taxonomy_dict_level['%s!%s' % (level, taxon)] = taxonomy_dict[level][taxon]

        data_keys_list = list(taxonomy_dict_level.keys())
        data_dict = {sample_name: taxonomy_dict_level}

        args.target_data_group = level
        T = TableForLayerAdditionalData(args, r=terminal.Run(verbose=False))

        # first remove previous entries for `data_group`
        data_keys_in_db, data_dict_in_db = T.get()
        if len(data_keys_in_db):
            T.remove(data_keys_list=data_keys_in_db)

        # then add the incoming data
        T.add(data_dict, data_keys_list)

        run.info_single("%d entries added into layers additional data as data group '%s'" % (len(data_keys_list), level))

        if taxa_excluded_count:
            run.warning("A total of %d taxon names, which made up %.2f%% of the input data altogether were excluded "
                        "from this import due to the min abundance criterion of %.1f%% for inclusion." % \
                                (len(taxa_excluded), taxa_excluded_count * 100 / total_count, args.min_abundance))


def get_args():
    parser_names = parsers.get_parser_names('taxonomy_layers')

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))

    parser.add_argument('--parser', default = None,
                        help="Parser to make sense of the input files. There are %d parsers readily available: %s." % (len(parser_names), parser_names))
    parser.add_argument('-i', '--input-files', metavar = 'FILE(S)', nargs='+', default=None, required=True,
                        help='Input file(s) for selected parser. Each parser (except "blank") requires input files to\
                              process that you generate before running anvio. Please see the documentation for details.')
    parser.add_argument('--min-abundance', metavar='PERCENTAGE', default=0.1, type=float,
                        help='Short read-based taxonomy can be extremely noisy. Therefore, here we have defeault minimum\
                              percentage cutoff of 0.1%% to eliminate any taxon that occurs less than that in a given input\
                              file.')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
