#!/usr/bin/env python
# -*- coding: utf-8
"""A script to import collections (and their colors)"""

import sys
import copy
from collections import Counter

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.collections import TablesForCollections


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'profile-db', 'pan-db', 'collection-txt',]
__provides__ = ['collection',]
__description__ = "Import an external binning result into anvi'o"
__resources__ = [("Another description as part of the metagenomic workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-collection")]




@terminal.time_program
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

    #################################################################################################
    # <SANITY CHECKS>
    #################################################################################################
    if not args.pan_or_profile_db:
        raise ConfigError("You must provide an anvi'o pan or profile database for this to work :(")

    utils.is_pan_or_profile_db(args.pan_or_profile_db)

    filesnpaths.is_output_file_writable(args.pan_or_profile_db)

    if args.contigs_db:
        utils.is_contigs_db(args.contigs_db)

    if args.pan_or_profile_db and utils.get_db_type(args.pan_or_profile_db) == 'pan' and args.contigs_db:
        raise ConfigError("There is no need to provide a contigs database when you are working with an anvi'o pan "
                           "database")

    if not args.contigs_db and args.contigs_mode:
        raise ConfigError("There is no reason for you to use the `--contigs-mode` flag when you have "
                           "not declared an anvi'o contigs database")

    # if there is a profile database, check whether there is a contigs database associated with the profile
    if args.pan_or_profile_db and utils.get_db_type(args.pan_or_profile_db) == 'profile':
        pan_or_profile_db = dbops.ProfileDatabase(args.pan_or_profile_db)

        if pan_or_profile_db.meta['contigs_db_hash'] and not args.contigs_db:
            raise ConfigError("The profile database you provided is associated with an anvi'o contigs database (i.e., "
                              "it is not an ad hoc profile database but a part of the cool dbs club). Which means, you "
                              "must provide a path for the contigs database, so anvi'o can double check your item names "
                              "in your collection. Sorry :/")

    if not args.collection_name:
        raise ConfigError("You must give a name for this collection.")

    if not args.contigs_db:
        run.warning("You did not provide a contigs database. Fine. So be it. But know this: anvi'o has no way to check "
                    "the consistency of names you provide in the input file. So if you made a mistake while generating "
                    "this collection, it probably will cause issues later on.")

    if not args.pan_or_profile_db:
        run.warning("Since you haven't provided an anvi'o profile database, this program will add your collection into "
                    "the contigs database you provided. If you use the same collection name later in one of your profile "
                    "databases that will be generated from this contigs database, things may go South, and anvi'o would "
                    "not even care.")

    utils.check_collection_name(args.collection_name)

    filesnpaths.is_file_tab_delimited(args.data, expected_number_of_fields = 2)
    if args.bins_info:
        filesnpaths.is_file_tab_delimited(args.bins_info, expected_number_of_fields = 3)

    num_occurences_of_entries = Counter([l.split('\t')[0] for l in open(args.data).readlines()])
    if max(num_occurences_of_entries.values()) != 1:
        raise ConfigError("Some %(item)s names occur more than once in the input file. A %(item)s cannot belong in two "
                           "bins, and neither there should be the same bin assignment for a given %(item)s. Long story "
                           "short, each name should appear only once in your input file, and it is not the case :/" \
                                                                        % {'item': 'contig' if args.contigs_mode else 'split'})
    #################################################################################################
    # </SANITY CHECKS>
    #################################################################################################

    # initiate the contigs database if it is present
    if args.contigs_db:
        contig_name_to_splits_dict = utils.get_contig_name_to_splits_dict(args.contigs_db)
        contig_names_in_db = contig_name_to_splits_dict.keys()

    # read the input file with split/contig - bin ID associations
    input_data_file_content = utils.get_TAB_delimited_file_as_dictionary(args.data, no_header = True, column_names = ['split_id', 'bin_name'])

    # populate bins_info_dict there is any information about bins
    bins_info_dict = {}
    if args.bins_info:
        try:
            bins_info_dict = utils.get_TAB_delimited_file_as_dictionary(args.bins_info, no_header = True, column_names = ['bin_name', 'source', 'html_color'])
        except Exception as e:
            raise ConfigError("Someone was not happy with the TAB-delimited bins info file you provided. Here "
                              "is the complaint: %s" % e)

    # convert contig names into split names if necessary so we can do everything else using the information in
    # pan or profile database item names (this is only relevant for a workflow that includes a contigs db)
    if args.contigs_mode:
        input_names_missing_from_contigs_db = [n for n in input_data_file_content if n not in contig_names_in_db]

        if len(input_names_missing_from_contigs_db):
            if len(input_names_missing_from_contigs_db) == len(input_data_file_content):
                raise ConfigError(f"None of the contig names in your input file matches to the contig names found in the contigs :( "
                                  f"Maybe those are split names and you are mistakenly using `--contigs-mode` flag? Well, here is an "
                                  f"example: {input_names_missing_from_contigs_db[0]}.")
            else:
                raise ConfigError(f"OK. Of the {len(input_data_file_content)} in your input file, {len(input_names_missing_from_contigs_db)} "
                                  f"have no match to any contig names found in your contigs database. Here is an example contig name from "
                                  f"your file that does not match to anything in the contigs db: {input_names_missing_from_contigs_db[0]}")

        # convert input data names to split names:
        contig_names = list(input_data_file_content.keys())
        for contig_name in contig_names:
            for split_name in contig_name_to_splits_dict[contig_name]:
                input_data_file_content[split_name] = copy.deepcopy(input_data_file_content[contig_name])
            input_data_file_content.pop(contig_name)

        run.info('Contig/split name conversion', '%d contig names converted into %d split names.' % (len(contig_names), len(input_data_file_content)))

    run.info('Item names in input', len(input_data_file_content))
    run.info('Num bins in input', len(set([e['bin_name'] for e in list(input_data_file_content.values())])))

    # learning about the input names like a pro
    input_names = set(input_data_file_content.keys())

    # here we attempt to make sure the names in the input file are relevant to the
    # names in the contigs database database. but clearly it is not relevant if there
    # is no contigs database is associated with the profile database, so if there is none,
    # we cheat :
    if args.pan_or_profile_db:
        if utils.get_db_type(args.pan_or_profile_db) == 'profile':
            associated_with_a_contigs_db = dbops.ProfileDatabase(args.pan_or_profile_db).meta['contigs_db_hash']
        else:
            associated_with_a_contigs_db = None

        if utils.is_blank_profile(args.pan_or_profile_db):
            if associated_with_a_contigs_db:
                db_names = utils.get_all_item_names_from_the_database(args.contigs_db)
            else:
                db_names = input_names

                run.warning("Since you are working with a blank proifle, anvi'o is not going to check whether the names item "
                            "names in your collections file matches to the item names in other databases of yours. It is all "
                            "fine for now, but this requires you to be even extra careful with your downstream analyses in "
                            "case stuff hits the fan later.")
        else:
            db_names = utils.get_all_item_names_from_the_database(args.pan_or_profile_db)

            if not len(db_names):
                raise ConfigError("Anvi'o got an empty list of names for `items` from this databaes. But it is impossible :/")

        run.info('Items in %s database' % utils.get_db_type(args.pan_or_profile_db), len(db_names), mc='green')
    elif args.contigs_db:
        db_names = utils.get_all_item_names_from_the_database(args.contigs_db)
        run.info('Items in contigs database', len(db_names), mc='green')
    else:
        db_names = input_names
        run.warning("Quite an improper setup (no pan, profile, or contigs databases). We are using %d "
                    "input names found in the collection file as item names. Because you are the "
                    "boss that's why." % len(db_names))

    item_names_shared = set.intersection(db_names, input_names)
    if not len(item_names_shared):
        raise ConfigError(f"There is no overlap between the item names found in your input file and the item names found "
                          f"in the database :( For instance, one of the names in your file looks like this: "
                          f"'{input_names.pop()}', and in contrast this is an example name from the database: '{db_names.pop()}'. "
                          f"Solving this may be as easy as adding the flag `--contigs-mode` to your command. Or it may also be "
                          f"the case that the content of the input file has nothing to do with the database you are trying "
                          f"to import these items into. Anvi'o is Jon Snow and hopes that you will figure this out and come "
                          f"back.")
    else:
        run.info("Item names shared between input and db", len(item_names_shared))

    # find names that are unique to the input file
    item_names_unique_to_input_file = input_names - db_names
    if len(item_names_unique_to_input_file):
        for item_name in item_names_unique_to_input_file:
            input_data_file_content.pop(item_name)

        input_names = set(input_data_file_content.keys())
        run.warning('%d item(s) that appeared only in the input file will be ignored (such as this one: %s). '
                    'Just so you know.' % (len(item_names_unique_to_input_file), item_names_unique_to_input_file.pop()))

    # items in db but missing from the input file:
    item_names_unique_to_db = db_names - input_names
    if len(item_names_unique_to_db):
        run.warning("%d item(s) that were in the database, but were not in the input file, will not be described by "
                    "any bin in the collection %s. That is totally fine, but anvi'o hopes that you are aware of that. "
                    "This means you have more things in your database than the number of things your input file "
                    "describes. Here is an example of something that is in your database but not in any bin in your "
                    "input file: %s." % (len(item_names_unique_to_db), args.collection_name, item_names_unique_to_db.pop()))

    data = {}

    # populate the data dictionary
    for entry_name in input_data_file_content:
        bin_name = input_data_file_content[entry_name]['bin_name']

        if bin_name not in data:
            data[bin_name] = set([])

        data[bin_name].add(entry_name)

    if args.pan_or_profile_db:
        collections = TablesForCollections(args.pan_or_profile_db)
    else:
        collections = TablesForCollections(args.contigs_db)

    collections.append(args.collection_name, data, bins_info_dict)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('data', metavar = "TAB DELIMITED FILE",
                        help = 'The input file that describes bin IDs for each split or contig.')

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    parser.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db', {'required': False}))
    parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name', {'required': True}))
    parser.add_argument(*anvio.A('bins-info'), **anvio.K('bins-info'))
    parser.add_argument(*anvio.A('contigs-mode'), **anvio.K('contigs-mode'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
