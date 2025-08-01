#!/usr/bin/env python
# -*- coding: utf-8
"""Allows you to merge multiple bins in a collection"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.ccollections as ccollections

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["pan-db", "profile-db", "collection", "bin"]
__provides = ["bin"]
__description__ = "Merge a given set of bins in an anvi'o collection"


run = terminal.Run()
progress = terminal.Progress()


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def run_program():
    args = get_args()

    utils.is_pan_or_profile_db(args.pan_or_profile_db)

    c = ccollections.Collections(r = run, p = progress)
    c.populate_collections_dict(args.pan_or_profile_db)

    if args.list_collections:
        c.list_collections()
        sys.exit()

    if not args.collection_name:
        raise ConfigError("You must declare a colection name. Use --list-collections flag "
                           "to list available ones.")

    c.sanity_check(args.collection_name)

    collection_info = c.collections_dict[args.collection_name]
    run.info('Collection found', '"%s" (%d items in %d bins)' % (args.collection_name, collection_info['num_splits'], collection_info['num_bins']))

    if args.list_bins:
        c.list_bins_in_collection(args.collection_name)
        sys.exit()

    if not args.bin_names_list:
        raise ConfigError("Anvi'o needs a list of bin names to merge :/")

    if not args.new_bin_name:
        raise ConfigError("Anvi'o needs a new bin name. Or it will make up one for you, and you "
                          "know what happens anvi'o makes stuff up.")

    bin_names_list = [b.strip() for b in args.bin_names_list.split(',')]
    run.info('Bins to merge', ', '.join(bin_names_list))

    utils.is_this_name_OK_for_database('bin name', args.new_bin_name, stringent=False)
    run.info('New bin name', args.new_bin_name)

    c.merge_bins(args.collection_name, args.new_bin_name, bin_names_list)

    if not len(bin_names_list):
        raise ConfigError("Your bin names list variable looks empty. Something is wrong..")


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("DB AND COLLECTION", "Simple enough. This guy needs a pan or profile database\
                                        and a collection name. You can get a list of available collections with another\
                                        flag down below.")
    groupA.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))

    groupB = parser.add_argument_group("BINS TO WORK WITH", "Here you need to define a list of bin names to merge, and the\
                                        new bin name for them to merge under. Your bin names should be comma-separated.\
                                        Both 'name_1, name_2, name_3' and name_1,name_2,name_3 will work. Your new bin name\
                                        better be a single word, meaningful name so anvi'o does not complain about it later.")
    groupB.add_argument(*anvio.A('bin-names-list'), **anvio.K('bin-names-list'))
    groupB.add_argument(*anvio.A('new-bin-name'), **anvio.K('new-bin-name'))

    groupC = parser.add_argument_group("SWEET FLAGS OF CONVENIENCE", "We gotchu.")
    groupC.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupC.add_argument(*anvio.A('list-bins'), **anvio.K('list-bins'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
