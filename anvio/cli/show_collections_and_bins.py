#!/usr/bin/env python
# -*- coding: utf-8
"""Takes an anvi'o profile or pan database and returns a quick text summary
of bins and collections found in it."""

import sys

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
__requires__ = ["pan-db", "profile-db"]
__description__ = "A script to display collections stored in an anvi'o profile or pan database"


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

    utils.is_pan_or_profile_db(args.pan_or_profile_db)

    progress.new('Accessing to the collections table')
    progress.update('...')
    collections = ccollections.Collections()
    collections.populate_collections_dict(args.pan_or_profile_db)
    progress.end()

    if not collections.collections_dict:
        raise ConfigError("There are no collections in this profile database. Consider making some "
                          "in interactive mode, or importing a collection with "
                          "`anvi-import-collection` (a tutorial on importing collections here: "
                          "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-collection)")

    for source in collections.collections_dict:
        num_contigs = 0
        for bin_name, split_list in collections.get_collection_dict(source).items():
            contig_set = set([x.split("_split")[0] for x in split_list])
            num_contigs += len(contig_set)

        run.warning('', header = 'Collection: "%s"' % source, lc = 'green')
        run.info('Collection ID', source)
        run.info('Number of bins', collections.collections_dict[source]['num_bins'])
        run.info('Number of splits described', collections.collections_dict[source]['num_splits'])
        run.info('Number of contigs described', num_contigs)
        run.info('Bin names', ', '.join(sorted(collections.collections_dict[source]['bin_names'].split(','))), nl_after = 2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
