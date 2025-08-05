#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.completeness import Completeness
from anvio.dbinfo import is_contigs_db


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db", "splits-txt", "hmm-source"]
__description__ = "A script to generate completeness info for a given list of _splits_"


def compute_completeness(args):
    run = terminal.Run()
    progress = terminal.Progress()

    is_contigs_db(args.contigs_db)

    completeness = Completeness(args.contigs_db, args.completeness_source)

    if args.list_completeness_sources:
        run.info('Available singlecopy sources', ', '.join(completeness.sources))
        sys.exit()

    contigs_db = dbops.ContigsDatabase(args.contigs_db)
    splits_in_db = set(contigs_db.db.get_table_as_dict(t.splits_info_table_name).keys())
    contigs_db.disconnect()

    if args.splits_of_interest:
        splits_in_users_list = set([s.strip() for s in open(args.splits_of_interest).readlines() if s.strip() and not s.startswith('#')])

        splits_of_interest = splits_in_db.intersection(splits_in_users_list)

        if len(splits_of_interest) != len(splits_in_users_list):
            if not len(splits_of_interest):
                run.warning('None of the split names you provided in %s matched split names in the database...' % args.splits_of_interest)
                sys.exit()
            else:
                run.warning('Only %d of %d split names you listed in "%s" matched split names in the database...'\
                                                % (len(splits_of_interest), len(splits_in_users_list), args.splits_of_interest))
    else:
        splits_of_interest = splits_in_db

    try:
        p_completion, p_redundancy, domain, domain_probabilities, info_text, results_dict = completeness.get_info_for_splits(splits_of_interest, min_e_value = args.min_e_value)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

    results = [list(v.values())[0] for v in results_dict.values()]

    run.warning('', header = 'Completeness for %d splits (p < %g)' % (len(splits_of_interest), args.min_e_value))
    for v in results:
        run.info('%s (%sl SCGs)' % (v['source'], v['domain']), '%.2f%% complete, %.2f%% redundant' % (v['percent_completion'], v['percent_redundancy']))

    print()


@terminal.time_program
def main():
    args = get_args()

    try:
        compute_completeness(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('splits-of-interest'), **anvio.K('splits-of-interest'))
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('min-e-value'), **anvio.K('min-e-value'))
    parser.add_argument(*anvio.A('list-completeness-sources'), **anvio.K('list-completeness-sources'))
    parser.add_argument(*anvio.A('completeness-source'), **anvio.K('completeness-source'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
