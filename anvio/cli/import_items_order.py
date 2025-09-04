#!/usr/bin/env python
# -*- coding: utf-8
"""A program to import an items order into an anvi'o database"""

import sys

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbinfo import is_blank_profile, is_pan_or_profile_db
from anvio.utils.database import get_db_type


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db', 'profile-db', 'misc-data-items-order-txt', 'dendrogram', 'phylogeny']
__provides__ = ['misc-data-items-order',]
__description__ = "Import a new items order into an anvi'o database"


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
    progress = terminal.Progress()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    input_file_path = A('input_order')
    db_path = A('db_path')
    order_name = A('name') or 'LAZY_USERS_UNKNOWN_ITEMS_ORDER'
    make_default = A('make_default')

    if not input_file_path or not db_path:
        raise ConfigError("Probably it will come as a surprise, but you *must* provide a file path for the "
                          "input order and the target database :/")

    filesnpaths.is_file_plain_text(input_file_path)
    is_pan_or_profile_db(db_path, genes_db_is_also_accepted=True)

    order_data = [l.strip() for l in open(input_file_path).readlines() if len(l.strip())]

    # Determine the nature of the incoming order data:
    if len(order_data) == 1:
        try:
            filesnpaths.is_proper_newick(order_data[0])
        except Exception as e:
            raise ConfigError("Your input file contained a single line and anvi'o thought that it shold be a "
                              "newick tree, but then ETE complained about it, so we are all upset here: '%s' "
                              ":(" % (e))

        order_data = order_data[0]
        order_data_type_newick = True
    else:
        order_data_type_newick = False
        order_data = ','.join(order_data)

    db_type = get_db_type(db_path)

    run.info("Target database", db_path)
    run.info("Database type", db_type, nl_after=1)
    run.info("Order file path", input_file_path)
    run.info("Order data type", 'newick' if order_data_type_newick else 'basic')
    run.info("Order name", order_name, mc='red', nl_after=1)

    # for blank profile databases, we don't ask for for the consistency of names
    # to be checked.
    if db_type == 'profile' and is_blank_profile(db_path):
        check_names_consistency = False
    else:
        check_names_consistency = True

    if order_data_type_newick:
        dbops.add_items_order_to_db(db_path, order_name, order_data, order_data_type_newick=order_data_type_newick,
                                    distance='NA', linkage='NA', dont_overwrite=True, make_default=make_default,
                                    check_names_consistency=check_names_consistency)
    else:
        dbops.add_items_order_to_db(db_path, order_name, order_data, order_data_type_newick=order_data_type_newick,
                                    dont_overwrite=True, make_default=make_default, check_names_consistency=check_names_consistency)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupI = parser.add_argument_group('CRITICAL INPUT', 'Basically the input file and the target database')
    groupI.add_argument('-i', '--input-order', metavar = 'FILE', help = "One of the two important things you must provide: the file\
                                that contains the items order. The format of this file is important. It can either contain a\
                                proper newick tree in it, or a complete list of 'items' in the target database where every line\
                                of the file is simply an item name. If you are providing a newick tree, the entire file should\
                                be a single line. I know it sounds hard, but you seriously can do this.")
    groupI.add_argument('-p', '--db-path', metavar = 'DB PATH', help = "An appropriate anvi'o database to import the items order. Currently\
                                it can be a profile, pan, or genes database. But you should try your chances with other kinds of\
                                databases for fun and games. Basically, if the database contains an items order table, then things\
                                will work. Otherwise, you will probably get angry errors back in the worst case scenario.")

    groupZ = parser.add_argument_group('NOT SO CRITICAL INPUT', 'Because not all parameters are created equal')
    groupZ.add_argument('--name', metavar='ORDER NAME', help="What should we call this order? Give it a concise, single-word name.")
    groupZ.add_argument('--make-default', default=False, action="store_true", help="You have the option to make this order the \
                                default order in the database. Which means, anvi'o will use this one when someone runs the\
                                program anvi-interactive and presses draw. Big responsibility. But if you have a 'default' state,\
                                it will not work because the default items order in the state file overwrites the one that comes\
                                from the database. So not that big of a responsibility.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
