#!/usr/bin/env python
# -*- coding: utf-8
"""A program to export an items order from an anvi'o database"""

import sys

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db', 'profile-db']
__provides__ = ['misc-data-items-order-txt', 'dendrogram', 'phylogeny']
__description__ = "Export an item order from an anvi'o database"


def main():
    args = get_args()
    run = terminal.Run()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    order_name = A('name')
    db_path = A('db_path')

    output_file_path = A('output_file') or 'unknown_items_order.txt'

    try:
        if not db_path:
            raise ConfigError("Probably it will come as a surprise, but you *must* provide an input database path :/")

        filesnpaths.is_output_file_writable(output_file_path)

        item_order_names, item_orders_dict = dbops.get_item_orders_from_db(db_path)

        if not len(item_order_names):
            raise ConfigError("There are no item orders in this database :/")

        if not order_name:
            run.warning("You must choose an order. Here is what you have in here:", header="Available item orders", lc='yellow')

            for item_order in item_order_names:
                item_order_name, item_order_distance, item_order_clustering = item_order.split(':')
                nl_after = 1 if item_order == item_order_names[-1] else 0
                if item_order_distance:
                    run.info_single("%s (newick; distance: %s, clustering: %s)." % (item_order_name, item_order_distance, item_order_clustering), nl_after=nl_after)
                else:
                    run.info_single("%s (list)." % (item_order_name), nl_after=nl_after)

        items_order_of_interest = None
        for item_order in item_order_names:
            item_order_name, item_order_distance, item_order_clustering = item_order.split(':')
            if order_name == item_order_name or order_name == item_order:
                items_order_of_interest = item_orders_dict[item_order]

        if not items_order_of_interest:
            raise ConfigError(f"The item order '{order_name}' is not one of the item orders in the database. Here "
                              f"is a comma-separated list of what you have in there (please note that you can, but do not need "
                              f"to include the distance/clustering types after the ':' character in each name): "
                              f"{', '.join(item_order_names)}")

        order_data_type_newick = items_order_of_interest['type'] == 'newick'
        run.info("Database", db_path)
        run.info("Database type", utils.get_db_type(db_path))
        run.info("Order name", order_name)
        run.info("Order data type", 'newick' if order_data_type_newick else 'basic')

        if order_data_type_newick:
            open(output_file_path, 'w').write('%s\n' % items_order_of_interest['data'])
        else:
            open(output_file_path, 'w').write('%s\n' % '\n'.join(items_order_of_interest['data']))

        run.info("Output file", output_file_path, mc='red')

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupI = parser.add_argument_group('INPUT', 'The database and the items order of interest')
    groupI.add_argument('-p', '--db-path', metavar = 'DB PATH', help = "An appropriate anvi'o database.")
    groupI.add_argument('--name', metavar='ORDER NAME', help="The name of the order you want to export. If you don't\
                                    provide an order name, anvi'o will show you the names of all available orders in\
                                    the database.")

    groupZ = parser.add_argument_group('OUPUT', 'Output file name and stuff')
    groupZ.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
