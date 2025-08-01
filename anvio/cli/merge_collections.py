#!/usr/bin/env python
# -*- coding: utf-8

"""A program to merge multiple binning results into a single additional data file"""

import os
import sys

import anvio
import anvio.utils as u
import anvio.terminal as terminal
from anvio.errors import ConfigError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__authors__ = ['meren']
__provides__ = []
__requires__ = ["contigs-db", "collection-txt"]
__description__ = "Generate an additional data file from multiple collections"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(1)


def run_program():
    args = get_args()
    run = terminal.Run()

    input_dicts = {}
    for input_file in args.input_files:
        if len(u.get_columns_of_TAB_delim_file(input_file)) != 1:
            raise ConfigError("Are you sure %s is an anvi'o collections file? "
                               "Because it doesn't look like one :/" % input_file)

        source = '.'.join(os.path.basename(os.path.abspath(input_file)).split('.')[:-1])

        input_dicts[source] = u.get_TAB_delimited_file_as_dictionary(input_file, no_header = True)

        run.info('New Source', '%s, w/ %d contigs' % (source, len(input_dicts[source])))

    contig_names = []
    for source in input_dicts:
        contig_names.extend(list(input_dicts[source].keys()))
    contig_names = set(contig_names)
    run.info('Final number of unique contigs', len(contig_names))

    sources = sorted(input_dicts.keys())

    contig_name_to_splits = u.get_contig_name_to_splits_dict(args.contigs_db)

    combined_dict = {}
    for contig_name in contig_names:
        if not contig_name in contig_name_to_splits:
            raise ConfigError("Oh. You have the wrong stuff. Probably. Because, the contig '%s' does not match to any of "
                               "the contig names in your database. Here is a random contig name you have in it in "
                               "comparison: '%s'." % (contig_name, list(contig_name_to_splits.keys())[0]))
        for split_name in contig_name_to_splits[contig_name]:
            combined_dict[split_name] = dict(list(zip(sources, [None] * len(sources))))

        for source in sources:
            if contig_name in input_dicts[source]:
                for split_name in contig_name_to_splits[contig_name]:
                    combined_dict[split_name][source] = input_dicts[source][contig_name]['column_00001']


    u.store_dict_as_TAB_delimited_file(combined_dict, args.output_file, headers = ['contig'] + sources)
    run.info('Output', args.output_file)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument('-i', '--input-files', metavar = 'FILE(S)', nargs='+', default = None, required = True,
                        help = 'Input file(s). TAB-delimited input files should have two columns, where the first\
                                               column holds the contig name, and the second one the bin id. This\
                                               is the standard ouptut of the program anvi-export-collection.')
    parser.add_argument('-o', '--output-file', metavar = 'OUTPUT_FILE', default = None, required = True,
                        help = 'Output file name.')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
