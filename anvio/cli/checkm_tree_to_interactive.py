#!/usr/bin/env python
# -*- coding: utf-8
import os
import re
import sys

from ete3 import Tree

import anvio
import anvio.utils as utils
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ozcan']
__requires__ = ["phylogeny"]
__provides__ = ["interactive"]
__description__ = ("A helper script to convert CheckM trees into anvio interactive with taxonomy information")


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
    input_tree = args.tree
    output_dir = args.output_dir

    run = terminal.Run()

    taxon_letter_to_name = dict((column[2], column[2:]) for column in tables.taxon_names_table_structure[1:])
    filesnpaths.is_file_exists(input_tree)
    filesnpaths.check_output_directory(output_dir)

    with open(input_tree, 'r') as f:
        newick = f.read()

    run.info('Input tree found', input_tree)

    uid_regex = re.compile(r"'?(\w+)\|(.*?)\|([\.eE\-0-9]+)?'?")
    uid_to_taxonomy = {}

    for match in re.finditer(uid_regex, newick):
        uid = match.group(1)
        taxonomy_list = match.group(2).split(';')

        for taxonomy_item in taxonomy_list:
            if len(taxonomy_item):
                # format of taxonomy_item is:
                # o__Planctomycetales
                taxon_letter, full_name = taxonomy_item[0], taxonomy_item[3:]

                if taxon_letter not in taxon_letter_to_name:
                    continue

                taxon = taxon_letter_to_name[taxon_letter]

                if uid not in uid_to_taxonomy:
                    uid_to_taxonomy[uid] = {}

                uid_to_taxonomy[uid][taxon] = full_name

    newick = re.sub(uid_regex, r'\1', newick)
    tree = Tree(newick, format=1)

    output_dict = {}
    for leaf in tree.get_leaves():
        output_dict[leaf.name] = dict((taxon, '') for taxon in taxon_letter_to_name.values())

    for node in tree.traverse("levelorder"):
        if node.name in uid_to_taxonomy:
            for leaf in node.get_leaves():
                for taxon in uid_to_taxonomy[node.name]:
                    output_dict[leaf.name][taxon] = uid_to_taxonomy[node.name][taxon]

    os.mkdir(output_dir)
    output_newick_path = os.path.join(output_dir, 'newick.tree')
    output_view_data_path = os.path.join(output_dir, 'view_data.txt')

    run.info('Output tree path', output_newick_path)
    run.info('Output view data path', output_view_data_path)

    utils.store_dict_as_TAB_delimited_file(output_dict, output_view_data_path,
        headers=['items'] + list(taxon_letter_to_name.values()))

    with open(output_newick_path, 'w') as f:
        f.write(newick)

    run.info_single("Successfully converted.", nl_before=1, nl_after=1)


def get_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument('-t', '--tree', required=True, metavar='CHECKM TREE',
                            help="Tree file generated by CheckM.")
    parser.add_argument('-o', '--output-dir', required=True, metavar='DIRECTORY', 
                            help="The directory name that output files will be stored.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
