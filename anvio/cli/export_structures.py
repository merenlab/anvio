#!/usr/bin/env python
# -*- coding: utf-8

import sys
import anvio
from anvio.argparse import ArgumentParser

import anvio.utils as utils
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.structureops import StructureDatabase


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__requires__ = ['structure-db']
__provides__ = ['protein-structure-txt']
__description__ = "Export .pdb structure files from a structure database"


def main():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    gene_caller_ids = A('gene_caller_ids')
    genes_of_interest_path = A('genes_of_interest')
    structure_db_path = A('structure_db')
    output_dir = A('output_dir')

    try:
        utils.is_structure_db(structure_db_path)
        structure_db = StructureDatabase(structure_db_path, ignore_hash=True, create_new=False)
        filesnpaths.check_output_directory(output_dir)

        if gene_caller_ids and genes_of_interest_path:
            raise ConfigError("Pick one of --gene-caller-ids and --genes-of-interest")
        elif genes_of_interest_path:
            filesnpaths.is_file_exists(args.genes_of_interest)
            genes_of_interest = set(int(g.strip()) for g in open(args.genes_of_interest, 'r').readlines())
        elif gene_caller_ids:
            genes_of_interest = set(int(g) for g in gene_caller_ids.split(','))
        else:
            genes_of_interest = structure_db.genes_with_structure

        genes_missing_from_structure_db = [gene for gene in genes_of_interest if gene not in structure_db.genes_with_structure]
        if genes_missing_from_structure_db:
            show_a_few = genes_missing_from_structure_db if len(genes_missing_from_structure_db) <= 10 else genes_missing_from_structure_db[:10]
            raise ConfigError("{} gene(s) were specified by you but don't exist in the structure database. Here are some of their IDs: {}".
                format(len(genes_missing_from_structure_db), ', '.join([str(x) for x in show_a_few])))

        structure_db.export_pdbs(genes_of_interest, output_dir)
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db'))
    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    parser.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    parser.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))

    return parser.parse_args()


if __name__ == '__main__':
    main()
