#!/usr/bin/env python
# -*- coding: utf-8
DESCRIPTION = """This program imports protein abundance data into a profile database."""

import pandas as pd

from sys import exit
from argparse import Namespace

import anvio.tables as tables

from anvio import A, K
from anvio.errors import ConfigError
from anvio.dbops import ProfileDatabase
from anvio import __version__ as VERSION
from anvio.argparse import ArgumentParser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__authors__ = ['semiller10']
__requires__ = ['profile-db']
__provides__ = []
__description__ = DESCRIPTION


def main() -> None:
    args = get_args()

    input_table = pd.read_csv(args.protein_profile_table, sep='\t', header=0)
    column_names = ['source', 'accession', 'sample', 'abundance']

    try:
        if sorted(input_table.columns) != sorted(column_names):
            raise ConfigError(
                "The protein profile table did not have the required header with the following column "
                f"names: {', '.join(column_names)}."
            )
        input_table = input_table.rename(
            {
                'source': 'reference_source',
                'accession': 'reference_id',
                'sample': 'sample_name',
                'abundance': 'abundance_value'
            },
            axis=1
        )

        # Assign unique protein IDs.
        protein_id = 0
        protein_id_table = []
        for source, accession in set(zip(input_table['reference_source'], input_table['reference_id'])):
            protein_id_table.append([protein_id, source, accession])
            protein_id += 1
        input_table = input_table.merge(
            pd.DataFrame(protein_id_table, columns=['protein_id', 'reference_source', 'reference_id']),
            on=['reference_source', 'reference_id'],
            how='left'
        )
        input_table = input_table[tables.protein_abundances_table_structure]

        pdb = ProfileDatabase(args.profile_db)
        pdb.db._exec(f'''DELETE from {tables.protein_abundances_table_name}''')
        pdb.db._exec_many(
            f'''INSERT INTO {tables.protein_abundances_table_name} VALUES
            ({','.join('?' * len(tables.protein_abundances_table_structure))})''',
            input_table.values
        )
        pdb.disconnect()
    except ConfigError as e:
        print(e)
        exit(-1)


def get_args() -> Namespace:
    parser = ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        '--protein-profile-table', required=True,
        help=
        """This tab-delimited table contains protein abundance data from different samples. It \
        should have the following named columns: 'source', 'accession', 'sample', and 'abundance'. \
        Each row corresponds to a distinct protein abundance measurement. 'source' is the source of \
        the protein accessions, and must be a gene function annotation source stored in the anvi'o \
        profile database (available sources can be found with the program, `anvi-db-info`). \
        'accession' is the protein ID in the annotation source. A contigs database built from a \
        GenBank file, for example, could contain the source, 'NCBI_PGAP', and the accession, \
        'WP_011862028.1'. 'sample' is the name of the sample in which the measurement was made. It \
        need not be the same as any nucleotide sequence samples stored in the profile database. \
        'abundance' is the protein abundance value itself."""
    )
    parser.add_argument(*A('profile-db'), **K('profile-db'))

    args = parser.get_args(parser)
    return args


if __name__ == '__main__':
    main()
