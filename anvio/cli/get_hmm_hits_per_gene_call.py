#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db", "hmm-source", "hmm-hits"]
__provides__ = ["functions-txt"]
__description__ = ("A simple script to generate a TAB-delimited file gene caller IDs and their "
                   "HMM hits for a given HMM source")


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
    args, unknown = get_args()
    run = terminal.Run()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    hmm_source = A('hmm_source')
    contigs_db_path = A('contigs_db')
    output_file_path = A('output_file')

    utils.is_contigs_db(contigs_db_path)
    filesnpaths.is_output_file_writable(output_file_path)

    if not hmm_source:
        raise ConfigError("You need to declare an HMM source for this to work :/")

    if not output_file_path:
        raise ConfigError("This will not work without an output file path.")

    contigs_db = dbops.ContigsDatabase(contigs_db_path)
    hmm_results_dict = utils.get_filtered_dict(contigs_db.db.get_table_as_dict(t.hmm_hits_table_name), 'source', set([hmm_source]))

    if not len(hmm_results_dict):
        raise ConfigError("Your contigs database does not have any HMM hits for the HMM source %s :/" % hmm_source)
    else:
        run.info('Num HMM hits', len(hmm_results_dict))

    header = ['unique_entry_id', 'gene_callers_id', 'source', 'gene_name', 'gene_hmm_id']
    utils.store_dict_as_TAB_delimited_file(hmm_results_dict, output_file_path, headers=header)

    run.info('Output', output_file_path)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupB = parser.add_argument_group("INPUT: ANVI'O CONTIGS DB")
    groupB.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupD = parser.add_argument_group('INPUT: HMM SOURCE')
    groupD.add_argument(*anvio.A('hmm-source'), **anvio.K('hmm-source'))

    groupD = parser.add_argument_group('OUTPUTTAH')
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True }))

    return parser.parse_known_args()


if __name__ == "__main__":
    main()
