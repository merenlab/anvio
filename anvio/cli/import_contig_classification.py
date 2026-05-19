#!/usr/bin/env python

import sys

import anvio
import anvio.utils as utils
import anvio.tables as t
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.contigclassification import TablesForContigClassification, VALID_CLASSES

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = []
__requires__ = ['contigs-db', 'contig-classification-txt']
__provides__ = ['contig-classification']
__description__ = ("Import contig-level classification data (domain, virus, plasmid, etc.) "
                   "into a contigs database from a user-prepared tab-delimited file")


@terminal.time_program
def main():
    args = get_args()

    try:
        run_program(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program(args):
    A = lambda x: args.__dict__[x] if x in args.__dict__ else None

    contigs_db_path = A('contigs_db')
    input_files = A('input_files')
    just_do_it = A('just_do_it')

    combined_data = {}
    for input_file in input_files:
        file_data = utils.get_TAB_delimited_file_as_dictionary(input_file,
                                                               expected_fields=t.contig_classification_table_structure,
                                                               column_mapping=[str, int, str, str, str],
                                                               only_expected_fields=True,
                                                               indexing_field=-1)
        offset = len(combined_data)
        for k, v in file_data.items():
            combined_data[k + str(offset)] = v

    invalid_classes = set(row['class'] for row in combined_data.values()) - VALID_CLASSES
    if invalid_classes:
        raise ConfigError(f"The 'class' column contains invalid value(s): "
                          f"{', '.join(str(v) for v in sorted(invalid_classes))}. "
                          f"Valid class integers are: 0 (non-eukaryotic), 1 (eukaryotic), 2 (virus), "
                          f"3 (plasmid), 4 (organelle), 5 (unclassified).")

    tables = TablesForContigClassification(contigs_db_path)
    tables.create(combined_data, just_do_it=just_do_it)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('REQUIRED')
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    group1.add_argument('-i', '--input-files', metavar='FILE', nargs='+', default=None, required=True,
                        help="One or more tab-delimited files containing contig classification data. "
                             f"Each file must have the following columns: "
                             f"{', '.join(t.contig_classification_table_structure)}. "
                             "The 'class' column must be an integer: 0=non-eukaryotic, "
                             "1=eukaryotic, 2=virus, 3=plasmid, 4=organelle, 5=unclassified. "
                             "The 'source' column identifies the tool that produced the classification. "
                             "The 'tool_classification' column holds the raw string from the tool "
                             "(use semicolons for multi-level entries, e.g. Viruses;Duplodnaviria). "
                             "The 'confidence' column holds the tool confidence score, or 'NA' if "
                             "not available. Multiple files may be provided and will be concatenated.")

    group2 = parser.add_argument_group('OPTIONAL')
    group2.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
