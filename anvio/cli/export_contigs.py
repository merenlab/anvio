#!/usr/bin/env python
# -*- coding: utf-8
"""A script to export a FASTA file of contigs (or splits) from a contigs database."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
import anvio.errors


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db']
__provides__ = ['contigs-fasta']
__description__ = "Export contigs (or splits) from an anvi'o contigs database"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        if args.contigs_of_interest:
            filesnpaths.is_file_tab_delimited(args.contigs_of_interest, expected_number_of_fields=1)
            seq_names_to_export = [line.strip() for line in open(args.contigs_of_interest).readlines()]
        else:
            seq_names_to_export = None

        utils.export_sequences_from_contigs_db(args.contigs_db,
                                               args.output_file,
                                               seq_names_to_export=seq_names_to_export,
                                               splits_mode=args.splits_mode,
                                               just_do_it=args.just_do_it,
                                               truncate=(not args.no_wrap),
                                               run=run)

        run.info('Export mode', 'splits' if args.splits_mode else 'contigs')
        run.info('Output FASTA', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('contigs-of-interest'), **anvio.K('contigs-of-interest'))
    parser.add_argument('--splits-mode', default=False, action="store_true", help="Export split\
                        sequences instead.")
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    parser.add_argument(*anvio.A('no-wrap'), **anvio.K('no-wrap'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
