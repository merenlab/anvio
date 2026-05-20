#!/usr/bin/env python

import sys

import anvio

from anvio.user_annotation import UserAnnotationDBSetup
from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program
from anvio import constants


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['l-gallucci', 'iwilkie']
__provides__ = ['user-annotation-db']
__description__ = ("Set up user-provided HMM profiles or protein FASTA files as custom functional "
                   "annotation databases for use with `anvi-run-user-annotation`")


@time_program
def main():
    args = get_args()

    try:
        setup = UserAnnotationDBSetup(args)
        setup.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("INPUT", "A tab-delimited file that lists the databases to set up.")
    groupA.add_argument('--input-tsv', metavar='FILE', default=None,
                        help="A tab-delimited text file with two or three columns: 'name', 'path', and "
                             "optionally 'companion_fasta'. Each row describes one database. 'name' is a "
                             "unique identifier, 'path' is the HMM profile or protein FASTA file, and the "
                             "optional 'companion_fasta' is a protein FASTA of the sequences that built "
                             "the HMM profiles (enables cross-validation with DIAMOND when running). "
                             "HMM profiles (.hmm or .hmm.gz) and protein FASTA files (.faa, .fasta, .fa, "
                             ".fas, .fna) are both accepted. Lines starting with '#' are treated as comments.")

    groupB = parser.add_argument_group("OUTPUT", "Where to store the prepared databases.")
    groupB.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {
        'required': False,
        'default': constants.default_user_annotation_data_dir,
        'help': f"Directory where prepared databases will be stored. "
                f"Default: %(default)s"
    }))

    groupC = parser.add_argument_group("OPTIONS")
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupC.add_argument(*anvio.A('reset'), **anvio.K('reset',
                        {'help': "Remove the existing output directory and start fresh. Use with caution: "
                                 "all previously prepared databases in the directory will be lost. "
                                 "Cannot be combined with --list or --remove."}))

    groupD = parser.add_argument_group("MANIFEST MANAGEMENT",
                                       "Inspect or modify the manifest of an existing annotation directory "
                                       "without re-running the full setup. Use these flags STANDALONE — "
                                       "they cannot be combined with --input-tsv or --reset.")
    groupD.add_argument('--list', default=False, action='store_true',
                        help="List all databases currently registered in the annotation directory and exit. "
                             "Standalone flag: do not combine with --input-tsv or --reset.")
    groupD.add_argument('--remove', default=None, metavar='NAME',
                        help="Remove a database from the manifest and delete its prepared files. "
                             "Pass the database name as it appears in the manifest. "
                             "The other databases in the directory are not affected. "
                             "Standalone flag: do not combine with --input-tsv or --reset.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
