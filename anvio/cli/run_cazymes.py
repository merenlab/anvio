#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.cazymes as cazymes

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mschecht']
__requires__ = ['contigs-db', 'cazyme-data']
__provides__ = ['functions']
__description__ = "Run dbCAN CAZymes on contigs-db"


@time_program
def main():
    args = get_args()

    try:
        p = cazymes.CAZyme(args)
        p.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument('--noise-cutoff-terms', default=None, help="Choose a threshold for the CAZyme HMMs (default --noise-cutoff-terms \"-E 1e-12\". \
                        If you are not sure what to pick run `hmmsearch -h` or `hmmscan -h` and check out your options.")
    parser.add_argument('--cazyme-data-dir', default=None, type=str, help="The directory path for your downloaded CAZyme database. Anvi'o will try to use the default path \
                        if you do not specify anything.")
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('hmmer-program'), **anvio.K('hmmer-program'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
