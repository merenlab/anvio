#!/usr/bin/env python

import sys

import anvio
import anvio.phrogs as phrogs

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['copilot']
__requires__ = ['contigs-db', 'phrogs-data']
__provides__ = ['functions']
__description__ = "Run PHROGs HMMs on contigs-db"


@time_program
def main():
    args = get_args()

    try:
        p = phrogs.PHROGs(args)
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
    parser.add_argument('--phrogs-data-dir', default=None, type=str, help="The directory path for your downloaded PHROGs database. "
                        "If omitted, anvi'o tries the default location.")
    parser.add_argument('--noise-cutoff-terms', default='-E 1e-5', help="Filtering options for HMMER. By default, `-E 1e-5` is used.")
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('hmmer-program'), **anvio.K('hmmer-program'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
