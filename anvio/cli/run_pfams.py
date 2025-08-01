#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.pfam as pfam

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'ozcan']
__requires__ = ['contigs-db', 'pfams-data',]
__provides__ = ['functions',]
__description__ = "Run Pfam on Contigs Database"


@time_program
def main():
    args = get_args()

    try:
        p = pfam.Pfam(args)
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
    parser.add_argument(*anvio.A('pfam-data-dir'), **anvio.K('pfam-data-dir'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('hmmer-program'), **anvio.K('hmmer-program'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
