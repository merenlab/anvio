#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.interacdome as interacdome

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__requires__ = ['contigs-db', 'interacdome-data']
__provides__ = ['binding-frequencies-txt', "misc-data-amino-acids"]
__description__ = "Run InteracDome on a contigs database"
__resources__ = [("Estimating per-residue binding frequencies with InteracDome", "http://merenlab.org/2020/07/22/interacdome/")]

@time_program
def main():
    args = get_args()
    try:
        p = interacdome.InteracDomeSuper(args)
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
    parser.add_argument(*anvio.A('interacdome-data-dir'), **anvio.K('interacdome-data-dir'))
    parser.add_argument(*anvio.A('interacdome-dataset'), **anvio.K('interacdome-dataset'))
    parser.add_argument(*anvio.A('min-binding-frequency'), **anvio.K('min-binding-frequency'))
    parser.add_argument(*anvio.A('min-hit-fraction'), **anvio.K('min-hit-fraction'))
    parser.add_argument(*anvio.A('information-content-cutoff'), **anvio.K('information-content-cutoff'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    parser.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'default': 'INTERACDOME'}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
