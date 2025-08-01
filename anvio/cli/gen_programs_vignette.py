#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.programs as programs

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ['Xabier Vázquez-Campos']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Generate a markdown summary (vignette) of anvi'o programs"


def main():
    args = get_args()

    try:
        vignette = programs.ProgramsVignette(args)
        vignette.generate()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'default': 'vignette-out.md'}))
    parser.add_argument('-p', '--program-names-to-focus', default=None, help="Comma-spearated list of program names to focus\
                         Mostly for debugging purposes.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
