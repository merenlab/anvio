#!/usr/bin/env python
# -*- coding: utf-8

"""The client for LinkMers class.

   See https://github.com/meren/anvio/issues/144 for details."""

import sys

import anvio

from anvio.errors import ConfigError, FilesNPathsError
from anvio.bamops import LinkMers


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['bam-file',]
__provides__ = ['linkmers-txt',]
__description__ = "Reports sequences stored in one or more BAM files that cover one of more specific nucleotide positions in a reference"


def main():
    args = get_args()

    try:
        r = LinkMers(args)
        r.process()
        r.report(args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument('-i', '--input-files', metavar = 'INPUT_BAM(S)', nargs='+', default = None, required = True,
                        help = 'Sorted and indexed BAM files to analyze. It is essential that all BAM files must be\
                                the result of mappings against the same contigs.')

    parser.add_argument(*anvio.A('contigs-and-positions'), **anvio.K('contigs-and-positions'))
    parser.add_argument(*anvio.A('only-complete-links'), **anvio.K('only-complete-links'))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))
    parser.add_argument(*anvio.A('list-contigs'), **anvio.K('list-contigs'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
