#!/usr/bin/env python
# -*- coding: utf-8
"""A program that computes rarefaction curves for a given pan-db"""

import sys

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.panops import RarefactionAnalysis
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2025, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ahenoch']
__resources__ = [("An example output of this program in the context of the Prochlorococcus metapangenome from Delmont and Eren 2018 is included in the pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/")]
__tags__ = ["pangenomics"]
__requires__ = ['pan-db']
__provides__ = ['rarefaction-curves']
__description__ = ("A program that computes rarefaction curves and Heaps' Law fit for a given pangenome")


@terminal.time_program
def main():
    args = get_args()

    try:
        analysis = RarefactionAnalysis(args)
        k, alpha = analysis.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('PAN DATABASE')
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))

    groupO = parser.add_argument_group('OUTPUT', "This program can generate a visualization of the rarefaction curves for you. If you provide nothing, it will simply print out"
                        "the Heaps' Law fit values for alpha and K for you to interpret. The format of the output file will be determined automatically by the extension. So if "
                        "you name your output file as `output.svg` you will get an SVG file, and if you name it as `output.pdf`, you will get a PDF.")
    groupO.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))

    groupE = parser.add_argument_group('PARAMETERS OF CONVENIENCE', "You want convenience and control? You have it (but not much).")
    groupE.add_argument('--iterations', default=100, type=int, help="How many subsampling iterations you wish to be performed? The default is %(default)d. Don't go below 10 :)")
    groupE.add_argument('--skip-output-files', default=False, action="store_true", help="Do not bother generating any output files, just report Heaps' fit values.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
