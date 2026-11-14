#!/usr/bin/env python

import sys

import anvio

from anvio.amrfinder import AMRFinder
from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']
__requires__ = ['contigs-db', 'amrfinder-data']
__provides__ = ['functions']
__description__ = ("You can use this to annotate genes in a contigs-db with NCBI's AMRFinderPlus to find antimicrobial "
                   "resistance genes, stress response genes (such as metal and biocide resistance), and virulence genes.")


@time_program
def main():
    args = get_args()

    try:
        amrfinder = AMRFinder(args)
        amrfinder.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('amrfinder-data-dir'), **anvio.K('amrfinder-data-dir'))
    parser.add_argument(*anvio.A('full-amrfinder-report'), **anvio.K('full-amrfinder-report'))
    parser.add_argument('--skip-amrfinder-e-values', default=False, action='store_true', help="By default, AMRFinderPlus "
                        "does not report e-values for its HMM hits because each model is manually curated and has its own "
                        "threshold. However, because it's nice to know, anvi'o recovers them from the hmmsearch table "
                        "AMRFinderPlus leaves behind in its temporary directory. We can't see into the future and it is "
                        "possible that future versions of AMRFinderPlus won't produce that temporary file anymore. In that "
                        "case anvi'o would be confused and tell you so. You can use this flag to tell anvi'o to go on with "
                        "annotating genes with AMRFinderPlus without snatching the e-values, in which case they will all be 0.")
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it', {'help': "Using this flag, anvi'o will overwrite "
                        "previous AMRFinderPlus annotations in your contigs db."}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
