#!/usr/bin/env python
# -*- coding: utf-8
"""Script to merge multiple profiles."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.merger as merger
import anvio.constants as constants

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['single-profile-db', 'contigs-db']
__provides__ = ['profile-db', 'misc-data-items-order']
__description__ = "Merge multiple anvio profiles"
__resources__ = [("Another description as part of the metagenomic workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile")]

def main():
    args = get_args()

    try:
        merger.MultipleRuns(args).merge()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('input', metavar = 'SINGLE_PROFILE(S)', nargs='+',
                        help = "Anvi'o single profiles to merge")

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    parser.add_argument(*anvio.A('sample-name'), **anvio.K('sample-name'))
    parser.add_argument(*anvio.A('description'), **anvio.K('description'))
    parser.add_argument(*anvio.A('skip-hierarchical-clustering'), **anvio.K('skip-hierarchical-clustering'))
    parser.add_argument(*anvio.A('enforce-hierarchical-clustering'), **anvio.K('enforce-hierarchical-clustering'))
    parser.add_argument(*anvio.A('distance'), **anvio.K('distance', {'default': None, 'help':
                      'The distance metric for the hierarchical clustering. If you do not use this flag,\
                       the default distance metric will be used for each clustering configuration\
                       which is "%s".' % constants.distance_metric_default}))
    parser.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'default': None, 'help':
                      'The same story with the `--distance`, except, the system default for this one\
                       is %s.' % constants.linkage_method_default}))
    parser.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
