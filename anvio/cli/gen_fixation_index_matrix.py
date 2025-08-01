#!/usr/bin/env python
# -*- coding: utf-8
"""Fetches information from the variable positions table"""

import sys

import anvio

from anvio.argparse import ArgumentParser
from anvio.variabilityops import variability_engines
from anvio.errors import ConfigError, FilesNPathsError
from anvio.variabilityops import VariabilityFixationIndex


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__requires__ = ['contigs-db', 'profile-db', 'structure-db', 'bin', 'variability-profile-txt', 'splits-txt']
__provides__ = ['fixation-index-matrix']
__resources__ = [("Utilizing fixation index to study SAR11 population structure", "http://merenlab.org/data/sar11-saavs/#generating-distance-matrices-from-fixation-index-for-saavs-and-snvs-data"),
                 ("Measuring Distances Between Genomes in the Infant Gut Tutorial", "http://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst")]
__description__ = "Generate a pairwise matrix of a fixation indices between samples"


def main():
    args = get_args()

    try:
        if args.engine not in variability_engines:
            raise ConfigError("You are doing something wrong :/ Focus '%s' does not correspond to an available engine." % args.engine)

        FST = VariabilityFixationIndex(args)
        FST.process()
        FST.report()

    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DATABASES', "Declaring relevant anvi'o databases. First things first. Some\
                                                     are mandatory, some are optional.")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {"required": False}))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {"required": False}))
    groupA.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db', {"required": False}))
    groupA.add_argument(*anvio.A('variability-profile'), **anvio.K('variability-profile', {"required": False}))

    groupB = parser.add_argument_group('FOCUS :: BIN', "You need to pick someting to focus. You can ask anvi'o\
                                                             to work with a bin in a collection.")
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupC = parser.add_argument_group('FOCUS :: SPLIT NAMES', 'Alternatively you can declare split names to focus.')
    groupC.add_argument(*anvio.A('splits-of-interest'), **anvio.K('splits-of-interest'))

    groupD = parser.add_argument_group('FOCUS :: GENE CALLER IDs', 'Alternatively you can declare gene caller IDs to focus.')
    groupD.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupD.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    groupD.add_argument(*anvio.A('only-if-structure'), **anvio.K('only-if-structure'))

    groupE = parser.add_argument_group('SAMPLES', "You can ask anvi'o to focus only on a subset of samples.")
    groupE.add_argument(*anvio.A('samples-of-interest'), **anvio.K('samples-of-interest'))

    groupF = parser.add_argument_group('ENGINE', "Set your engine. This is important as it will define the output\
                                                  profile you will get from this program. The engine can focus on\
                                                  nucleotides (NT), codons (CDN), or an amino acids (AA).")
    groupF.add_argument(*anvio.A('engine'), **anvio.K('engine'))

    groupG = parser.add_argument_group('FILTERS', 'Parameters that will help you to do a very precise analysis.\
                                                  If you declare nothing from this bunch, you will get "everything"\
                                                  to play with, which is not necessarily a good thing...')
    groupG.add_argument(*anvio.A('min-coverage-in-each-sample'), **anvio.K('min-coverage-in-each-sample'))

    groupH = parser.add_argument_group('OUTPUT', 'Output file and style')
    groupH.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'default': 'fixation_indices.txt', 'metavar': 'FIXATION_INDICES'}))

    groupI = parser.add_argument_group('EXTRAS', 'Because why not be extra?')
    groupI.add_argument('--keep-negatives', action='store_true',
                        help = 'Negative numbers are theoretically possible, and are sometimes \
                                interpreted as out-breeding. By default, we set negative numbers to \
                                0 so the results are reflective of a standard distance metric. \
                                Provide this flag if you would prefer otherwise.')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
