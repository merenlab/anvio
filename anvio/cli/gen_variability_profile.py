#!/usr/bin/env python
# -*- coding: utf-8
"""Fetches information from the variable positions table"""

import sys
from anvio.argparse import ArgumentParser

import anvio

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program
from anvio.variabilityops import variability_engines


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl', 'meren', 'ShaiberAlon']
__resources__ = [("All about SNVs, SCVs, and SAAVs", "http://merenlab.org/2015/07/20/analyzing-variability/"), ("This program in action in the anvi'o structure tutorial", "http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#supplying-anvi-display-structure-with-sequence-variability")]
__requires__ = ['contigs-db', 'profile-db', 'structure-db', 'bin', 'variability-profile', 'splits-txt']
__provides__ = ['variability-profile-txt']
__description__ = ("Generate a table that comprehensively summarizes the variability of nucleotide, "
                   "codon, or amino acid positions. We call these single nucleotide variants (SNVs), "
                   "single codon variants (SCVs), and single amino acid variants (SAAVs), respectively")

@time_program
def main():
    args = get_args()

    try:
        if args.engine not in variability_engines:
            raise ConfigError("You are doing something wrong :/ Focus '%s' does not correspond to an available engine." % args.engine)

        variability_engine = variability_engines[args.engine](args)
        variability_engine.process()
        variability_engine.report()
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
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupA.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db', {"required": False}))

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
    groupG.add_argument(*anvio.A('num-positions-from-each-split'), **anvio.K('num-positions-from-each-split'))
    groupG.add_argument(*anvio.A('min-departure-from-reference'), **anvio.K('min-departure-from-reference'))
    groupG.add_argument(*anvio.A('max-departure-from-reference'), **anvio.K('max-departure-from-reference'))
    groupG.add_argument(*anvio.A('min-departure-from-consensus'), **anvio.K('min-departure-from-consensus'))
    groupG.add_argument(*anvio.A('max-departure-from-consensus'), **anvio.K('max-departure-from-consensus'))
    groupG.add_argument(*anvio.A('min-occurrence-of-variable-positions'), **anvio.K('min-occurrence-of-variable-positions'))
    groupG.add_argument(*anvio.A('min-coverage-in-each-sample'), **anvio.K('min-coverage-in-each-sample'))
    groupG.add_argument(*anvio.A('quince-mode'), **anvio.K('quince-mode'))
    groupG.add_argument(*anvio.A('kiefl-mode'), **anvio.K('kiefl-mode'))

    groupH = parser.add_argument_group('OUTPUT', 'Output file and style')
    groupH.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'default': 'variability.txt', 'metavar': 'VARIABILITY_PROFILE'}))
    groupH.add_argument(*anvio.A('include-contig-names'), **anvio.K('include-contig-names'))
    groupH.add_argument(*anvio.A('include-split-names'), **anvio.K('include-split-names'))
    groupH.add_argument(*anvio.A('include-additional-data'), **anvio.K('include-additional-data'))
    groupH.add_argument(*anvio.A('include-site-pnps'), **anvio.K('include-site-pnps'))
    groupH.add_argument(*anvio.A('compute-gene-coverage-stats'), **anvio.K('compute-gene-coverage-stats'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
