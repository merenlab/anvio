#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.kegg as kegg

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'semiller10']
__requires__ = ["contigs-db", "kegg-data",]
__provides__ = ["kegg-functions", "functions"]
__description__ = "Run KOfam HMMs on an anvi'o contigs database"

@time_program
def main():
    args = get_args()

    try:
        p = kegg.RunKOfams(args)
        p.process_kofam_hmms()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupR = parser.add_argument_group('REQUIRED INPUT', 'The stuff you need for this to work.')
    groupR.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupO = parser.add_argument_group('OPTIONAL INPUT', "Optional params for a custom experience.")
    groupO.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir'))
    groupO.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupO.add_argument(*anvio.A('hmmer-program'), **anvio.K('hmmer-program'))
    groupO.add_argument(*anvio.A('include-stray-KOs'), **anvio.K('include-stray-KOs'))
    groupO.add_argument(*anvio.A('keep-all-hits'), **anvio.K('keep-all-hits'))
    groupO.add_argument(*anvio.A('log-bitscores'), **anvio.K('log-bitscores'))
    groupO.add_argument(*anvio.A('skip-brite-hierarchies'), **anvio.K('skip-brite-hierarchies'))
    groupO.add_argument(*anvio.A('no-hmmer-prefiltering'), **anvio.K('no-hmmer-prefiltering'))
    groupO.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    groupH = parser.add_argument_group('BITSCORE RELAXATION HEURISTIC', "Sometimes, KEGG-provided bitscore thresholds "
                                                "are too stringent, causing us to miss valid annotations. We apply the "
                                                "following heuristic to relax those thresholds and annotate genes that "
                                                "would otherwise miss a valid annotations: for every gene call that is "
                                                "not annotated, examine hits to that gene with an evalue <= X and a bitscore "
                                                " > Y * KEGG's threshold (where Y is a float from 0 to 1). If those hits "
                                                "are all to a unique KO profile, annotate the gene with that KO. Long story "
                                                "over, you can set X and Y using the parameters below.")
    groupH.add_argument(*anvio.A('skip-bitscore-heuristic'), **anvio.K('skip-bitscore-heuristic'))
    groupH.add_argument(*anvio.A('heuristic-e-value'), **anvio.K('heuristic-e-value'))
    groupH.add_argument(*anvio.A('heuristic-bitscore-fraction'), **anvio.K('heuristic-bitscore-fraction'))


    return parser.get_args(parser)


if __name__ == '__main__':
    main()
