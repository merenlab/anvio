#!/usr/bin/env python
# -*- coding: utf-8

import sys
import anvio
from anvio.argparse import ArgumentParser

import anvio.structureops as structops

from anvio.errors import ConfigError, FilesNPathsError, ModellerError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__requires__ = ["contigs-db", "structure-db"]
__description__ = ("Add or re-run genes from an already existing structure database. All settings used "
                   "to generate your database will be used in this program")


def main():
    args = get_args()

    try:
        structops.StructureSuperclass(args, create=False)._run()
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)
    except ModellerError as e:
        print(e)
        sys.exit(3)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupD = parser.add_argument_group('DATABASES', 'Declaring relevant anvi\'o databases. First things first.')
    groupG = parser.add_argument_group('GENES', 'Specify which genes you want to be modelled. If a gene already '
                                                'exists in the DB, it will be overwritten if --overwrite is set. '
                                                'Otherwise, an error will be raised.')
    groupO = parser.add_argument_group('OUTPUT', 'Output file and output style.')
    groupM = parser.add_argument_group('MODELLER PARAMS', 'Parameters for MODELLER\'s homology modeling.')
    groupE = parser.add_argument_group('EXTRA', 'Everything else.')

    groupD.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupD.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db'))
    groupG.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupG.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    groupG.add_argument(*anvio.A('external-structures'), **anvio.K('external-structures'))
    groupO.add_argument(*anvio.A('dump-dir'), **anvio.K('dump-dir'))
    groupM.add_argument('--list-modeller-params', action='store_true', help='Since you are updating an existing DB, modeller params are set in '
                                                                            'place. You can have this program list them by providing this flag')

    groupE.add_argument("--rerun-genes", action='store_true', help = \
                        """Supply if you would like to rerun structural modelling for your genes of
                        interest if they are already present in your DB""")
    groupE.add_argument("--modeller-executable", type=str, help = \
                        """The MODELLER program to use. For example, `mod9.19`. Anvi'o will try and find
                        it if not provided.""")

    groupE.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    params = {
        'default': 25,
        'help': anvio.K('write-buffer-size-per-thread')['help'] + \
                ' If --num-threads is 1, this parameter is ignored because the DB is written to after each gene'
    }
    groupE.add_argument(*anvio.A('write-buffer-size-per-thread'), **anvio.K('write-buffer-size-per-thread', params))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
