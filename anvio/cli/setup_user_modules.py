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
__authors__ = ['ivagljiva']
__requires__ = ["user-modules-data"]
__provides__ = ["modules-db", "user-modules-data"]
__description__ = "Set up user-defined metabolic pathways into an anvi'o-compatible database"


@time_program
def main():
    args = get_args()

    args.kegg_archive = None
    args.kegg_snapshot = None
    args.download_from_kegg = None

    try:
        setup = kegg.KeggSetup(args)
        setup.setup_user_data()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    groupU = parser.add_argument_group('USER DATA SETUP', "If you have defined your own set "
                                       "of metabolic modules, you can use this program to parse these "
                                       "files into the format required for metabolism estimation.")
    groupU.add_argument(*anvio.A('user-modules'), **anvio.K('user-modules'))

    groupE = parser.add_argument_group('EXTRAS', "Extras for the extra.")
    groupE.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir', {'help': "You may need to provide the location of "
                                                                               "your existing KEGG data so that we can properly "
                                                                               "sanity check your data. Use this parameter to do so."}))
    groupE.add_argument(*anvio.A('reset'), **anvio.K('reset'))
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
