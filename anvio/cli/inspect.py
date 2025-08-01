#!/usr/bin/env python
# -*- coding: utf-8
"""This program runs inspect page directly"""

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.interactive as interactive

from anvio.argparse import ArgumentParser
from anvio.bottleroutes import BottleApplication
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ozcan']
__requires__ = ['profile-db', 'contigs-db', 'bin']
__provides__ = ['interactive', 'contig-inspection']
__description__ = "Start an anvi'o inspect interactive interface"
__resources__ = [("Visualizing contig coverages", "https://merenlab.org/2019/11/25/visualizing-coverages/")]


def main():
    args = get_args()
    run = terminal.Run()

    try:
        if not args.split_name and not args.just_do_it:
            raise ConfigError("You didn't provide a split name to `anvi-inspect`. If you don't care and want to start the "
                              "interactive interface with a random split from the profile database, please use "
                              "the flag `--just-do-it`")

        args.mode = 'inspect'
        d = interactive.Interactive(args)
        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)

    app = BottleApplication(d)
    app.run_application(args.ip_address, args.port_number)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DEFAULT INPUTS', "The interactive interface can be started with anvi'o databases.")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group('VISUALS RELATED', "Parameters that give access to various adjustements regarding\
                                                           the interface.")
    groupA.add_argument(*anvio.A('split-name'), **anvio.K('split-name'))
    groupB.add_argument(*anvio.A('hide-outlier-SNVs'), **anvio.K('hide-outlier-SNVs'))

    groupC = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupC.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupC.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupC.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))

    groupD = parser.add_argument_group('GENERAL CONVENIENCE', "From anvi'o developers to you.")
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
