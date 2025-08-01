#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.interactive as interactive
from anvio.bottleroutes import BottleApplication

from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva']
__requires__ = ['contigs-db', "kegg-data", "kegg-functions", "profile-db", "collection", "bin"]
__provides__ = ['interactive']
__description__ = "Start the anvi'o interactive interactive for viewing KEGG metabolism data"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_all_npm_packages_installed()

        args.mode = 'metabolism'
        d = interactive.MetabolismInteractive(args)
        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)

        app = BottleApplication(d)
        app.run_application(args.ip_address, args.port_number)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except DictIOError as e:
        print(e)
        sys.exit(-3)


def get_args():
    parser = ArgumentParser(description=__description__)
    groupI = parser.add_argument_group('INPUT', "The minimum you must provide this program is a contigs database. In which case "
                                                "anvi'o will attempt to estimate and display metabolism for all contigs in it, assuming that "
                                                "the contigs database represents a single genome. If the contigs database is actually "
                                                "a metagenome, you should use the `--metagenome` flag to explicitly declare that.")
    groupI.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    groupI.add_argument(*anvio.A('metagenome-mode'), **anvio.K('metagenome-mode'))
    groupI.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir'))

    groupP = parser.add_argument_group('ADDITIONAL INPUT', "If you also provide a profile database AND a collection name, anvi'o will "
                                                           "estimate metabolism separately for each bin in your collection. You can also limit "
                                                           "those estimates to a specific bin or set of bins in the collection using the parameters "
                                                           "`--bin-id` or `--bin-ids-file`, respectively.")
    groupP.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupP.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupP.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupP.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupC = parser.add_argument_group('OUTPUT', "Parameters for controlling estimation output. The output will be TAB-delimited files which by "
                                                 "default are prefixed with 'kegg-metabolism', but you can of course change that name here.")
    groupC.add_argument(*anvio.A('module-completion-threshold'), **anvio.K('module-completion-threshold'))

    groupB = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupB.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupB.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupB.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupB.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupB.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
