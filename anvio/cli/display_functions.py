#!/usr/bin/env python
# -*- coding: utf-8
"""A program to display the distribution of functions across genomes."""

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
__authors__ = ['meren']
__requires__ = ['functions', 'genomes-storage-db', 'internal-genomes', 'external-genomes', 'groups-txt']
__provides__ = ['interactive', 'functional-enrichment-txt']
__description__ = "Start an anvi'o interactive display to see functions across genomes"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_all_npm_packages_installed()
        args.mode = 'functional'
        d = interactive.Interactive(args)
        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except DictIOError as e:
        print(e)
        sys.exit(-3)

    if args.dry_run:
        run.info_single('Dry run? Kthxbai.', nl_after=1, nl_before=1)
        sys.exit()

    app = BottleApplication(d)
    app.run_application(args.ip_address, args.port_number)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('GENOMES', "Tell anvi'o where your genomes are.")
    groupA.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))

    groupB = parser.add_argument_group('GROUPS', "If you want, you can also tell anvi'o how to group your genomes "
                                "so it can also compute functional enrichment between them.")
    groupB.add_argument(*anvio.A('groups-txt'), **anvio.K('groups-txt'))
    groupB.add_argument('--print-genome-names-and-quit', default=False, action='store_true', help="Sometimes, "
                                "especially when you are interested in creating a groups file for your genomes you "
                                "gather from multiple different sources, it may be difficult to know every single "
                                "genome name that will go into your analysis. If you declare this flag, after "
                                "initializing everyghing, anvi'o will print out every genome name it found and quit, "
                                "so you can actually put together a groups file for them.")

    groupC = parser.add_argument_group('FUNCTIONS', "Tell anvi'o which functional annotation source you like above all, and other "
                                "important details you like about your analysis.")
    groupC.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source', {'required': True}))
    groupC.add_argument(*anvio.A('aggregate-based-on-accession'), **anvio.K('aggregate-based-on-accession'))
    groupC.add_argument(*anvio.A('aggregate-using-all-hits'), **anvio.K('aggregate-using-all-hits'))
    groupC.add_argument('--min-occurrence', metavar="NUM GENOMES", default=1, help=("The minimum number of occurrence of any "
                                "given function accross genomes. If you set a value, those functions that occur in less number "
                                "of genomes will be excluded."), type=int)

    groupD = parser.add_argument_group('PROFILE DB', "To store visuals state, collections, and such. It will be AUTOMATICALLY "
                                "generated for you, and you can't use an existing profile for this. But then once it is generated "
                                "you can use that profile with `anvi-interactive`. It is actually objectively very cool.")
    groupD.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))

    groupE = parser.add_argument_group('GENES', "By default, anvi'o will look for genes in contigs databases that are identified "
                                "by `pyrodigal-gv`. But if you have generated your contigs databse with external gene calls, or have "
                                "otherwise used another gene caller than the default, you can explicitly ask anvi'o to use that "
                                "one to recover your genes.")
    groupE.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller'))

    groupF = parser.add_argument_group('VISUALS RELATED', "Parameters that give access to various adjustements regarding\
                                                           the interface.")
    groupF.add_argument(*anvio.A('title'), **anvio.K('title'))
    groupF.add_argument(*anvio.A('state-autoload'), **anvio.K('state-autoload'))
    groupF.add_argument(*anvio.A('collection-autoload'), **anvio.K('collection-autoload'))
    groupF.add_argument(*anvio.A('export-svg'), **anvio.K('export-svg'))

    groupG = parser.add_argument_group('SWEET PARAMS OF CONVENIENCE', "Parameters and flags that are not quite essential (but\
                                                                       nice to have).")
    groupG.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))
    groupG.add_argument(*anvio.A('skip-news'), **anvio.K('skip-news'))

    groupH = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupH.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupH.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupH.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupH.add_argument(*anvio.A('read-only'), **anvio.K('read-only'))
    groupH.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupH.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))
    groupH.add_argument(*anvio.A('user-server-shutdown'), **anvio.K('user-server-shutdown'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
