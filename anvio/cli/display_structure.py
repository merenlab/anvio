#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.interactive as interactive
from anvio.bottleroutes import BottleApplication

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl', 'ozcan']
__requires__ = ['structure-db', 'variability-profile-txt', 'contigs-db', 'profile-db', 'splits-txt']
__provides__ = ['interactive']
__resources__ = [("The overview page from the release", "http://merenlab.org/software/anvio-structure/"),
                 ("The section of the Infant Gut Tutorial focused on anvi-display-structure", "http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures"),
                 ("Integrating sequence variants and predicted protein structures", "http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/")]
__description__ = "Interactively visualize sequence variants on protein structures"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_all_npm_packages_installed()

        args.mode = 'structure'

        d = interactive.StructureInteractive(args)
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

    groupS = parser.add_argument_group("STRUCTURE", "Information related to the structure database, which can be \
                                                     created with anvi-gen-structure-database.")
    groupV = parser.add_argument_group("VARIABILITY", "We can overlay codon and amino acid variability in your \
                                                       metagenomes but we need a data source of this variability. \
                                                       Most simply, anvi'o can learn this information when you \
                                                       provide both your profile (-p) and contigs (-c) databases. \
                                                       Alternatively, you can provide a variability table output \
                                                       (-V) from the program anvi-gen-variability-profile. If you \
                                                       don't want to visualize variants, this is the wrong tool for \
                                                       the job. Instead, export the PDB files with \
                                                       anvi-export-structures, and open with a \
                                                       more comprehensive protein viewing software.")
    groupR = parser.add_argument_group("REFINING PARAMETERS", "Which samples, genes, and contigs etc. are you \
                                                               interested in? Define that stuff here.")
    groupP = parser.add_argument_group("SERVER CONFIGURATION", "For power users.")

    min_dfc_dict = {
        "help": "Takes a value between 0 and 1, where 1 is maximum divergence from the consensus. "
                "it can be an expensive operation to display every variable position, and so the "
                "default is 0.05. To display every variable position, set this parameter to 0.",
        "default": 0.05
    }

    groupS.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db'))
    groupV.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {"required": False}))
    groupV.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {"required": False}))
    groupV.add_argument(*anvio.A('variability-profile'), **anvio.K('variability-profile'))
    groupR.add_argument(*anvio.A('splits-of-interest'), **anvio.K('splits-of-interest'))
    groupR.add_argument(*anvio.A('samples-of-interest'), **anvio.K('samples-of-interest'))
    groupR.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupR.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    groupR.add_argument(*anvio.A('min-departure-from-consensus'), **anvio.K('min-departure-from-consensus', min_dfc_dict))
    groupR.add_argument('--SAAVs-only', required=False, action='store_true', help = 'If provided, variability will \
                                                                     be generated for single amino \
                                                                     acid variants (SAAVs) and not \
                                                                     for single codon variants \
                                                                     (SCVs).  This could save you \
                                                                     some time if you\'re only \
                                                                     interested in SAAVs.')
    groupR.add_argument('--SCVs-only', required=False, action='store_true', help = 'If provided, variability will \
                                                                     be generated for single codon \
                                                                     variants (SCVs) and not \
                                                                     for single amino acid variants \
                                                                     (SAAVs).  This could save you \
                                                                     some time if you\'re only \
                                                                     interested in SCVs.')
    groupP.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupP.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupP.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupP.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupP.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
