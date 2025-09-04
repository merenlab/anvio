#!/usr/bin/env python
# -*- coding: utf-8
"""Get codon or amino acid frequency statistics from genomes, genes, and functions."""

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
__authors__ = ['meren', 'semiller10']
__requires__ = ['contigs-db', 'profile-db', 'collection', 'bin']
__provides__ = ['interactive', 'svg']
__description__ = ("Display codon frequency statistics across genes in a given genome in the anvi'o interactive interface.")


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()


    try:
        utils.is_all_npm_packages_installed()
        args.mode = 'codon-frequencies'
        d = interactive.Interactive(args)

        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)

    if args.dry_run:
        run.info_single('Dry run? Kthxbai.', nl_after=1, nl_before=1)
        sys.exit()

    app = BottleApplication(d)
    app.run_application(args.ip_address, args.port_number)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('GENOME', "Show us where to find the genome to Display codons "
        "across genes in a single genome. It could also be a bin stored in a collection.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupA.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupA.add_argument(
        '--gene-caller-ids', type=int, nargs='+', help="Select genes by ID, space-separated.")

    groupC = parser.add_argument_group('OUTPUT PARAMETERS')
    groupC.add_argument('--infinity-to-zero', default=False, action='store_true', help="Replace "
             "NA (empty) values in output with 0.0. NA occurs with `--synonymous` when all "
             "codons for an amino acid are absent in a gene or function, resulting in 0/0, "
             "reported as NA. Use with caution.")

    groupD = parser.add_argument_group('FREQUENCY STATISTICS', "How should codon frequencies be computed.")
    groupD.add_argument(*anvio.A('relative'), **anvio.K('relative'))
    groupD.add_argument(*anvio.A('synonymous'), **anvio.K('synonymous'))
    groupD.add_argument(*anvio.A('return-amino-acids'), **anvio.K('return-amino-acids'))
    groupD.add_argument(*anvio.A('encodings-txt'), **anvio.K('encodings-txt'))
    groupD.add_argument(*anvio.A('sum'), **anvio.K('sum'))
    groupD.add_argument(*anvio.A('average'), **anvio.K('average'))

    groupF = parser.add_argument_group('FILTER GENES, FUNCTIONS, CODONS',
        "Genes/functions can be filtered by the number of codons they contain, e.g., ignore genes "
        "shorter than 300 codons. Codons can be selected a priori, e.g., ignore Ala codons, or "
        "rarer codons can be excluded, e.g., ignore amino acids that are decoded by ≤3 codons in "
        "≥90%% of genes. Filters can improve the statistical utility of codon relative frequency "
        "data.")
    groupF.add_argument(*anvio.A('gene-min-codons'), **anvio.K('gene-min-codons'))
    groupF.add_argument(*anvio.A('function-min-codons'), **anvio.K('function-min-codons'))
    groupF.add_argument(*anvio.A('exclude-amino-acids'), **anvio.K('exclude-amino-acids'))
    groupF.add_argument(*anvio.A('include-amino-acids'), **anvio.K('include-amino-acids'))
    groupF.add_argument(*anvio.A('sequence-min-amino-acids'), **anvio.K('sequence-min-amino-acids'))
    groupF.add_argument(*anvio.A('pansequence-min-amino-acids'), **anvio.K('pansequence-min-amino-acids'))
    groupF.add_argument(*anvio.A('min-codon-filter'), **anvio.K('min-codon-filter'))

    groupG = parser.add_argument_group('VISUALS RELATED', "Parameters that give access to various adjustements regarding\
                                                           the interface.")
    groupG.add_argument(*anvio.A('title'), **anvio.K('title'))
    groupG.add_argument(*anvio.A('state-autoload'), **anvio.K('state-autoload'))
    groupG.add_argument(*anvio.A('collection-autoload'), **anvio.K('collection-autoload'))
    groupG.add_argument(*anvio.A('export-svg'), **anvio.K('export-svg'))

    groupH = parser.add_argument_group('SWEET PARAMS OF CONVENIENCE', "Parameters and flags that are not quite essential (but\
                                                                       nice to have).")
    groupH.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))
    groupH.add_argument(*anvio.A('skip-news'), **anvio.K('skip-news'))


    groupX = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupX.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupX.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupX.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupX.add_argument(*anvio.A('read-only'), **anvio.K('read-only'))
    groupX.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupX.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))
    groupX.add_argument(*anvio.A('user-server-shutdown'), **anvio.K('user-server-shutdown'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
