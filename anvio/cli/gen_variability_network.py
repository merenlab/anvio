#!/usr/bin/env python
# -*- coding: utf-8
"""Creates a network file for a given variability profile output"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.variabilityops import VariabilityNetwork


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["variability-profile-xml"]
__requires__ = ["variability-profile-txt"]
__resources__ = [("An applicatio of this program in the  Infant Gut Tutorial", "https://merenlab.org/tutorials/infant-gut/#visualizing-snv-profiles-as-a-network")]
__description__ = ("Generate a network description from an anvi'o variability profile.")


def main():
    args, VN = get_args()

    if args.list_competing_NT_calculators:
        run = terminal.Run()

        run.warning("Following are the available options for competing NT calculations:",
                    header="COMPETING NT CALCULATORS", lc="green")
        for calculator_option in VN.competing_NT_calculators:
            run.info(f"'{calculator_option}'", VN.competing_NT_calculators[calculator_option]['description'], nl_after=1)

        sys.exit()

    try:
        variable_nt_positions = VariabilityNetwork(args)
        variable_nt_positions.generate()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    import argparse
    VN = VariabilityNetwork(argparse.Namespace())
    competing_NT_calculator_options = [f"'{k}'" for k in VN.competing_NT_calculators]

    groupA = parser.add_argument_group('INPUT', 'What goes in..')
    groupA.add_argument('-i', '--input-file', metavar='VARIABILITY_PROFILE', required=True,
                        help="The anvi'o variability profile. Please see `anvi-gen-variability-profile` to "
                             "generate one.")

    groupB = parser.add_argument_group('OUTPUT', 'What comes out..')
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))
    groupB.add_argument('--as-matrix', default=False, action="store_true", help ="Report the output as a table, rather than "
                            "an network")


    groupC = parser.add_argument_group('ADDITIONAL PARAMS', 'More options for the wise.')
    groupC.add_argument(*anvio.A('max-num-unique-positions'), **anvio.K('max-num-unique-positions'))
    groupC.add_argument('--include-competing-NTs', metavar='OPTION', default=False, required=False,
                        help =f"Whether anvi'o should use competing nucleotides, rather than unique positions for "
                              f"variable nucleotides (which is the default behavior), to define the nodes in the "
                              f"network. If you would like that, then you will need to explicitly ask for a calculator "
                              f"for competing nuclotides. The available options are the following: {', '.join(competing_NT_calculator_options)}. "
                              f"If you would like to learn about these options before choosing one of them (as "
                              f"you should), please first use the flag `--list-competing-NT-calculators` (in which "
                              f"case the program will list descriptions for these options and exit.")
    groupC.add_argument('--edge-variable', metavar='OPTION', default="departure_from_reference", required=False,
                        help ="In the output anvi'o will connect SNVs to samples using edges, weights of which will be "
                            "defined by this variable.")
    groupC.add_argument('--list-competing-NT-calculators', default=False, action="store_true",
                        help ="List options and their descriptions that can be used with `--include-competing-NTs`.")

    return parser.get_args(parser), VN


if __name__ == '__main__':
    main()
