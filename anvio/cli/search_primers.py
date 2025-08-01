#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.errors import ConfigError, FilesNPathsError
from anvio.sequencefeatures import PrimerSearch

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["short-reads-fasta"]
__requires__ = ["samples-txt", "primers-txt"]
__description__ = ("You provide this program with FASTQ files for one or more samples AND one or more primer sequences, "
                   "and it collects reads from FASTQ files that matches to your primers. This tool can be "
                   "most powerful if you want to collect all short reads from one or more metagenomes that are "
                   "downstream to a known sequence. Using the comprehensive output files you can analyze the "
                   "diversity of seuqences visually, manually, or using established strategies such as oligotyping.")


def main():
    args = get_args()

    try:
        s = PrimerSearch(args)
        s.process()
        s.print_summary()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES', "Here you are expected to declare your FASTQ files and sequences "
                            "which you are interested to find in those FASTQ files. Each file should have at least one "
                            "entry")
    groupA.add_argument(*anvio.A('samples-txt'), **anvio.K('samples-txt'))
    groupA.add_argument('--primers-txt', required=True, metavar='FILE', help="A two-column file that contains one "
                            "or more primer sequences. See the primers-txt artifact for details.")

    groupB = parser.add_argument_group('PARAMETERS OF LIFE OR DEATH', "Here you are expected to set appropriate "
                            "parameters for your search (or you can choose to go with the defaults)")
    groupB.add_argument('-m', '--min-remainder-length', metavar='INT', type=int, default=60,
                        help="Minimum length of the remainder of the read after a match. If your short read "
                              "is XXXMMMMMMYYYYYYYYYYYYYY, where Ms indicate the primer sequence, "
                              "min remainder length is equal to the length of nucleotide matching Y. Default is %(default)d.")
    groupB.add_argument('-x', '--only-report-primer-matches', default=False, action="store_true",
                        help="For a given sequence with a primer match, report only the part of it that matches to the primer. "
                             "In other words, if your short read is XXXMMMMMMYYYYYYYYYYYYYY, where Ms indicate part of your "
                             "sequence that matches to the primer, setting this flag will remove all Xs and Ys. Setting "
                             "this flag will automatically set the `--min-remainder-length` to 0. If you have a problem with "
                             "that, let us know.")
    groupB.add_argument('-R', '--only-report-remainders', default=False, action="store_true",
                        help="For a given sequence with a primer match, report only the part that follows the primer. "
                             "In other words, if your short read is XXXMMMMMMYYYYYYYYYYYYYY, where Ms indicate part of your "
                             "sequence that matches to the primer, setting this flag will only report Ys. This flag "
                             "is obviously incompatible with `--only-report-primer-matches` flag.")

    groupC = parser.add_argument_group('TESTING?', "Be our guest.")
    groupC.add_argument('--stop-after', metavar='INT', type=int, default=0, help="Stop after X number of hits because "
                               "who needs data.")

    groupD = parser.add_argument_group('OUTPUT', "Tell anvi'o where to put your thingies")
    groupD.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
