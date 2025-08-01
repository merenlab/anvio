#!/usr/bin/env python
# -*- coding: utf-8
"""A program to find palindromes in sequences."""

import sys

import anvio

from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.sequencefeatures import Palindromes
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__requires__ = ["dna-sequence", "fasta", "contigs-db"]
__provides__ = ['palindromes-txt']
__description__ = "A program to find palindromes in sequences"


@time_program
def main():
    args, unknown = get_args()

    try:
        if not (args.dna_sequence or args.contigs_db or args.fasta_file):
            raise ConfigError("You should send an input sequence source for this program to work.")

        if not args.output_file:
            args.verbose = True

        p = Palindromes(args)

        if args.contigs_db or args.fasta_file:
            p.process()
        elif args.dna_sequence:
            p.find(args.dna_sequence)
            p.report()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SEQUENCE SOURCE', "Where should anvi'o find your sequences?")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file'))
    groupA.add_argument(*anvio.A('dna-sequence'), **anvio.K('dna-sequence'))

    groupB = parser.add_argument_group('SEARCH ALGORITHM', "If your purpose is to search for palindromes in general, you "
                    "shouldn't change the default algorithm, which is based on BLAST. But there is an exceedingly faster "
                    "numba-based sequence search algorithm as an alternative for shorter sequences. There will be some "
                    "differences between BLAST and numba, and some of the parameters you will set for palindrome parameters "
                    "below will only be taken into consideration by BLAST.")
    groupB.add_argument(*anvio.A('palindrome-search-algorithm'), **anvio.K('palindrome-search-algorithm', {'default': 'BLAST'}))

    groupC = parser.add_argument_group('PALINDROME PROPERTIES', "Some essential stuff here about your palindromes. Please also "
                    "see the search sensitivity section below.")
    groupC.add_argument(*anvio.A('min-palindrome-length'), **anvio.K('min-palindrome-length'))
    groupC.add_argument(*anvio.A('max-num-mismatches'), **anvio.K('max-num-mismatches'))
    groupC.add_argument(*anvio.A('min-distance'), **anvio.K('min-distance'))
    groupC.add_argument(*anvio.A('min-mismatch-distance-to-first-base'), **anvio.K('min-mismatch-distance-to-first-base'))

    groupD = parser.add_argument_group('BLAST SENSITIVITY & PERFORMANCE', "These parameters are only relevant for BLAST.")
    groupD.add_argument('--blast-word-size', metavar="INT", type=int, default=10, help="This parameter is passed to "
                    "blastn as the `-word_size`, which literally means in the BLAST world the length of best perfect "
                    "among your alignments. The shorter the word size, the more short palindromes you will find, but it "
                    "will also influence your ability to find palindromes with mismatches. For instance, if the word size "
                    "is 10, then you will not find a palindrome that is 18 nt long but have a mismatch right at the 9th "
                    "nucleotide (because individual perfect matches in the alignment will have word sizes of less than 10). "
                    "So if you want to search for palindromes with a lot of mismatches, then you would like to keep your "
                    "word size small. But smaller word sizes will impact your performance negatively, and very small ones "
                    "will simply make it impossible to finish running for especially long contigs. If you want to find "
                    "palindromes with no mismatches, then you can safely match the word size to the minimum palindrome "
                    "length. If you want to do some test run, take a DNA sequence (say about 1,000 nts) that contains a "
                    "palindrome that looks like the kinds of palindromes you will be interested in finding, and run the "
                    "program `anvi-search-palindromes` with `--dna` parameter and `--verbose` flag. As you play with the "
                    "word size, minimum number of mismatches, and the minimum palindrome length, the output messages will "
                    "help you determine your best parameters.")
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    groupE = parser.add_argument_group('OUTPUT', "Output options.")
    groupE.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupE.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))

    return parser.parse_known_args()

if __name__ == "__main__":
    main()
