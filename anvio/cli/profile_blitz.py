#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.profiler import BAMProfilerQuick

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["bam-stats-txt"]
__requires__ = ["bam-file", "contigs-db"]
__description__ = ("FAST profiling of BAM files to get contig- or gene-level coverage and detection stats. "
                   "Unlike `anvi-profile`, which is another anvi'o program that can profile BAM files, this "
                   "program is designed to be very quick and only report long-format files for various "
                   "read recruitment statistics per item. Plase also see the program "
                   "`anvi-script-get-coverage-from-bam` for recovery of data from BAM files without an "
                   "anvi'o contigs database")


@terminal.time_program
def main():
    args = get_args()

    try:
        p = BAMProfilerQuick(args)
        p.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('bam_files', metavar = 'BAM_FILE(S)', nargs='+',
                        help = "One or more indexed BAM files")

    groupA = parser.add_argument_group('INPUT DB', "You will need to give this program an anvi'o contigs database.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))

    groupB = parser.add_argument_group('GENES?', "You can work with genes instead of contigs")
    groupB.add_argument(*anvio.A('gene-mode'), **anvio.K('gene-mode', {'help': ("This program "
                        "by default will summarize coverage and detection stats for contigs found "
                        "in your contigs database. Declaring this flag will change that behavior "
                        "and report coverage and detection stats for each gene. Brace yourself for "
                        "a huge file for large contigs databases lol :(")}))
    groupB.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller'))

    groupC = parser.add_argument_group('OUTPUT', "How do you want to store your output data.")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))
    groupC.add_argument(*anvio.A('report-minimal'), **anvio.K('report-minimal', {'help': "Using this flag, you can "
                        "ask anvi'o to report minimum amount of data about your genes or contigs (such as mean "
                        "coverage and detection) rather than a full blown output file with as much information as "
                        "anvi'o can offer (such as, mean coverage, detection, Q2Q3 coverage, standard deviation of "
                        "coverage, min/max values of coverage, GC-content and length of items, etc). Using this flag "
                        "can cut your processing time in half. See the help docs for example output files for contigs "
                        "and gene mode."}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

