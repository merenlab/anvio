#!/usr/bin/env python
# -*- coding: utf-8

"""A client to access short reads in BAM files.

   See https://github.com/meren/anvio/issues/173 for details."""

import sys

import anvio
import anvio.terminal as terminal

from anvio.bamops import GetReadsFromBAM
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__requires__ = ['profile-db', 'contigs-db', 'bin', 'bam-file']
__provides__ = ['short-reads-fasta']
__description__ = ("Get short reads back from a BAM file with options for compression, splitting of "
                   "forward and reverse reads, etc")


@terminal.time_program
def main():
    args = get_args()

    try:
        r = GetReadsFromBAM(args)
        r.store_short_reads()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument('input_bams', metavar = 'BAM FILE[S]', nargs='+',
                        help = 'BAM file(s) to access to recover short reads')

    output_file_kwargs = {'help': 'File path(s) to store results. Multiple files should be separated by commas (no spaces).'}

    groupA = parser.add_argument_group('INPUT OPTION #1', "Work with good'ol anvi'o databases, collections, and bins, and get all the reads from "
                                            "one or more BAM files that match to the contigs matching to your selections.")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupA.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupA.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupB = parser.add_argument_group('INPUT OPTION #2', "Target a specific contig and/or specific start/stop positions to get all the matching "
                                            "reads from your BAM file(s).")
    groupB.add_argument(*anvio.A('target-contig'), **anvio.K('target-contig'))
    groupB.add_argument(*anvio.A('target-region-start'), **anvio.K('target-region-start'))
    groupB.add_argument(*anvio.A('target-region-end'), **anvio.K('target-region-end'))

    groupC = parser.add_argument_group('FILTERS', 'Choose which reads to work (or not to work) with like a pro.')
    groupC.add_argument(*anvio.A('fetch-filter'), **anvio.K('fetch-filter'))

    groupD = parser.add_argument_group('OUTPUT', 'Output options')
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupD.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    groupD.add_argument(*anvio.A('gzip-output'), **anvio.K('gzip-output'))
    groupD.add_argument(*anvio.A('split-R1-and-R2'), **anvio.K('split-R1-and-R2'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
