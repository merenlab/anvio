#!/usr/bin/env python
# -*- coding: utf-8

import sys
import numpy as np
import pandas as pd

import anvio
import anvio.dbops as dbops
import anvio.bamops as bamops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = []
__requires__ = ["bam-file"]
__description__ = ("Report the distribution of template lengths from a BAM file. The purpose of "
                   "this is to get an idea about the insert size distribution in a BAM file "
                   "rapidly by summarizing distances between each paired-end read in a given "
                   "read recruitment experiment.")


@terminal.time_program
def main():
    args = get_args()

    try:
        bamops.PairedEndTemplateDist(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('bam_file', metavar = 'BAM_FILE', help = "An indexed BAM file")

    groupA = parser.add_argument_group('INPUT OPTIONS', "Things you don't care but should.")
    groupA.add_argument('--min-template-length-frequency', type=int, default=10, help="How many times a template "
                        "lenght should be observed to be considered as a viable template lenght? If "
                        "this number is zero, you will have an extremely large number template lengths "
                        "that will likely be due to noise. This number can be best set as a function of "
                        "the total number of mapped reads in a bam file. The default value is set with "
                        "the assumption that you have millions of reads in the BAM file.")
    groupA.add_argument('--max-template-length-to-consider', type=int, default=500000, help="Some paired end reads will "
                        "map incredibly distant locations given the linear sequences you have used for "
                        "read recruitment if the DNA templates were circular in the sample used for "
                        "sequencing. Which is a beautiful thing, but can ruin your histogram. Using this "
                        "parameter, you can discard paired end reads that are extremely distant from "
                        "one another.")


    groupB = parser.add_argument_group('OUTPUT OPTIONS', "Output options.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupB.add_argument('--plot-data', action='store_true', default=False, help="In addition to providing "
                        "you with a TAB-delimited output file, anvi'o can also try to plot a summary "
                        "histogram for the template length distribution across ALL contigs in a given BAM "
                        "file")

    args = parser.get_args(parser)


if __name__ == '__main__':
    main()
