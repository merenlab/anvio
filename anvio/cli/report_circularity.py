#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.bamops import CircularityPredictor
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2025, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["bam-file"]
__provides__ = ["contig-circularity-report-txt"]
__description__ = ("Predict contig circularity from paired-end read alignments in a given BAM file. "
                   "This program samples insert sizes, looks for RF pairs spanning "
                   "junctions, and reports per-contig circularity statistics.")


@terminal.time_program
def main():
    args = get_args()

    try:
        predictor = CircularityPredictor(args)
        predictor.process()
        predictor.predict_all(min_contig_length=args.min_contig_length,
                              output_file=args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('bam_file', metavar='BAM_FILE', help="Sorted and indexed BAM file to analyze.")
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))

    groupA = parser.add_argument_group('INSERT SIZE ESTIMATION', "Control how FR pairs are sampled to estimate insert sizes.")
    groupA.add_argument('--max-num-pairs-for-is-est', type=int, default=None,
                        help="Maximum number of FR pairs to sample when estimating insert size (default: 100000).")
    groupA.add_argument('--min-pairs-for-stats', type=int, default=None,
                        help="Minimum number of FR pairs required to compute insert size statistics (default: 1000).")

    groupB = parser.add_argument_group('CIRCULARITY PREDICTION', "Parameters used to decide whether contigs are circular. "
                        "They are indeed confusing, but if you read the documentation, you will probably be fine :)")
    groupB.add_argument('--insert-tolerance-factor', type=float, default=None,
                        help="MAD multiplier around the median insert size for acceptable circular junction inserts (default: 3.0).")
    groupB.add_argument('--min-supporting-pairs', type=int, default=None,
                        help="Absolute minimum number of RF pairs required to call a contig circular (default: 5).")
    groupB.add_argument('--expected-fraction-threshold', type=float, default=None,
                        help="Minimum fraction of expected RF pairs to call a contig circular (default: 0.3).")
    groupB.add_argument('--circularity-confidence-threshold', type=float, default=None,
                        help="Minimum circularity support value to call a contig circular (default: 0.5).")
    groupB.add_argument('--min-contig-length', type=int, default=None,
                        help="Only evaluate contigs at least this long (default: evaluate all contigs).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
