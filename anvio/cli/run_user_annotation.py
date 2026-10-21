#!/usr/bin/env python

import sys

import anvio

from anvio.user_annotation import UserAnnotationRunner
from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program
from anvio import constants


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['l-gallucci', 'iwilkie']
__requires__ = ['contigs-db', 'user-annotation-db']
__provides__ = ['hmm-hits', 'functions']
__description__ = ("Run custom functional annotation using databases prepared by "
                   "`anvi-setup-user-annotation-db`. HMM sources are searched with "
                   "hmmscan/hmmsearch and DIAMOND sources with diamond blastp. All "
                   "results are stored in the gene_functions table.")


@time_program
def main():
    args = get_args()

    try:
        runner = UserAnnotationRunner(args)
        runner.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("REQUIRED", "The contigs database and the prepared annotation databases.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupA.add_argument('--annotation-db-dir', metavar='DIR',
                        default=constants.default_user_annotation_data_dir,
                        help="Path to the directory created by `anvi-setup-user-annotation-db`. "
                             "It must contain a manifest.json file that lists the prepared databases. "
                             "Default: %(default)s")

    groupB = parser.add_argument_group("DATABASE SELECTION",
                                       "By default all databases in the annotation directory are used. "
                                       "Use `anvi-setup-user-annotation-db --list` to see what is registered.")
    groupB.add_argument('--database', default=None, metavar='NAME[,NAME,...]',
                        help="A comma-separated list of database names (as defined in the input TSV "
                             "given to `anvi-setup-user-annotation-db`) to run. If omitted, all "
                             "databases in the annotation directory are used. Use 'all' to run "
                             "every database explicitly.")

    groupC = parser.add_argument_group("PERFORMANCE")
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupC.add_argument(*anvio.A('hmmer-program'), **anvio.K('hmmer-program'))

    groupChmm = parser.add_argument_group("HMM OPTIONS", "Applies only to HMM-based databases.")
    groupChmm.add_argument('--cut-tc', default=None, metavar='FILE',
                           help="Path to a file with per-model custom trusted cutoff (TC) values. "
                                "Two formats accepted: tab-delimited (model_name<TAB>seq_tc[<TAB>dom_tc]) or "
                                "colon format (model_name: seq_tc), e.g. '[FeFe]: 15.9'. "
                                "Domain score (dom_tc) defaults to seq_tc when omitted. "
                                "Models listed here override their embedded TC (or fallback evalue) and are "
                                "searched with --cut_tc. All other models continue to use their embedded "
                                "TC/GA/NC annotations or the evalue fallback. Models not found in any loaded "
                                "database are silently ignored.")

    groupD = parser.add_argument_group("DIAMOND OPTIONS", "Applies only to DIAMOND-based databases.")
    groupD.add_argument('--evalue', default=anvio.K('min-e-value')['default'], type=float, metavar='FLOAT',
                        help="E-value cutoff for DIAMOND blastp searches. "
                             "Default inherits the anvio standard: %(default)g. "
                             "For HMM-based databases the threshold is set at setup time from the "
                             "profile annotations (TC/GA/NC) and cannot be changed here.")
    groupD.add_argument('--min-pident', default=None, type=float, metavar='FLOAT',
                        help="Minimum percent identity for DIAMOND hits (0–100). Hits below this "
                             "threshold are discarded before being stored in the gene_functions "
                             "table. Applies only to DIAMOND searches.")
    groupD.add_argument('--qcov', default=None, type=float, metavar='FLOAT',
                        help="Minimum query coverage for DIAMOND blastp hits (0–100 percent). "
                             "Hits where the aligned region covers less than this fraction of the "
                             "query sequence are discarded. Applies only to DIAMOND searches.")
    groupD.add_argument('--max-hsps', default=None, type=int, metavar='INT',
                        help="Maximum number of HSPs per target sequence per query for DIAMOND blastp. "
                             "Applies only to DIAMOND searches. Default: DIAMOND's built-in default.")
    groupD.add_argument('--diamond-sensitivity', default=None, metavar='MODE',
                        choices=['fast', 'mid-sensitive', 'sensitive', 'very-sensitive', 'ultra-sensitive'],
                        help="DIAMOND sensitivity mode. Options (slowest = most sensitive): "
                             "fast, mid-sensitive, sensitive, very-sensitive, ultra-sensitive. "
                             "Default: DIAMOND's built-in default (between fast and mid-sensitive). "
                             "Applies only to DIAMOND searches.")

    groupE = parser.add_argument_group("AUTHORITY")
    groupE.add_argument('--force-overwrite', default=False, action='store_true',
                        help="By default, databases that have already been annotated in the contigs "
                             "database are skipped (a warning is printed). Use this flag to drop "
                             "existing annotations and re-annotate. Each database is evaluated "
                             "independently, so you can safely mix already-annotated and new databases "
                             "in one run — only the ones needing re-annotation require this flag.")
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
