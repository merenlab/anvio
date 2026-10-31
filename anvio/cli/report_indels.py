#!/usr/bin/env python

import sys

import anvio
import anvio.indels as indels

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program, Run

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet']
__requires__ = ["contigs-db", "profile-db"]
__provides__ = []
__description__ = ("Aggregate large indels stored in one or more profile-db `indels` tables "
                   "across samples. Cluster the indels by sequence similarity, group events "
                   "by insertion site, and emit a loci x samples element-presence matrix "
                   "ready for `anvi-interactive --manual`. Useful for tracking mobile genetic "
                   "elements (ICEs, transposases, etc.) across metagenomes.")

run = Run()


@time_program
def main():
    args = get_args()

    try:
        I = indels.Indels(args)
        I.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group(
        'INPUT DATA',
        "One contigs-db and one or more profile-dbs. Use -p multiple times for repeated profile-dbs, "
        "OR pass a file listing them with --profile-dbs-file. A merged profile-db with multiple samples "
        "is also fully supported (just pass it once with -p).")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupA.add_argument('-p', '--profile-db', dest='profile_db', action='append', metavar='PROFILE_DB',
                        help="Path to an anvi'o profile database. Can be a merged profile-db (multiple samples "
                             "in one db) or a single-sample profile-db. You can pass this flag multiple times to "
                             "include several profile-dbs.")
    groupA.add_argument('--profile-dbs-file', metavar='FILE',
                        help="Alternative to repeating -p: a plain-text file with one profile-db path per line "
                             "(blank lines and lines starting with '#' are ignored). Mutually exclusive with -p.")

    groupB = parser.add_argument_group(
        'EVENT FILTERS',
        "Which indel rows from the profile-db `indels` table to consider. Default focuses on MGE-sized events.")
    groupB.add_argument('--min-length', type=int, default=500, metavar='INT',
                        help="Minimum indel length (in bp) to consider. Default: %(default)s. Drop this to also "
                             "report shorter structural variants; raise it to focus on really big elements.")
    groupB.add_argument('--max-length', type=int, default=None, metavar='INT',
                        help="Optional upper bound on indel length (in bp). Default: no upper bound.")
    groupB.add_argument('--event-type', choices=['INS', 'DEL', 'both'], default='both',
                        help="Which event types to include. INS = insertion (extra sequence in reads), "
                             "DEL = deletion (missing reference sequence in reads). Default: %(default)s. "
                             "Note that DEL sequences are reconstructed on the fly from the contigs-db, since "
                             "they are not stored in the indels table.")

    groupC = parser.add_argument_group(
        'CLUSTERING',
        "Indel sequences are clustered with mmseqs2 easy-cluster to group homologous events together. "
        "By default, identity is length-tiered: stricter for short indels, looser for long ones, since "
        "longer alignments tolerate more divergence.")
    groupC.add_argument('--min-identity', type=float, default=None, metavar='FLOAT',
                        help="Override the default length-tiered identity thresholds with a single global value "
                             "(e.g. 0.9). Must be in (0, 1]. Default: tiered (0.95 for 500 bp-1 kb, 0.90 for "
                             "1-5 kb, 0.85 for >5 kb).")
    groupC.add_argument('--min-coverage', type=float, default=0.8, metavar='FLOAT',
                        help="mmseqs2 -c / coverage threshold (bidirectional). Default: %(default)s.")
    groupC.add_argument('--position-tolerance', type=int, default=50, metavar='INT',
                        help="When collapsing events within a cluster into a single locus, two events on the same "
                             "contig are merged if their positions are within this many bp of the running mean. "
                             "Default: %(default)s.")
    groupC.add_argument('--min-samples', type=int, default=2, metavar='INT',
                        help="Drop loci where fewer than this many distinct samples carry the element with "
                             "presence > 0 in the final matrix. This counts samples on the presence side "
                             "(non-zero matrix cells) regardless of whether the signal came from an INS event, "
                             "a DEL event with partial allele frequency, or an auxiliary-coverage detection at "
                             "a DEL locus. INS-dominant and DEL-dominant loci are treated symmetrically. "
                             "Default: %(default)s. Set to 1 to keep loci where the element appears in only one "
                             "sample (useful for very small test datasets).")
    groupC.add_argument('--min-reads-per-element', type=int, default=5, metavar='INT',
                        help="Minimum number of reads (summed across all events and samples) that must "
                             "support a clustered element (= MGE family) for it to be reported. Counts both "
                             "INS-supporting and DEL-supporting reads from any sample. Use this to suppress "
                             "low-occurrence elements that the data doesn't really see. Default: %(default)s. "
                             "Set to 1 to disable this filter.")

    groupD = parser.add_argument_group('OUTPUT', "Where the matrix + decoration files will be written.")
    groupD.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
