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
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False,
                                                                         'action': 'append',
                                                                         'help': "Path to an anvi'o profile database. Can be a merged profile-db "
                                                                                 "(multiple samples in one db) or a single-sample profile-db. You "
                                                                                 "can pass this flag multiple times to include several profile-dbs."}))
    groupA.add_argument(*anvio.A('profile-dbs-file'), **anvio.K('profile-dbs-file'))

    groupB = parser.add_argument_group(
        'EVENT FILTERS',
        "Which indel rows from the profile-db `indels` table to consider. Default focuses on MGE-sized events.")
    groupB.add_argument(*anvio.A('report-indels-min-length'), **anvio.K('report-indels-min-length'))
    groupB.add_argument(*anvio.A('report-indels-max-length'), **anvio.K('report-indels-max-length'))
    groupB.add_argument(*anvio.A('report-indels-event-type'), **anvio.K('report-indels-event-type'))

    groupC = parser.add_argument_group(
        'CLUSTERING',
        "Indel sequences are clustered with mmseqs2 easy-cluster to group homologous events together. "
        "By default, identity is length-tiered: stricter for short indels, looser for long ones, since "
        "longer alignments tolerate more divergence.")
    groupC.add_argument(*anvio.A('report-indels-min-identity'), **anvio.K('report-indels-min-identity'))
    groupC.add_argument(*anvio.A('report-indels-min-coverage'), **anvio.K('report-indels-min-coverage'))
    groupC.add_argument(*anvio.A('report-indels-position-tolerance'), **anvio.K('report-indels-position-tolerance'))
    groupC.add_argument(*anvio.A('report-indels-min-samples'), **anvio.K('report-indels-min-samples'))
    groupC.add_argument(*anvio.A('report-indels-min-reads-per-element'), **anvio.K('report-indels-min-reads-per-element'))

    groupD = parser.add_argument_group('OUTPUT', "Where the matrix + decoration files will be written.")
    groupD.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
