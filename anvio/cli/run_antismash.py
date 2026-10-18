#!/usr/bin/env python

import sys

import anvio
import anvio.antismash as antismash

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError, CommandError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']
__requires__ = ['contigs-db']
__can_use__ = ['external-genomes']
__provides__ = ['functions']
__description__ = ("Run antiSMASH on one or many anvi'o contigs databases and import biosynthetic gene "
                   "cluster (BGC) annotations into the gene functions table")


@time_program
def main():
    args = get_args()

    try:
        antismash.Antismash(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)
    except CommandError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("INPUT", "The contigs database(s) to annotate. Provide either a single "
                             "one with -c, or many at once through an external-genomes file with -e.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': False}))

    groupB = parser.add_argument_group("ANTISMASH", "Options forwarded to antiSMASH.")
    groupB.add_argument(*anvio.A('antismash-dir'), **anvio.K('antismash-dir'))
    groupB.add_argument(*anvio.A('antismash-data-dir'), **anvio.K('antismash-data-dir'))
    groupB.add_argument(*anvio.A('antismash-taxon'), **anvio.K('antismash-taxon'))
    groupB.add_argument(*anvio.A('min-contig-length'), **anvio.K('min-contig-length', {'default': 3000,
                        'help': "antiSMASH cannot meaningfully call clusters on short fragments, so anvi'o "
                                "filters out contigs shorter than this before the BGC search. The default is "
                                "%(default)d bp, which is also antiSMASH's own default."}))
    groupB.add_argument(*anvio.A('antismash-add-on'), **anvio.K('antismash-add-on'))

    groupC = parser.add_argument_group("OUTPUT", "anvi-run-antismash by default annotates the contigs-db with "
                             "smCOGs and identifies BGC clusters (antiSMASH genes and antiSMASH regions) and also "
                             "saves a browsable HTML report plus a per-region summary. You can use this flag:")
    groupC.add_argument(*anvio.A('include-detailed-antismash-output'), **anvio.K('include-detailed-antismash-output'))
    groupC.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'help': "Used to define where to save "
                        "the antiSMASH report and summary files. It is an optional flag and if undefined, "
                        "anvi-run-antismash stores the output in the same directory as the contigs-db as "
                        "<contigs-db basename>-ANTISMASH. If you run antiSMASH on several contig-dbs using an "
                        "external-genomes file, anvi-run-antismash creates a directory for each contig-db "
                        "automatically."}))
    groupC.add_argument(*anvio.A('antismash-regions-output'), **anvio.K('antismash-regions-output'))
    groupC.add_argument(*anvio.A('antismash-work-dir'), **anvio.K('antismash-work-dir'))

    groupD = parser.add_argument_group("RESOURCES")
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
