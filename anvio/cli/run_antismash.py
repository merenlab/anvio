#!/usr/bin/env python

import sys

import anvio
import anvio.antismash as antismash

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']
__requires__ = ['contigs-db', 'external-genomes']
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


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("INPUT", "The contigs database(s) to annotate. Provide either a single "
                             "one with -c, or many at once through an external-genomes file with -e.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': False}))

    groupB = parser.add_argument_group("ANTISMASH", "Options forwarded to antiSMASH.")
    groupB.add_argument("--antismash-dir", default=None,
                        help="The directory antiSMASH is installed in (its conda environment), e.g. "
                             "~/miniforge3/envs/antismash. You usually do NOT need this: anvi'o first looks for "
                             "'antismash' on your PATH, and failing that automatically searches your conda "
                             "environments for it. If anvi'o cannot find antiSMASH on its own, or if you have "
                             "several antiSMASH installations and would like to point to a specific one you can "
                             "use this flag.")
    groupB.add_argument("--antismash-data-dir", default=None,
                        help="This points to the antiSMASH database and you only need this flag if you did not "
                             "store the database in the default directory, in which case the program does not know "
                             "where to look.")
    groupB.add_argument("--taxon", default="bacteria",
                        choices=["bacteria", "fungi"],
                        help="Source organism kingdom. Default: bacteria. You can change it to fungi.")
    groupB.add_argument("--min-contig-length", type=int, default=3000,
                        help="antiSMASH uses a default contig length of min 3000 bp and filters out shorter contigs "
                             "before BGCs. You can manually change that with this flag.")
    groupB.add_argument("--antismash-add-on", default="",
                        help="Using this flag you can ask anvi-run-antismash to include additional analysis "
                             "available in antiSMASH, including --fullhmmer, --clusterhmmer, --tigrfam, --pfam2go, "
                             "--tfbs and --cassis.")

    groupC = parser.add_argument_group("OUTPUT", "anvi-run-antismash by default annotates the contigs-db with "
                             "smCOGs and identifies BGC clusters (antiSMASH genes and antiSMASH regions) and also "
                             "saves a browsable HTML report plus a per-region summary. You can use this flag:")
    groupC.add_argument("--include-detailed-output", default=False, action="store_true",
                        help="To include ClusterBlast and MIBiG analysis, smCOG trees as well as the per-region "
                             "GenBank (.gbk) files (useful as input to tools like BiG-SCAPE), the raw JSON, and "
                             "the comparison tabs. This slows down the analysis and generates quite a few files, "
                             "but gives you a more comprehensive output of which we can't judge if it is useful "
                             "to you or not.")
    groupC.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'help': "Used to define where to save "
                        "the antiSMASH report and summary files. It is an optional flag and if undefined, "
                        "anvi-run-antismash stores the output in the same directory as the contigs-db as "
                        "<contigs-db basename>-ANTISMASH. If you run antiSMASH on several contig-dbs using an "
                        "external-genomes file, anvi-run-antismash creates a directory for each contig-db "
                        "automatically."}))
    groupC.add_argument("--regions-output", default=None,
                        help="Can be used to set the path for the per-region dataframe. The default is to include "
                             "it as antismash_regions.txt inside the report directory. This is only relevant if you "
                             "run anvi-run-antismash on a single database, if you use an external-genomes.txt "
                             "anvi-run-antismash creates individual directories automatically.")
    groupC.add_argument("--work-dir", default=None,
                        help="You can use this directory to store the intermediates instead of a tempdir "
                             "(helpful for debugging). This is only relevant if you run anvi-run-antismash on a "
                             "single database.")

    groupD = parser.add_argument_group("RESOURCES")
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
