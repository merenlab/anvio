#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.kegg as kegg

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva']
__requires__ = ["contigs-db", "kegg-data", "kegg-functions", "profile-db", "collection", "bin",
                "external-genomes", "internal-genomes", "metagenomes", "user-modules-data",
                "enzymes-txt", "pan-db", "genomes-storage-db"]
__provides__ = ["kegg-metabolism","user-metabolism"]
__description__ = "Reconstructs metabolic pathways and estimates pathway completeness for a given set of contigs"


@time_program
def main():
    args = get_args()

    try:
        if args.metagenomes or args.external_genomes or args.internal_genomes:
            m = kegg.KeggMetabolismEstimatorMulti(args)
        else:
            m = kegg.KeggMetabolismEstimator(args)

        if args.list_available_modes:
            m.list_output_modes()
            sys.exit()
        elif args.list_available_output_headers:
            m.list_output_headers()
            sys.exit()
        else:
            m.estimate_metabolism()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupD = parser.add_argument_group('METABOLISM DATA LOCATION', "What modules database(s) do you want to use?")
    groupD.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir'))
    groupD.add_argument(*anvio.A('user-modules'), **anvio.K('user-modules', {'help': "Directory location where user-defined module "
                                                                      "files (and database) are kept. You must have previously run "
                                                                      "`anvi-setup-user-modules` on this directory for this to work."}))
    groupD.add_argument(*anvio.A('only-user-modules'), **anvio.K('only-user-modules'))

    groupI = parser.add_argument_group('INPUT #1 - ESTIMATION ON SINGLE GENOMES OR METAGENOMES',
                                                   "The minimum you must provide this program is a contigs database. In which case "
                                                   "anvi'o will attempt to estimate metabolism for all contigs in it, assuming that "
                                                   "the contigs database represents a single genome. If the contigs database is actually "
                                                   "an unbinned metagenome and you want per-contig estimates instead, you can use the "
                                                   "`--metagenome` flag to explicitly declare that. It is also acceptable to run this "
                                                   "program on metagenomes without using metagenome mode if you want community-level "
                                                   "metabolism estimates (ie, combining information across all contigs from different "
                                                   "populations in the community).")
    groupI.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupI.add_argument(*anvio.A('metagenome-mode'), **anvio.K('metagenome-mode'))

    groupP = parser.add_argument_group('INPUT #2 - ESTIMATION ON BINS', "If you also provide a profile database AND a collection name, anvi'o will "
                                                   "estimate metabolism separately for each bin in your collection. You can also limit "
                                                   "those estimates to a specific bin or set of bins in the collection using the parameters "
                                                   "`--bin-id` or `--bin-ids-file`, respectively.")
    groupP.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupP.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupP.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupP.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupM = parser.add_argument_group('INPUT #3 - MULTI-MODE', "If you have multiple contigs databases to work with, you can put them all into a file. "
                                                   "Then anvi'o will run estimation separately on each database and generate a single output file for all. "
                                                   "There are 3 types of input files to choose from depending on whether you want to treat your DBs as "
                                                   "single genomes (external), genomes in collections (internal), or unbinned metagenomes (for contig-level "
                                                   "estimates, i.e. 'metagenome mode') in your contigs DBs.")
    groupM.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupM.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupM.add_argument(*anvio.A('metagenomes'), **anvio.K('metagenomes', {'help': "Same format as an external genomes file, but choosing this option "
                                                                                "ensures each DB in the input file is analyzed "
                                                                                "with 'metagenome mode' for per-contig estimates."}))

    groupX = parser.add_argument_group('INPUT #4 - ESTIMATION ON A LIST OF ENZYMES',
                                                   "If all you have is a list of enzymes, you can use them to estimate metabolism "
                                                   "by creating an enzymes-txt file which describes at least 1) some sort of unique identifier "
                                                   "for the associated gene (gene name, gene caller id, etc), 2) the enzyme accession, and 3) "
                                                   "the annotation source. Provide this file as input, and this program will pretend that "
                                                   "all enzymes in the file come from a single genome and will give back to you a summary of the "
                                                   "metabolic capabilities encoded by this list of enzymes.")
    groupX.add_argument(*anvio.A('enzymes-txt'), **anvio.K('enzymes-txt', {'required': False}))

    groupN = parser.add_argument_group('INPUT #5 - ESTIMATION ON A PANGENOME', 
                                                   "If you have a collection of binned gene clusters in a pangenome, you can estimate "
                                                   "metabolism on each bin. You will need to provide the collection name (and optionally, "
                                                   "bin name(s)) -- see INPUT #2 section above. This program will use the most common "
                                                   "annotation (for each annotation source) from each gene cluster. Note that this input "
                                                   "option is not compatible with the `--add-copy-number` or `--add-coverage` flags.")
    groupN.add_argument('--pan-db', **anvio.K('pan-db', {'required': False}))
    groupN.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))

    groupC = parser.add_argument_group('OUTPUT - GENERAL OPTIONS', "Parameters for controlling estimation output of any type. The output will be "
                                                 "TAB-delimited files which by default are prefixed with 'kegg-metabolism', but you can of course "
                                                 "change that name here.")
    groupC.add_argument(*anvio.A('module-completion-threshold'), **anvio.K('module-completion-threshold'))
    groupC.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    groupC.add_argument(*anvio.A('include-zeros'), **anvio.K('include-zeros'))
    groupC.add_argument(*anvio.A('only-complete'), **anvio.K('only-complete'))
    groupC.add_argument(*anvio.A('add-copy-number'), **anvio.K('add-copy-number'))
    groupC.add_argument(*anvio.A('include-kos-not-in-kofam'), **anvio.K('include-kos-not-in-kofam'))
    groupC.add_argument(*anvio.A('include-stray-KOs'), **anvio.K('include-stray-KOs', {'help': "'Stray KOs' are what we call KEGG Orthlogs "
                                                            "that KEGG does not provide a bit score threshold for. Anvi'o can estimate "
                                                            "thresholds (and sometimes updates the HMMs) for these KOs, and you can annotate "
                                                            "them in your data if you run `anvi-run-kegg-kofams` with the `--include-stray-KOs` "
                                                            "flag. If you did that, and you now want to include those annotations in metabolism "
                                                            "estimates, then you should apply the same flag here. This is kind of similar to the "
                                                            "`--include-kos-not-in-kofam` flag, but it applies only to the stray KOs that you can "
                                                            "annotate with `anvi-run-kegg-kofams --include-stray-KOs`."}))
    groupC.add_argument(*anvio.A('ignore-unknown-KOs'), **anvio.K('ignore-unknown-KOs'))
    groupC.add_argument(*anvio.A('exclude-dashed-reactions'), **anvio.K('exclude-dashed-reactions'))

    groupL = parser.add_argument_group('OUTPUT - LONG-FORMAT OPTIONS', "Parameters for controlling long-format output (the default).")
    groupL.add_argument(*anvio.A('output-modes'), **anvio.K('output-modes'))
    groupL.add_argument(*anvio.A('list-available-modes'), **anvio.K('list-available-modes'))
    groupL.add_argument(*anvio.A('custom-output-headers'), **anvio.K('custom-output-headers'))
    groupL.add_argument(*anvio.A('list-available-output-headers'), **anvio.K('list-available-output-headers'))
    groupL.add_argument(*anvio.A('add-coverage'), **anvio.K('add-coverage'))

    groupM = parser.add_argument_group('OUTPUT - MATRIX OPTIONS', "Parameters for controlling matrix output. Use --matrix-format to request this "
                                                    "type of output.")
    groupM.add_argument(*anvio.A('matrix-format'), **anvio.K('matrix-format', {'help': "If you want to generate the output in several sparse matrices instead "
                                                                                   "of one file, use this flag. In each matrix, contigs DBs will be arranged in "
                                                                                   "columns and KEGG modules in rows. This output option is especially "
                                                                                   "appropriate for input option #3."}))
    groupM.add_argument(*anvio.A('include-metadata'), **anvio.K('include-metadata'))
    groupM.add_argument(*anvio.A('module-specific-matrices'), **anvio.K('module-specific-matrices'))
    groupM.add_argument(*anvio.A('no-comments'), **anvio.K('no-comments'))

    groupD = parser.add_argument_group('DEBUG', "Parameters to use if you think something fishy is going on or otherwise want to exert more control. Go for it.")
    groupD.add_argument(*anvio.A('get-raw-data-as-json'), **anvio.K('get-raw-data-as-json'))
    groupD.add_argument(*anvio.A('store-json-without-estimation'), **anvio.K('store-json-without-estimation'))
    groupD.add_argument(*anvio.A('estimate-from-json'), **anvio.K('estimate-from-json'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
