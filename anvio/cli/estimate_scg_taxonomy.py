#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal
import anvio.taxonomyops.scg as scgtaxonomyops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'qclayssen']
__requires__ = ['profile-db', 'contigs-db', 'scgs-taxonomy', 'collection', 'bin', 'metagenomes']
__provides__ = ['genome-taxonomy', 'genome-taxonomy-txt']
__resources__ = [("Usage examples and warnings", "http://merenlab.org/scg-taxonomy")]
__description__ = ("Estimates taxonomy at genome and metagenome level. This program is the entry point to "
                   "estimate taxonomy for a given set of contigs (i.e., all contigs in a contigs database, "
                   "or contigs described in collections as bins). For this, it uses single-copy core "
                   "gene sequences and the GTDB database")


@terminal.time_program
def main():
    args = get_args()

    try:
        if args.metagenomes or args.external_genomes:
            t = scgtaxonomyops.SCGTaxonomyEstimatorMulti(args)
        else:
            t = scgtaxonomyops.SCGTaxonomyEstimatorSingle(args)

        t.estimate()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT #1', "The minimum you must provide this program is a contigs database. In which case\
                                                    anvi'o will attempt to estimate taxonomy for all contigs in it, assuming that\
                                                    the contigs database represents a single genome. If the contigs database is actually\
                                                    a metagenome, you should use the `--metagenome` flag to explicitly declare that.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('metagenome-mode'), **anvio.K('metagenome-mode', {'help': "Treat a given contigs database as an \
                                                    unbinned metagenome rather than treating it as a single genome."}))

    groupB = parser.add_argument_group('INPUT #2', "In addition, you can also point out a profile database. In which case you also must\
                                                    provide a collection name. When you do that anvi'o will offer taxonomy estimates for\
                                                    each bin in your collection.")
    groupB.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))

    groupC = parser.add_argument_group('INPUT #3', "You can also work with a external genomes or metagenomes file, assuming that you\
                                                    have multiple (meta)genomes with or without associated mapping results, and anvi'o would\
                                                    generate a single output file for all.")
    groupC.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupC.add_argument(*anvio.A('metagenomes'), **anvio.K('metagenomes'))

    groupD = parser.add_argument_group('OUTPUT AND FORMATTING', "Anvi'o will do its best to offer you some fancy output tables for your viewing \
                                                                 pleasure by default. But in addition to that, you can ask the resulting information \
                                                                 to be stored in a TAB-delimited file (which is a much better way to include the \
                                                                 results in your study as supplementary information, or work with these results \
                                                                 using other analysis tools such as R). Depending on the mode you are running this \
                                                                 program, anvi'o may ask you to use an 'output file prefix' rather than an 'output \
                                                                 file path'.")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupD.add_argument(*anvio.A('per-scg-output-file'), **anvio.K('per-scg-output-file'))
    groupD.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    groupD.add_argument(*anvio.A('taxonomic-level'), **anvio.K('taxonomic-level', {'default': None}))
    groupD.add_argument(*anvio.A('matrix-format'), **anvio.K('matrix-format', {'help': "If you want the reports to look like sparse matrices whenever "
                                                                                   "possible, declare this flag. Matrices are especially good to use "
                                                                                   "when you are working with internal/external genomes since they can "
                                                                                   "show you quickly the distribution of each taxon across all metagenomes "
                                                                                   "in programs like EXCEL. WELL TRY IT AND SEE."}))
    groupD.add_argument(*anvio.A('raw-output'), **anvio.K('raw-output'))

    groupE = parser.add_argument_group('PERFORMANCE', "We are not sure if allocating more threads for this operation will change anything.\
                                                       But hey. One can try.")
    groupE.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    groupE = parser.add_argument_group('AUTHORITY', "Assert your dominance.")
    groupE.add_argument(*anvio.A('scg-name-for-metagenome-mode'), **anvio.K('scg-name-for-metagenome-mode'))
    groupE.add_argument(*anvio.A('report-scg-sequences-file-prefix'), **anvio.K('report-scg-sequences-file-prefix'))
    groupE.add_argument(*anvio.A('report-scg-frequencies'), **anvio.K('report-scg-frequencies'))
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    groupF = parser.add_argument_group('ADVANCED', "Very pro-like stuff.")
    groupF.add_argument(*anvio.A('simplify-taxonomy-information'), **anvio.K('simplify-taxonomy-information'))
    groupF.add_argument(*anvio.A('compute-scg-coverages'), **anvio.K('compute-scg-coverages'))
    groupF.add_argument(*anvio.A('update-profile-db-with-taxonomy'), **anvio.K('update-profile-db-with-taxonomy'))

    groupG = parser.add_argument_group('BORING', "Options that you will likely never need.")
    groupG.add_argument(*anvio.A('taxonomy-database'), **anvio.K('taxonomy-database'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
