#!/usr/bin/env python
# -*- coding: utf-8
"""Generate summary output from anvi'o datbases"""

import os
import sys

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['profile-db', 'contigs-db', 'collection', 'pan-db', 'genomes-storage-db',]
__provides__ = ['summary']
__description__ = ("Summarizer for anvi'o pan or profile db's. Essentially, this program takes a collection id along with either "
                   "a profile database and a contigs database or a pan database and a genomes storage and generates a static HTML "
                   "output for what is described in a given collection. The output directory will contain almost everything any "
                   "downstream analysis may need, and can be displayed using a browser without the need for an anvi'o installation. "
                   "For this reason alone, reporting summary outputs as supplementary data with publications is a great idea for "
                   "transparency and reproducibility")
__resources__ = [("anvi-summarize in the metagenomic workflow tutorial", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-summarize"), ("anvi-summarize in the pangenomic workflow tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/#summarizing-an-anvio-pan-genome")]


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()

    # bad code but a good check early on. Due to multiple inheriteance downstream, there is actually no good
    # place to perform this check before its too much CPU time has already passed (see #1430 for an actual
    # user complaint):
    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    if A('output_dir'):
        filesnpaths.check_output_directory(A('output_dir'))

    utils.is_all_npm_packages_installed()

    # k, gud. lets move on with business.
    db_type = utils.get_db_type(args.pan_or_profile_db)

    if db_type == 'pan':
        args.pan_db = args.pan_or_profile_db
        summary = summarizer.PanSummarizer(args)
    elif db_type == 'profile':
        args.profile_db = args.pan_or_profile_db
        profile_db = dbops.ProfileDatabase(args.profile_db)

        if not profile_db.meta['contigs_db_hash']:
            raise ConfigError("This profile database is not affiliated with any contigs database. Anvi'o currently does not "
                              "summarize such stuff :/")
        profile_db.disconnect()

        if not args.contigs_db:
            raise ConfigError("You must provide a contigs database when you summarize anvi'o profiles. True story.")

        # set this global arg so that all functional annotations are reflected in the summary output
        anvio.RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE = True

        summary = summarizer.ProfileSummarizer(args)
    else:
        raise ConfigError("Well. '%s' is neither an anvi'o pan database, nor an anvi'o profile database. There is nothing this "
                           "program can't do for you if you feed it with the right stuff. Just sayin'" % args.pan_or_profile_db)

    summary.process()


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('PROFILE', "The profile. It could be a standard or pan profile database.")
    groupB = parser.add_argument_group('PROFILE TYPE SPECIFIC PARAMETERS', "If you are summarizing a collection stored in\
                                        a standard anvi'o profile, you will need a contigs database to go with it. If you\
                                        are working with a pan profile, then you will need to provide a genomes storage.\
                                        Don't worry too much, because anvi'o will warn you gently if you make a mistake.")
    groupC = parser.add_argument_group('STANDARD PROFILE SPECIFIC PARAMS', "Parameters that are only relevant to standard\
                                        profile summaries (declaring or not declaring them will not change anything if you\
                                        are summarizing a pan profile).")
    groupD = parser.add_argument_group('PAN PROFILE SPECIFIC PARAMS', "Parameters that are only relevant to pan\
                                        profile summaries (declaring or not declaring them will not change anything if you\
                                        are summarizing a standard profile).")
    groupE = parser.add_argument_group('COMMONS', "Common parameters for both pan and standard profile summaries.")
    groupF = parser.add_argument_group('EXTRA', "Extra stuff because you're extra.")

    groupA.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))

    groupB.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupB.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))

    groupC.add_argument(*anvio.A('init-gene-coverages'), **anvio.K('init-gene-coverages'))
    groupC.add_argument(*anvio.A('calculate-Q2Q3-carefully'), **anvio.K('calculate-Q2Q3-carefully'))
    groupC.add_argument(*anvio.A('reformat-contig-names'), **anvio.K('reformat-contig-names'))
    groupC.add_argument('--report-aa-seqs-for-gene-calls', default=False, action='store_true', help="You can use this flag if\
                                  you would like amino acid AND dna sequences for your gene calls in the genes output\
                                  file. By default, only dna sequences are reported.")

    groupD.add_argument(
        *anvio.A('report-DNA-sequences'),
        **anvio.K('report-DNA-sequences', {'help': anvio.K('report-DNA-sequences')['help'] + ' Also note, since gene clusters are \
                                                   aligned via amino acid sequences, using this flag removes alignment information \
                                                   manifesting in the form of gap characters (`-` characters) that would \
                                                   be present if amino acid sequences were reported. Read the warnings during \
                                                   runtime for more detailed information.'})
    )

    groupE.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupE.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupE.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupE.add_argument(*anvio.A('cog-data-dir'), **anvio.K('cog-data-dir'))
    groupE.add_argument(*anvio.A('quick-summary'), **anvio.K('quick-summary'))

    groupF.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))


    return parser.get_args(parser)


if __name__ == '__main__':
    main()
