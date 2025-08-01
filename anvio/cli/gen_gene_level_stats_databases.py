#!/usr/bin/env python
# -*- coding: utf-8
"""A program to generate genes databases all at once"""

import os
import sys

from argparse import Namespace

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.ccollections as ccollections

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.genelevelcoverages import TableForGeneLevelCoverages


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['profile-db', 'contigs-db', 'collection', 'bin',]
__provides__ = ['genes-db',]
__description__ = ("A program to compute genes databases for a ginen set of bins stored in "
                   "an anvi'o collection. Genes databases store gene-level coverage and "
                   "detection statistics, and they are usually computed and generated "
                   "automatically when they are required (such as running anvi-interactive "
                   "with `--gene-mode` flag). This program allows you to pre-compute them "
                   "if you don't want them to be done all at once")


run = terminal.Run()
progress = terminal.Progress()


def main():
    args = get_args()
    parameters = {'min_cov_for_detection': 0,
                  'outliers_threshold': args.outliers_threshold,
                  'zeros_are_outliers': args.zeros_are_outliers}

    try:
        splits_dict = ccollections.GetSplitNamesInBins(args).get_dict()
        bin_names = sorted(list(splits_dict.keys()))

        genes_dbs_found = [bin_name for bin_name in splits_dict if os.path.exists(utils.get_genes_database_path_for_bin(args.profile_db, args.collection_name, bin_name))]
        if genes_dbs_found:
            if args.just_do_it:
                # well then.
                pass
            else:
                raise ConfigError("It seems gene databases for some of the bins were already generated. As anvi'o will "
                                  "go through these bin names, it will makes sure those that already exist are clean "
                                  "and can be kept because the parameters that were used to generate them are identical "
                                  "to parameters set now. If they don't anvi'o will remove them and create new ones. But "
                                  "since you may have stored important information into these databases (such as states or "
                                  "collections), anvi'o wants you to give the chance to reconsider. If you don't care about "
                                  "them, use the `--just-do-it` flag, and anvi'o will go ahead and destroy anything that "
                                  "comes between you and your analyses.")


        run.info('Collection', args.collection_name)
        run.info('Bins to work with', ', '.join(bin_names))

        for parameter in parameters:
            run.info(parameter, parameters[parameter])

        for i in range(0, len(bin_names)):
            bin_name = bin_names[i]
            # recover the expected path
            genes_db_path = utils.get_genes_database_path_for_bin(args.profile_db, args.collection_name, bin_name)

            run.warning(None, header="WORKING ON '%s' (%d of %d)" % (bin_name, i + 1, len(bin_names)), lc='green')

            # if there is already one, we will see whether we want to keep it or remove it.
            if os.path.exists(genes_db_path):
                genes_db_is_clean = True
                try:
                    if args.inseq_stats:
                        # update the parameters dict to make sure gene-level stats databases are
                        # recognized distinctly
                        parameters['mode'] = 'INSEQ'
                        run.info('Mode', 'INSEQ', mc="red", nl_after=1)
                    else:
                        parameters['mode'] = 'STANDARD'
                        run.info('Mode', 'STANDARD', mc="red", nl_after=1)

                    TableForGeneLevelCoverages(genes_db_path, parameters, mode=parameters['mode'], run=terminal.Run(verbose=False)).check_params()
                    TableForGeneLevelCoverages(genes_db_path, parameters, mode=parameters['mode'], split_names=splits_dict[bin_name], run=terminal.Run(verbose=False)).check_split_names()
                except ConfigError as e:
                    run.info_single("A genes database found, but the parameters or split names associated with that "
                                    "database do not match to what we want now :( Anvi'o will regenerate it. This was "
                                    "the actual error, by the way: '%s'." % e.clear_text(), nl_after=1)

                    os.remove(genes_db_path)
                    genes_db_is_clean = False

                # if the genes db is clean, we don't need ot do anything
                if genes_db_is_clean:
                    run.info_single('A clean genes database found, moving on.', nl_after=1)
                    continue

            # we don't have a genes database, so we will create one.
            profile = dbops.ProfileSuperclass(Namespace(profile_db=args.profile_db,
                                                        contigs_db=args.contigs_db,
                                                        collection_name=args.collection_name,
                                                        bin_id=bin_name,
                                                        inseq_stats=args.inseq_stats,),
                                              r=terminal.Run(verbose=False))

            profile.init_gene_level_coverage_stats_dicts(outliers_threshold=args.outliers_threshold,
                                                         zeros_are_outliers=args.zeros_are_outliers)

            run.info('Genes database for %s' % bin_name, genes_db_path)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT DATABASES', "Which anvi'o databases do you wish to work today?")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))

    groupB = parser.add_argument_group('BIN(S) AND COLLECTION', "You can select a bin, multiple bins, or you can simply\
                                        focus on every bin in a collection by providing only a collection name. Once you\
                                        are done with your selection, anvi'o will generate an individual genes database\
                                        for each of the bin it finds.")
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupB.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupC = parser.add_argument_group('ADDITIONAL PARAMETERS', "These parameters are those that are critical to identify\
                                        outlier nucleotide positions and how to define what should be included in those\
                                        calculations. In most cases you can leave them as is, and things are going to be\
                                        alright.")
    groupC.add_argument(*anvio.A('zeros-are-outliers'), **anvio.K('zeros-are-outliers'))
    groupC.add_argument(*anvio.A('outliers-threshold'), **anvio.K('outliers-threshold'))

    groupD = parser.add_argument_group('PARAMETERS OF CONVENIENCE', "They say they save lives.")
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    groupE = parser.add_argument_group('INSEQ DATA', "When analyzing INSeq/Tn-Seq data")
    groupE.add_argument(*anvio.A('inseq-stats'), **anvio.K('inseq-stats'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
