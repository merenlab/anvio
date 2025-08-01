#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.tables.miscdata as miscdata

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mahmoudyousef98']
__requires__ = ["pan-db", "genomes-storage-db"]
__provides__ = []
__resources__ = [("The role of gene cluster homogeneity described in the Anvi'o pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters")]
__description__ = "Compute homogeneity for gene clusters"


@terminal.time_program
def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    if args.gene_cluster_id and args.gene_cluster_ids_file:
        raise ConfigError('You should either declare a single gene cluster name, or gene cluster names in a file')

    if (args.gene_cluster_id or args.gene_cluster_ids_file) and args.collection_name:
        raise ConfigError('You can either declare specific list of gene clusters to work with (through `--gene-cluster-id` or `--gene-cluster-ids-file`) or '
                          'go the collection way using parameters `--collection-name` and `--bin-name`. Those are not to be '
                          'mixed.')

    if not args.output_file and not args.store_in_db:
        if args.just_do_it:
            run.warning("Lol. You are making anvi'o compute stuff that will not be stored anywhere. Weird you :(", lc="green")
        else:
            raise ConfigError("But why no select a way to report this stuff? :/ You can ask anvi'o to store your results as a "
                              "a TAB-delimited file. Or add them to directly to the database. If you insiste that you just want "
                              "to run this analysis withour really storing anything anywhere you can always use `--just-do-it` "
                              "(because why not, it is your computer).")

    if args.output_file:
        filesnpaths.is_output_file_writable(args.output_file)

    if args.store_in_db and (args.gene_cluster_id or args.gene_cluster_ids_file or args.collection_name):
        if args.just_do_it:
            run.warning("You asked anvi'o to store your results into the database, but you are not running this analysis for all "
                        "gene clusters available in your pan database. Which is OK, but you will likely end up having some gene "
                        "clusters without any homogeneity estimates. WE DON'T CARE.", lc="green")
        else:
            raise ConfigError("You ask anvi'o to store results into the database, but you are also specifying a set of gene clusters "
                              "by using parameters listed under SELECTION (see the help menu for a refresher). Anvi'o is reluctant "
                              "to do this, becasue this may end up overwriting gene cluster homogeneity data in the database with "
                              "a small number of gene clusters. On the other hand you may know what you are doing, in which case "
                              "you should use `--just-do-it`, and anvi'o will be like 'K, sure'.")

    gene_cluster_ids = set([])
    if args.collection_name:
        progress.new('Initializing gene clusters')
        progress.update('...')
        pan = summarizer.PanSummarizer(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        progress.end()

        if not args.bin_id:
            raise ConfigError("When you use a collection name, you must also declare a bin id :/")

        pan.collections.is_bin_in_collection(collection_name=args.collection_name, bin_name=args.bin_id)
        collection_dict = pan.collections.get_collection_dict(args.collection_name)
        gene_cluster_ids = set(collection_dict[args.bin_id])

        run.info('Mode', 'Reporting homogeneity for gene clusters in the collection %s and bin %s.' % (args.collection_name, args.bin_id))
    elif (args.gene_cluster_id or args.gene_cluster_ids_file):
        if args.gene_cluster_id:
            gene_cluster_ids = set([args.gene_cluster_id])
            run.info('Mode', 'Reporting homogeneity for a single gene cluster.')
        else:
            columns = utils.get_columns_of_TAB_delim_file(args.gene_cluster_ids_file, include_first_column=True)
            if len(columns) != 1:
                raise ConfigError("The input file for gene cluster IDs must contain a single column. It seems yours has %d :/" % len(columns))

            gene_cluster_ids = set([p.strip('\n') for p in open(args.gene_cluster_ids_file).readlines()])
            run.info('Mode', 'Reporting homogeneity for a list of gene clusters from an input file.')
    else:
        pan = dbops.PanSuperclass(args, r=terminal.Run(verbose=False))
        gene_cluster_ids = pan.gene_cluster_names

        run.warning("You did not specify which gene clusters to work with, so anvi'o will work with all %d gene "
                    "clusters found in the pan database." % len(gene_cluster_ids))

    pan = dbops.PanSuperclass(args)
    pan.init_gene_clusters(gene_cluster_ids)
    d = pan.compute_homogeneity_indices_for_gene_clusters(gene_cluster_ids, num_threads=args.num_threads)

    if d is None:
        raise ConfigError('The homogeneity index algorithm has probably alerted you that an error has occured. Please '
                          'refer to that error to understand why this failed')

    if args.store_in_db:
        args.just_do_it = True
        miscdata.TableForItemAdditionalData(args).add(d, ['functional_homogeneity_index', 'geometric_homogeneity_index', 'combined_homogeneity_index'])

    if args.output_file:
        utils.store_dict_as_TAB_delimited_file(d, args.output_file)
        run.info("Output file", args.output_file, mc="green")


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES', "Input files from the pangenome analysis.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))

    groupB = parser.add_argument_group('REPORTING', "How do you want results to be reported? Anvi'o can produce a TAB-delimited output file for\
                                                     you (for which you would have to provide an output file name). Or the results can be stored\
                                                     in the pan database directly, for which you would have to explicitly ask for it. You can get\
                                                     both as well in case you are a fan of redundancy and poor data analysis practices. Anvi'o\
                                                     does not judge.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupB.add_argument(*anvio.A('store-in-db'), **anvio.K('store-in-db'))

    groupC = parser.add_argument_group('SELECTION', "Which gene clusters should be analyzed. You can ask for a single gene cluster,\
                                       or multiple ones listed in a file, or you can use a collection and bin name to list gene clusters\
                                       of interest.")
    groupC.add_argument(*anvio.A('gene-cluster-id'), **anvio.K('gene-cluster-id'))
    groupC.add_argument(*anvio.A('gene-cluster-ids-file'), **anvio.K('gene-cluster-ids-file'))
    groupC.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupC.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupD = parser.add_argument_group('OPTIONAL',"Optional stuff available for you to use")
    groupD.add_argument('--quick-homogeneity', default=False, action='store_true', help="By default, anvi'o will use a homogeneity\
                        algorithm that checks for horizontal and vertical geometric homogeneity (along with functional). With this\
                        flag, you can tell anvi'o to skip horizontal geometric homogeneity calculations. It will be less accurate but quicker.")
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
