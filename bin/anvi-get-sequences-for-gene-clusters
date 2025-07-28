#!/usr/bin/env python
# -*- coding: utf-8
"""Export aligned sequences from anvi'o pan genomes"""

import os
import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ["Ryan Bartelme", "Jay Osvatic"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db', 'genomes-storage-db']
__provides__ = ['genes-fasta', 'concatenated-gene-alignment-fasta', 'misc-data-items']
__resources__ = [("In action in the Anvi'o pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/#scrutinizing-phylogenomics")]
__description__ = "Do cool stuff with gene clusters in anvi'o pan genomes"


@terminal.time_program
def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except DictIOError as e:
        print(e)
        sys.exit(-3)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    sanity_check(args)

    gene_cluster_ids = set([])

    if args.collection_name or args.list_collections or args.list_bins:
        progress.new('Initializing')
        progress.update('...')
        pan = summarizer.PanSummarizer(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        progress.end()

        if not args.bin_id:
            raise ConfigError("When you use a collection name, you must also declare a bin id :/ You don't know what bin you "
                              "want? Use the flag `--list-bins`.")

        pan.collections.is_bin_in_collection(collection_name=args.collection_name, bin_name=args.bin_id)
        collection_dict = pan.collections.get_collection_dict(args.collection_name)
        gene_cluster_ids = set(collection_dict[args.bin_id])

        run.info('Mode', 'Working with gene clusters in collection %s and bin %s.' % (args.collection_name, args.bin_id))
    elif (args.gene_cluster_id or args.gene_cluster_ids_file):
        if args.gene_cluster_id:
            gene_cluster_ids = set([args.gene_cluster_id])
            run.info('Mode', 'Working with a single gene cluster.')
        else:
            columns = utils.get_columns_of_TAB_delim_file(args.gene_cluster_ids_file, include_first_column=True)
            if len(columns) != 1:
                raise ConfigError("The input file for gene cluster IDs must contain a single column. It seems yours has %d :/" % len(columns))

            gene_cluster_ids = set([p.strip('\n') for p in open(args.gene_cluster_ids_file, 'r').readlines()])
            run.info('Mode', 'Reporting alignments for a list of gene clusters from an input file.')
    else:
        run.info('Mode', 'Working with all gene clusters.')

    pan = dbops.PanSuperclass(args, r=terminal.Run(verbose=False))
    pan.init_gene_clusters(gene_cluster_ids)
    pan.run = run

    # this will simply return pan.gene_clusters if there are no filters requested
    filtered_gene_clusters_dict = pan.filter_gene_clusters_dict(args)

    num_gene_clusters = len(filtered_gene_clusters_dict)
    num_genes = sum([len(gc.keys()) for gc in filtered_gene_clusters_dict.values()])
    run.warning("Your filters resulted in %d gene clusters that contain a total of %d genes. "
                "for downstream analyses. Just so you know." % (num_gene_clusters, num_genes), header="INFO", lc="green")

    if args.dry_run:
        run.info_single("Dry run, eh? Well. Your dry run is done!")
        sys.exit(0)

    if args.output_file:
        if args.concatenate_gene_clusters:
            pan.write_sequences_in_gene_clusters_for_phylogenomics(gene_clusters_dict=filtered_gene_clusters_dict,
                                                                   output_file_path=args.output_file,
                                                                   report_DNA_sequences=args.report_DNA_sequences,
                                                                   align_with=args.align_with,
                                                                   separator=args.separator,
                                                                   partition_file_path=args.partition_file)
        else:
            pan.write_sequences_in_gene_clusters_to_file(gene_clusters_dict=filtered_gene_clusters_dict,
                                                         output_file_path=args.output_file,
                                                         report_DNA_sequences=args.report_DNA_sequences)
    elif args.split_output_per_gene_cluster:
        if args.concatenate_gene_clusters:
            for gene_cluster_name in filtered_gene_clusters_dict:
                output_file_path = f"{args.output_file_prefix}{gene_cluster_name}.fa"
                pan.write_sequences_in_gene_clusters_for_phylogenomics(gene_clusters_dict={gene_cluster_name: filtered_gene_clusters_dict[gene_cluster_name]},
                                                                       output_file_path=output_file_path,
                                                                       report_DNA_sequences=args.report_DNA_sequences,
                                                                       align_with=args.align_with,
                                                                       separator=args.separator,
                                                                       partition_file_path=args.partition_file)
        else:
            for gene_cluster_name in filtered_gene_clusters_dict:
                output_file_path = f"{args.output_file_prefix}{gene_cluster_name}.fa"
                pan.write_sequences_in_gene_clusters_to_file(gene_clusters_dict={gene_cluster_name: filtered_gene_clusters_dict[gene_cluster_name]},
                                                             output_file_path=output_file_path,
                                                             report_DNA_sequences=args.report_DNA_sequences)


def sanity_check(args):
    if args.gene_cluster_id and args.gene_cluster_ids_file:
        raise ConfigError('You should either declare a single gene cluster name or a set of gene cluster names in a file, but not both :/')

    if args.gene_cluster_id and args.split_output_per_gene_cluster:
        raise ConfigError("If you are providing a single gene cluster id via `--gene-cluster-id`, then you don't "
                          "need to use the flag `--split-output-per-gene-cluster`.")

    if (args.gene_cluster_id or args.gene_cluster_ids_file) and args.collection_name:
        raise ConfigError('You can either declare specific list of gene clusters to work with (through `--gene-cluster-id` or `--gene-cluster-ids-file`) or '
                          'go the collection way using parameters `--collection-name` and `--bin-name`. Those are not to be '
                          'mixed. If you need to know what collections are available in the pan database, use the flag '
                          '`--list-collections`.')

    if args.output_file and args.split_output_per_gene_cluster:
        raise ConfigError("The parameter `--output-file` is incompatible with the flag `--split-output-per-gene-cluster`. Why? "
                          "Because when you declare `--split-output-per-gene-cluster`, anvi'o reports many output files and "
                          "assigns the name of each ")

    if args.split_output_per_gene_cluster and args.concatenate_gene_clusters and args.partition_file:
        raise ConfigError("You can't ask anvi'o to produce split output files per concatenated gene cluster, and offer an output for "
                          "a partition file. Why? Because anvi'o developers are lazy, and thought now one would ever need it really. "
                          "But if you managed ot run into this error, please let an anvi'o developer know online and they will "
                          "implement something for it.")

    if args.split_output_per_gene_cluster:
        if not args.output_file_prefix:
            args.output_file_prefix = os.getcwd() + '/'
        else:
            if args.output_file_prefix.startswith('/'):
                pass
            elif args.output_file_prefix.endswith('/'):
                pass
            elif args.output_file_prefix.find('/') > -1:
                args.output_file_prefix = os.path.abspath('/'.join(args.output_file_prefix.split('/')[:-1])) + '/' + args.output_file_prefix.split('/')[-1] + '_'
            else:
                args.output_file_prefix = os.getcwd() + '/' + args.output_file_prefix + '_'

        filesnpaths.is_output_dir_writable(os.path.dirname(args.output_file_prefix))
    else:
        if not args.output_file:
            args.output_file = os.path.join(os.path.abspath(os.getcwd()), 'GENE-CLUSTER-ALIGNMENTS.fa')

        filesnpaths.is_output_file_writable(args.output_file, ok_if_exists=False)

    if not args.concatenate_gene_clusters and args.align_with:
        raise ConfigError("If you are not asking for concatenated gene clusters, then you really don't need to specify an aligner. "
                          "Yes, anvi'o could ignore that flag and move on, but who would we become if we get used to ignoring those "
                          "little signs of 'I am not sure what I am doing'?")

    if not args.concatenate_gene_clusters and args.partition_file:
        raise ConfigError("Partition files are only relevant when you use the flag `--concatenate-gene-clusters`.")

    filesnpaths.is_output_file_writable(args.partition_file) if args.partition_file else None

    return args


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES', "Input files from the pangenome analysis.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))

    groupB = parser.add_argument_group('OUTPUT OPTION #1: SINGLE FILE', "Here you get to chose an output file name to report things. LUCKY YOU.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'metavar':'FASTA'}))

    groupC = parser.add_argument_group('OUTPUT OPTION #2: MULTIPLE FILES', "Here you can ask anvi'o to report multiple FASTA files, where each gene cluster "
                                       "goes into its own unique FASTA file. If you declare teh flag `--split-output-per-gene-cluster`, anvi'o "
                                       "will store each gene cluster into your work directory, using the gene cluster name as the output file name. "
                                       "You can use the `--output-file-prefix` parameter to store them in a different directory and/or with specific "
                                       "file name prefixes.")
    groupC.add_argument(*anvio.A('split-output-per-gene-cluster'), **anvio.K('split-output-per-gene-cluster'))
    groupC.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))

    groupD = parser.add_argument_group('SELECTION', "Which gene clusters should be reported. You can ask for a single gene cluster,\
                                       or multiple ones listed in a file, or you can use a collection and bin name to list gene clusters\
                                       of interest. If you give nothing, this program will export alignments for every single gene cluster\
                                       found in the profile database (and this is called 'customer service').")
    groupD.add_argument(*anvio.A('gene-cluster-id'), **anvio.K('gene-cluster-id'))
    groupD.add_argument(*anvio.A('gene-cluster-ids-file'), **anvio.K('gene-cluster-ids-file'))
    groupD.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupD.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupE = parser.add_argument_group('ADVANCED FILTERS', "If you are here you must be looking for ways to specify\
                                        exactly what you want from that database of gene clusters. These filters will\
                                        be applied to what your previous selections reported.")
    groupE.add_argument(*anvio.A('min-num-genomes-gene-cluster-occurs'), **anvio.K('min-num-genomes-gene-cluster-occurs'))
    groupE.add_argument(*anvio.A('max-num-genomes-gene-cluster-occurs'), **anvio.K('max-num-genomes-gene-cluster-occurs'))
    groupE.add_argument(*anvio.A('min-num-genes-from-each-genome'), **anvio.K('min-num-genes-from-each-genome'))
    groupE.add_argument(*anvio.A('max-num-genes-from-each-genome'), **anvio.K('max-num-genes-from-each-genome'))
    groupE.add_argument(*anvio.A('max-num-gene-clusters-missing-from-genome'), **anvio.K('max-num-gene-clusters-missing-from-genome'))
    groupE.add_argument(*anvio.A('min-functional-homogeneity-index'), **anvio.K('min-functional-homogeneity-index'))
    groupE.add_argument(*anvio.A('max-functional-homogeneity-index'), **anvio.K('max-functional-homogeneity-index'))
    groupE.add_argument(*anvio.A('min-geometric-homogeneity-index'), **anvio.K('min-geometric-homogeneity-index'))
    groupE.add_argument(*anvio.A('max-geometric-homogeneity-index'), **anvio.K('max-geometric-homogeneity-index'))
    groupE.add_argument(*anvio.A('min-combined-homogeneity-index'), **anvio.K('min-combined-homogeneity-index'))
    groupE.add_argument(*anvio.A('max-combined-homogeneity-index'), **anvio.K('max-combined-homogeneity-index'))
    groupE.add_argument(*anvio.A('add-into-items-additional-data-table'), **anvio.K('add-into-items-additional-data-table'))

    groupF = parser.add_argument_group('OTHER STUFF', "Yes. Stuff that are not like the ones above.")
    groupF.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupF.add_argument(*anvio.A('list-bins'), **anvio.K('list-bins'))

    groupG = parser.add_argument_group('PHYLOGENOMICS', "Get separately aligned and concatenated sequences for phylogenomics.")
    groupG.add_argument(*anvio.A('concatenate-gene-clusters'), **anvio.K('concatenate-gene-clusters'))
    groupG.add_argument(*anvio.A('partition-file'), **anvio.K('partition-file'))
    groupG.add_argument(*anvio.A('separator'), **anvio.K('separator'))
    groupG.add_argument(*anvio.A('align-with'), **anvio.K('align-with'))
    groupG.add_argument(*anvio.A('list-aligners'), **anvio.K('list-aligners'))

    groupH = parser.add_argument_group('LIFE SAVERS', "Just when you need them.")
    groupH.add_argument(*anvio.A('report-DNA-sequences'), **anvio.K('report-DNA-sequences'))
    groupH.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupH.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))


    return parser.get_args(parser)


if __name__ == '__main__':
    main()
