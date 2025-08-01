#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.genomestorage as genomestorage

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import TableForItemAdditionalData

__citation__ = "'mOTUlizer' by Buck et al (https://doi.org/10.1101/2021.06.25.449606)"

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['moritzbuck']
__requires__ = ["pan-db", "genomes-storage-db"]
__provides__ = ["bin"]
__resources__ = [("GitHub repository for the mOTUPan", "https://github.com/moritzbuck/mOTUlizer/")]
__description__ = "Runs mOTUpan on your gene clusters to estimate whether they are core or accessory"


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

    try:
        from mOTUlizer.classes.mOTU import mOTU
    except ModuleNotFoundError:
        raise ConfigError(f"{__citation__} does not seem to be installed on your system as anvi'o was "
                          f"unable to import its Python libraries :( In theory you should be able to install "
                          f"mOTUlizer by running the following command in your anvi'o environment: pip install "
                          f"mOTUlizer. Please visit the URL https://github.com/moritzbuck/mOTUlizer to read most "
                          f"up-to-date installation instructions or to report issues. Thank you for your patience!")


    if not args.output_file and not args.store_in_db:
        if args.just_do_it:
            run.warning("Lol. You are making anvi'o compute stuff that will not be stored anywhere. Weird you :(", lc="green")
        else:
            raise ConfigError("But why no select a way to report this stuff? :/ You can ask anvi'o to store your results as a "
                              "a TAB-delimited file. Or add them to directly to the database. If you insist that you just want "
                              "to run this analysis without really storing anything anywhere you can always use `--just-do-it` "
                              "(because why not, it is your computer).")

    if args.output_file:
        filesnpaths.is_output_file_writable(args.output_file)

    pan = dbops.PanSuperclass(args, r=terminal.Run(verbose=False))

    # get a pan instance
    pan = dbops.PanSuperclass(args)

    # complain if there is no genomes storage
    if not pan.genomes_storage_is_available:
        raise ConfigError("The anvi'o pan class does not see a genomes storage. No genomes storage no cake.")

    # get a genome storage instance
    genome_storage = genomestorage.GenomeStorage(pan.genomes_storage_path, run=terminal.Run(verbose=False))

    # recover genome completeness values
    genome_completeness_dict = {g: v.get('percent_completion', None) for g, v in genome_storage.genomes_info.items()}
    genome_lengths_dict = {g: v.get('total_length', None) for g, v in genome_storage.genomes_info.items()}

    # a heuristic if there are genomes missing completeness estimates
    genomes_with_missing_completeness_data = [g for g in genome_completeness_dict if genome_completeness_dict[g] in [None, -1]]
    if genomes_with_missing_completeness_data:
        run.warning(f"Please read this carefully as {len(genomes_with_missing_completeness_data)} of your {len(genome_completeness_dict)} "
                    f"of your genomes are missing completion estimates. This may be due to multiple reasons: you may have not run the anvi'o "
                    f"program `anvi-run-hmms` on your contigs databases before generating the genome storage database. Or, some of your "
                    f"genomes may be too incomplete to have proper copmletion estimates. To continue this workflow, this prgram will something "
                    f"a bit tacky: it will find the genomes with highest completion value, and assign a completion value for the genomes in "
                    f"in your collection that are lacking completion estimates based on how their lenghts compare to the lenght of the genome "
                    f"with a known completion estimate. This indeed is a pretty brutal approximation to an optimal solution, but this is quite "
                    f"an unoptimal situation at the first place. Here are the list of genomes that makes you read this warning: "
                    f"{', '.join(genomes_with_missing_completeness_data)}.", header="PLEASE NOTE MISSING COMPLETION VALUES")

        highest_completion_value = sorted(genome_completeness_dict.items(), key=lambda x:x[1], reverse=True)[0][1]

        if not highest_completion_value or highest_completion_value < 50:
            raise ConfigError('Bad news :( Your most highly complete genome is less than 50% complete. This is futile.')

        # here we turn this value into an integer, so we can recruit more genomes. I.e., if a genome is 99.5 and another
        # is 99.1, it is best to take the average of their lenghts rather than just using the most complete one since
        # all these thigns are so arbitrary
        highest_completion_value = int(highest_completion_value)

        # subest the highly complete genomes nad learn their average lenghts:
        highly_complete_genomes = [g for g in genome_completeness_dict if genome_completeness_dict[g] > highest_completion_value]
        highly_complete_genome_lengths = [genome_lengths_dict[g] for g in highly_complete_genomes]
        avg_length_for_highly_complete_genomes = sum(highly_complete_genome_lengths) / len(highly_complete_genome_lengths)

        # update missing completion estimates
        for genome_name in genomes_with_missing_completeness_data:
            new_completion_value = genome_lengths_dict[genome_name] / avg_length_for_highly_complete_genomes * highest_completion_value
            genome_completeness_dict[genome_name] = 100 if new_completion_value > 100 else new_completion_value

        # let the user know what just happened
        run.warning(f"Missing completion values for {len(genomes_with_missing_completeness_data)} genomes are now "
                    f"recovered based on the lengths and completion values of these genomes: "
                    f"{', '.join(highly_complete_genomes)}", header="MISSING COMPLETION VALUES RECOVERED", lc="green")

    # initialize gene clusters
    pan.init_gene_clusters()

    # build the gene clusters dict for mOTU
    gene_clusters_dict = {g : set() for g in pan.genome_names}
    for gc, hits in pan.gene_clusters.items():
        for genome, genes in hits.items():
            if len(genes) > 0:
                gene_clusters_dict[genome].add(gc)

    # bleep bloop
    run.warning(f"Anvi'o is now using {__citation__} to identify gene clusters in your pangenome based on whether they are "
                f"core or accessory through a Bayesian calculation. If you have used the `--store-in-db` flag, this program "
                f"will add two new layers to your pan database: 'Gene_cluster_type' (i.e., whether a given gene cluster is "
                f"accessor or core) and 'Gene_cluster_type_LLR' (the log likelihood ratio of that assessment). When you "
                f"publish your findings, please do not forget to properly credit this work.", lc='green', header="CITATION")

    motu = mOTU(name="mOTUpan_core_prediction",
                faas={},
                gene_clusters_dict=gene_clusters_dict,
                genome_completion_dict=genome_completeness_dict,
                max_it=100,
                threads=args.num_threads,
                precluster=False,
                method='default')

    if args.num_bootstraps :
        motu.roc_values(boots=args.num_bootstraps)
        rocs = motu.roc_values(0)
        run.warning("You ran the bootstrap option for the pangenome predition, if you have no output file the results are only shown here :", lc = "gray")
        for k,v in rocs.items() :
            run.info(k,"{:.2f}".format(v) if k != 'nb_bootstraps' else v, lc = "gray")

    if args.store_in_db:
        data2add = {g: {'Gene_cluster_type' : 'CORE' if g in motu.core else 'ACCESSORY', 'Gene_cluster_type_LLR' : motu.likelies[g]} for g in pan.gene_cluster_names}
        panditional_data = TableForItemAdditionalData(args)
        panditional_data.add(data2add, data_keys_list=['Gene_cluster_type', 'Gene_cluster_type_LLR'])

    if args.output_file:
        with open(args.output_file, "w") as handle:
            print(motu.pretty_pan_table(), file = handle)
        run.info("Output file", args.output_file, mc="green", nl_before=1)


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

    groupC = parser.add_argument_group('OPTIONAL',"Optional stuff available for you to use")
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupC.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupC.add_argument('-B', '--num-bootstraps' , metavar = "NUM_BOOTSTRAPS", default = 0, type = int, help = "number of boostraps run on the partitioning to evaluate it's quality"  )

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
