#!/usr/bin/env python
# -*- coding: utf-8

import sys
import numpy as numpy

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db", "profile-db", "collection", "bin"]
__provides__ = ["view-data", "misc-data-items-txt"]
__resources__ = [("This program in action as part of the metapangenomic workflow", "http://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes")]
__description__ = ("Quantify the detection of genes in genomes in metagenomes to identify the "
                   "environmental core. This is a helper script for anvi'o metapangenomic workflow")


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

    fraction_of_median_coverage = float(args.fraction_of_median_coverage)
    min_detection = float(args.min_detection)

    run.info('Fraction of median coverage for core genes', fraction_of_median_coverage)
    run.info('Min detection in at last one metagenome', min_detection)

    progress.new('Initializing internal genomes')

    genome_name = args.bin_id

    if not genome_name:
        raise ConfigError("You forgot to provide a '--bin' parameter for the bin you wish to focus on :/")

    progress.update('Computing gene presence data ...')
    gene_presence_in_the_environment_dict = {}

    summary = summarizer.ProfileSummarizer(args)
    summary.init()
    summary.init_gene_level_coverage_stats_dicts()

    # recover the detection statistic
    summary.init_collection_profile(args.collection_name)
    detection_across_metagenomes = summary.collection_profile[genome_name]['detection']

    num_metagenomes_above_min_detection = [m for m in detection_across_metagenomes if detection_across_metagenomes[m] > min_detection]

    num_metagenomes = len(detection_across_metagenomes)
    average_detection_of_genome = sum(detection_across_metagenomes.values()) / num_metagenomes
    max_detection_of_genome = max(detection_across_metagenomes.values())
    if not len(num_metagenomes_above_min_detection):
        not_enough_detection = True
        run.warning("This genome is not detected in any of the %d metgenomes you have more than %.2f of its fraction. In fact, its "
                    "average detection is %.2f across all metagenomes, %.2f being the maximum. Well, you can take a look at the "
                    "anvi-summarize output if you would like to learn more. The bottom line is this: the metagenomes you have in this "
                    "project are not the best ones to investigate the occurrence of the genes in this genome in the environment. So. "
                    "This program will give you quite useless output files, but anvi'o hopes that you will enjoy them anyway :(" %\
                            (num_metagenomes, min_detection, average_detection_of_genome, max_detection_of_genome))
    else:
        not_enough_detection = False

    gene_presence_in_the_environment_dict = {}

    genome_bin_summary = summarizer.Bin(summary, genome_name)

    gene_coverages_across_samples = {}
    for gene_callers_id in genome_bin_summary.gene_level_coverage_stats_dict:
        gene_coverages_across_samples[gene_callers_id] = {}
        for sample_name in genome_bin_summary.gene_level_coverage_stats_dict[gene_callers_id]:
            gene_coverages_across_samples[gene_callers_id][sample_name] = genome_bin_summary.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['mean_coverage']

    # at this point we have all the genes in the genome bin. what we need is to characterize their detection. first,
    # summarize the coverage of each gene in all samples:
    sum_gene_coverages_across_samples = dict([(gene_callers_id, sum(gene_coverages_across_samples[gene_callers_id].values())) for gene_callers_id in gene_coverages_across_samples])

    # now we will identify the median coverage
    median_coverage_across_samples = numpy.median(list(sum_gene_coverages_across_samples.values()))

    # now we will store decide whether a gene found in this genome is also found in the environment, and store that
    # information into `gene_presence_in_the_environment_dict`, and move on to the next stage.
    for gene_caller_id in sum_gene_coverages_across_samples:
        if not_enough_detection:
            gene_presence_in_the_environment_dict[gene_caller_id] = 'UNKNOWN'
        elif sum_gene_coverages_across_samples[gene_caller_id] < median_coverage_across_samples * fraction_of_median_coverage:
            gene_presence_in_the_environment_dict[gene_caller_id] = 'NOT DETECTED'
        else:
            gene_presence_in_the_environment_dict[gene_caller_id] = 'DETECTED'

    progress.update('Storing the output for gene presence data')

    utils.store_dict_as_TAB_delimited_file(gene_coverages_across_samples, '%s-GENE-COVs.txt' % genome_name)

    output = open(genome_name + '-ENV-DETECTION.txt', 'w')
    output.write('gene_callers_id\tdetection\n')
    for gene_callers_id in gene_presence_in_the_environment_dict:
        output.write('%s\t%s\n' % (gene_callers_id, gene_presence_in_the_environment_dict[gene_callers_id]))
    output.close()

    progress.end()


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    parser.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    parser.add_argument('--min-detection', metavar="FLOAT", default=0.50, type=float, help="For this entire thing to work, the\
                        genome you are focusing on should be detected in at least one metagenome. If that is not the case, it would\
                        mean that you do not have any sample that represents the niche for this organism (or you do not have enough\
                        depth of coverage) to investigate the detection of genes in the environment. By default, this script requires\
                        at least '0.5' of the genome to be detected in at least one metagenome. This parameter allows you to change\
                        that. 0 would mean no detection test required, 1 would mean the entire genome must be detected.")
    parser.add_argument('--fraction-of-median-coverage', metavar="FLOAT", default=0.25, type=float, help="The value set here\
                        will be used to remove a gene if its total coverage across environments is less than the median coverage\
                        of all genes multiplied by this value. The default is 0.25, which means, if the median total coverage of\
                        all genes across all samples is 100X, then, a gene with a total coverage of less than 25X across all\
                        samples will be assumed not a part of the 'environmental core'.")


    return parser.get_args(parser)


if __name__ == '__main__':
    main()
