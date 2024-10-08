#!/usr/bin/env python
# -*- coding: utf-8
"""Code for genome similarity calculation"""

import os
import shutil
import argparse
import pandas as pd

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.genomedescriptions as genomedescriptions

from itertools import combinations

from anvio.errors import ConfigError
from anvio.drivers import pyani, sourmash, fastani
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.tables.miscdata import TableForLayerOrders

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Mahmoud Yousef"
__email__ = "mahmoudyousef@uchicago.edu"

run = terminal.Run()
progress = terminal.Progress()

J = lambda *args: os.path.join(*args)


class Dereplicate:
    def __init__(self, args):
        self.args = args

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x

        # input
        self.internal_genomes = A('internal_genomes', null)
        self.external_genomes = A('external_genomes', null)
        self.fasta_text_file = A('fasta_text_file', null)
        self.ani_dir = A('ani_dir', null)
        self.mash_dir = A('mash_dir', null)
        # mode
        self.program_name = A('program', null)
        self.representative_method = A('representative_method', null)
        self.similarity_threshold = A('similarity_threshold', null)
        # fastANI specific
        self.fastani_kmer_size = A('fastani_kmer_size', null)
        self.fragment_length = A('fragment_length', null)
        self.min_num_fragments = A('min_num_fragments', null)
        # pyANI specific
        self.min_alignment_fraction = A('min_alignment_fraction', null)
        self.significant_alignment_length = A('significant_alignment_length', null)
        self.min_percent_identity = A('min_percent_identity', null)
        self.use_full_percent_identity = A('use_full_percent_identity', null)
        # sourmash specific
        self.kmer_size = A('kmer_size', null)
        self.scale = A('scale', null)
        # output
        self.output_dir = A('output_dir', null)
        self.report_all = A('report_all', null)
        self.skip_fasta_report = A('skip_fasta_report', null)
        self.just_do_it = A('just_do_it', null)

        self.import_previous_results = False
        self.sequence_source_provided = False

        self.sanity_check()
        self.program_info = self.get_program_specific_info()

        self.clusters = {}
        self.cluster_report = {}
        self.genomes_info_dict = {}
        self.cluster_to_representative = {}
        self.genome_name_to_cluster_name = {}


    def get_program_specific_info(self):
        if self.program_name == 'pyANI':
            metric_name = 'percentage_identity' if not self.use_full_percent_identity else 'full_percentage_identity'
            necessary_reports = [metric_name, 'alignment_coverage']
        elif self.program_name == 'sourmash':
            metric_name = 'mash_similarity'
            necessary_reports = [metric_name]
        elif self.program_name == 'fastANI':
            metric_name = 'ani'
            necessary_reports = [metric_name]

        return {
            'metric_name': metric_name,
            'necessary_reports': necessary_reports
        }


    def is_genome_names_compatible_with_similarity_matrix(self, similarity_matrix, genome_names):
        missing_in_matrix = [n for n in genome_names if n not in similarity_matrix]
        missing_in_names = [n for n in similarity_matrix if n not in genome_names]

        should_raise = False
        if len(missing_in_matrix):
            should_raise = True
            raise_info = ('sequence source names', 'similarity matrix',
                          missing_in_matrix[:5] if len(missing_in_matrix) < 6 else missing_in_matrix)
        elif len(missing_in_names):
            should_raise = True
            raise_info = ('similarity matrix names', 'sequence sources',
                          missing_in_names[:5] if len(missing_in_names) < 6 else missing_in_names)

        if should_raise:
            raise ConfigError("At least one of your %s does not appear in the %s. This could be due "
                              "to many different reasons. Probably, the imported similarity results were generated from sequence sources "
                              "that are somehow different from the ones you provided here. If you imported your results (e.g "
                              "with --ani-dir, --mash-dir, etc), we recommend re-running with standard inputs. "
                              "See `INPUT OPTIONS` within the help menu. Here are some that are missing: %s." % (*raise_info, ))


    def sanity_check(self):
        if any([self.fasta_text_file, self.external_genomes, self.internal_genomes]):
            self.sequence_source_provided = True

        if not self.sequence_source_provided and self.report_all:
            raise ConfigError("You didn't provide a sequence source, but want to `--report-all sequences`. Hmmm.")

        if self.report_all and self.skip_fasta_sequences:
            raise ConfigError("You want to report all sequences (--report-all), but you also want to skip "
                              "reporting all sequences (--skip-fasta-report). Hmmm.")

        if not self.program_name:
            additional_text = ", even when you are importing results from a previous run"
            raise ConfigError("You must provide a program name to dereplicate your genomes using the `--program` "
                              "parameter%s. Chocies include %s." % (additional_text if self.ani_dir else '', ', '.join(list(program_class_dictionary.keys()))))

        if not any([self.program_name, self.ani_dir, self.mash_dir]):
            raise ConfigError("Anvi'o could not determine how you want to dereplicate "
                             "your genomes. Please take a close look at your parameters: either --program needs to be "
                             "set, or an importable directory (e.g. --ani-dir, --mash-dir, etc) needs to be provided.")

        if self.program_name not in list(program_class_dictionary.keys()):
            raise ConfigError("Anvi'o is impressed by your dedication to dereplicate your genomes through %s, but "
                             "%s is not compatible with `anvi-dereplicate-genomes`. Anvi'o can only work with pyANI "
                             "and sourmash separately." % (self.program_name, self.program_name))

        if self.ani_dir and self.mash_dir:
            raise ConfigError("Anvi'o cannot currently dereplicate using both ANI and mash similarity results at "
                             "the same time. Please pick one to use for dereplication")

        if self.ani_dir or self.mash_dir:
            self.import_previous_results = True

            if self.sequence_source_provided:
                additional_msg = ''

                if not self.just_do_it:
                    raise ConfigError("You provided any of {--external-genomes, --internal-genomes, --fasta-text-file} "
                                      "*alongside* a results directory. This requires that the external, internal, "
                                      "and fasta inputs are exactly those used to generate the results directory. "
                                      "This is sketchy and anvi'o doesn't recommend it. If you insist, re-run with "
                                      "--just-do-it. Otherwise, have the results regenerated here by removing '%s' "
                                      "as an input." % (self.ani_dir or self.mash_dir))
            else:
                additional_msg = ("In addition, no FASTAs will be generated since you did not provide any sequence "
                                  "sources for anvi'o.")

            run.warning("You chose to work with an already existing results folder. Please keep in mind that you "
                        "are now burdened with the responsibility of knowing what parameters you used to generate "
                        "these results.%s" % additional_msg)

        if self.ani_dir and not self.program_name in ['pyANI', 'fastANI']:
            raise ConfigError("You provided a pre-existing directory of ANI results (--ani-dir), but also provided a program "
                              "name ('%s') that was not compatible with ANI." % self.program_name)

        if self.mash_dir and not self.program_name in ['sourmash']:
            raise ConfigError("You provided a pre-existing directory of mash results (--mash-dir), but also provided a program "
                              "name ('%s') that was not compatible with mash." % self.program_name)

        if self.min_alignment_fraction < 0 or self.min_alignment_fraction > 1:
            if self.program_name == "pyANI":
                raise ConfigError("Alignment coverage is a value between 0 and 1. Your cutoff alignment coverage "
                                 "value of %.2f doesn't fit in these boundaries" % self.min_alignment_fraction)

        if self.similarity_threshold < 0 or self.similarity_threshold > 1:
            raise ConfigError("When anvi'o collapses %s's output into a similarity matrix, all values are reported as "
                             "similaritys between 0 and 1. %.2f can't be used to determine redundant genomes"\
                              % (self.program_name, self.similarity_threshold))

        if self.representative_method == "Qscore" and not (self.external_genomes or self.internal_genomes):
                self.representative_method = "centrality"
                run.warning("Anvi'o must be provided either an external or internal genome collection (or both) to be "
                            "used with 'Qscore', since this is the only way for anvi'o to learn about completion and "
                            "redundancy scores. Anvi'o will switch to 'centrality'")

        if self.representative_method == "length" and not self.sequence_source_provided:
                self.representative_method = "centrality"
                run.warning("Anvi'o must be provided a sequence source (external genomes, internal genomes, fasta "
                            "text file, or a combination thereof) to be used with 'length', since this is the only "
                            "way for anvi'o to learn about genome lengths. Anvi'o will switch to 'centrality'")


    def init_genome_similarity(self):
        args = self.args
        args.output_dir_exists = True

        self.similarity = program_class_dictionary[self.program_name](args)


    def get_genome_names(self):
        """
        genome_names are learned from the GenomeSimilarity class in self.similarity.__init__. But if similarity scores are
        imported and sequence sources were not provided, GenomeSimilarity obviously knows nothing of the genome_names. So
        in this fringe case we grab genome names from the results matrix
        """
        return (self.similarity.genome_names
                if self.similarity.genome_names
                else set(self.similarity_matrix.keys()))


    def get_similarity_matrix(self):
        return (self.import_similarity_matrix()
                if self.import_previous_results
                else self.gen_similarity_matrix())


    def gen_similarity_matrix(self):
        self.similarity.process(self.temp_dir)

        try:
            similarity_matrix = self.similarity.results[self.program_info['metric_name']]
        except KeyError:
            # With sourmash you don't always know the metric name, you only be sure of what it
            # contains. This is because the kmer is a part of the result name. This is my fault but
            # I'm too lazy to fix the design because sourmash is not appropriate for genome
            # comparison anyways.
            for result_name in self.similarity.results:
                if self.program_info['metric_name'] in result_name:
                    similarity_matrix = self.similarity.results[result_name]
                    break

        run.info('%s similarity metric' % self.program_name, 'calculated')

        if anvio.DEBUG:
            import json
            for report in self.similarity.results:
                run.warning(None, header=report)
                print(json.dumps(self.similarity.results[report], indent=2))

        return similarity_matrix


    def import_similarity_matrix(self):
        dir_name, dir_path = ('--ani-dir', self.ani_dir) if self.program_name in ['pyANI', 'fastANI'] else ('--mash-dir', self.mash_dir)

        if filesnpaths.is_dir_empty(dir_path):
            raise ConfigError("The %s you provided is empty. What kind of game are you playing?" % dir_name)
        files_in_dir = os.listdir(dir_path)

        for report in self.program_info['necessary_reports']:
            report_name = report + ".txt"
            matching_filepaths = [f for f in files_in_dir if report_name in f]

            if self.program_info['metric_name'] == 'percentage_identity':
                # FIXME very very bad block of code here. Why should this method know about
                # percentage_identity or full_percentage_identity?
                matching_filepaths = [f for f in matching_filepaths if 'full_percentage_identity' not in f]

            if len(matching_filepaths) > 1:
                raise ConfigError("Your results directory contains multiple text files for the matrix %s. "
                                 "Please clean up your directory and make sure that only one text file of this "
                                 "report exists" % report)
            elif not len(matching_filepaths):
                raise ConfigError("Your results directory does not have a text file for the report %s. "
                                 "Anvi'o cannot dereplicate genomes from prevous results without this report" % report)

            self.similarity.results[report] = utils.get_TAB_delimited_file_as_dictionary(J(dir_path, matching_filepaths[0]))

        run.info('%s results directory imported from' % self.program_name, dir_path)

        return self.similarity.results[self.program_info['metric_name']]


    def clean(self):
        if self.temp_dir:
            if anvio.DEBUG:
                run.warning("The temp directory, %s, is kept. Please don't forget to clean it up "
                            "later" % self.temp_dir, header="Debug")
            else:
                run.info_single('Cleaning up the temp directory (you can use `--debug` if you would '
                                'like to keep it for testing purposes)', nl_before=1, nl_after=1)

                shutil.rmtree(self.temp_dir)
                self.temp_dir = None


    def report(self):
        if self.sequence_source_provided and not self.skip_fasta_report:
            self.populate_genomes_dir()
            utils.store_dataframe_as_TAB_delimited_file(self.gen_fasta_report_output(), self.FASTA_REPORT_path)

        if not self.import_previous_results:
            self.similarity.report(self.SIMILARITY_SCORES_dir, ok_if_exists=True)

        utils.store_dataframe_as_TAB_delimited_file(self.gen_cluster_report_output(), self.CLUSTER_REPORT_path)


    def populate_genomes_dir(self):
        if self.report_all:
            temp_paths = {name: path for name, path in self.similarity.name_to_temp_path.items()}
        else:
            temp_paths = {name: path for name, path in self.similarity.name_to_temp_path.items() if name in self.cluster_to_representative.values()}

        # populate a path dictionary for gen_fasta_report_output
        self.output_fasta_paths = {}

        for name, temp_path in temp_paths.items():
            output_path = J(self.GENOMES_dir, name + '.fa')
            shutil.copy(src = temp_path, dst = output_path)
            self.output_fasta_paths[name] = output_path.lstrip(self.output_dir + '/')


    def init_output_dir(self):
        """
        DIRECTORY STRUCTURE:
        ===================
            output_dir/
            ├── GENOMES/
            ├── SIMILARITY_SCORES/
            ├── FASTA_REPORT.txt
            └── CLUSTER_REPORT.txt

        GENOMES/
            A folder with genomes. Each genome is a fasta file. Not provided if no sequence sources
            are provided, or if --skip-fasta-report.
        SIMILARITY_SCORES/
            A folder containing the output of similarity scores. Not provided if previous similarity
            scores are imported.
        FASTA_REPORT.txt
            A text file detailing the fasta paths of the output. Not provided if no sequence sources
            are provided, or if --skip-fasta-report.
        CLUSTER_REPORT.txt
            A text file detailing which genomes were determined to be redundant with one another.
        """
        filesnpaths.check_output_directory(self.output_dir, ok_if_exists=False)
        os.mkdir(self.output_dir)

        self.CLUSTER_REPORT_path = J(self.output_dir, 'CLUSTER_REPORT.txt')

        if self.sequence_source_provided and not self.skip_fasta_report:
            self.GENOMES_dir = J(self.output_dir, 'GENOMES')
            self.FASTA_REPORT_path = J(self.output_dir, 'FASTA_REPORT.txt')
            os.mkdir(self.GENOMES_dir)
        else:
            self.GENOMES_dir = None

        if not self.import_previous_results:
            self.SIMILARITY_SCORES_dir = J(self.output_dir, 'SIMILARITY_SCORES')
            os.mkdir(self.SIMILARITY_SCORES_dir)
        else:
            self.SIMILARITY_SCORES_dir = None


    def process(self):
        run.info('Run mode', self.program_name)

        self.init_output_dir()

        # inits self.similarity
        self.init_genome_similarity()

        # will hold a directory of fasta files to calculate similarity matrix
        self.temp_dir = self.similarity.get_fasta_sequences_dir() if self.sequence_source_provided else None

        self.similarity_matrix = self.get_similarity_matrix()
        self.genome_names = self.get_genome_names()

        self.is_genome_names_compatible_with_similarity_matrix(self.similarity_matrix, self.genome_names)

        run.info('Number of genomes considered', len(self.genome_names))

        self.init_clusters()
        self.populate_genomes_info_dict()
        self.dereplicate()


    def init_clusters(self):
        """Each genome starts in its own cluster, i.e. as many clusters as genomes"""
        for count, genome_name in enumerate(self.genome_names, start=1):
            cluster_name = 'cluster_%06d' % count
            self.genome_name_to_cluster_name[genome_name] = cluster_name
            self.clusters[cluster_name] = set([genome_name])


    def populate_genomes_info_dict(self):
        full_dict = {}

        if self.similarity.genome_desc:
            full_dict = self.similarity.genome_desc.genomes

        if self.similarity.fasta_txt:
            fastas = utils.get_TAB_delimited_file_as_dictionary(self.similarity.fasta_txt, expected_fields=['name', 'path'], only_expected_fields=True)

            for name in fastas:
                full_dict[name] = {}
                full_dict[name]['percent_completion'] = None
                full_dict[name]['percent_redundancy'] = None
                full_dict[name]['total_length'] = sum(utils.get_read_lengths_from_fasta(fastas[name]['path']).values())

        if self.representative_method == 'Qscore':
            missing_completion = False
            for genome in full_dict:
                if 'percent_completion' not in full_dict[genome] or 'percent_redundancy' not in full_dict[genome]:
                    self.representative_method = 'centrality'
                    missing_completion = genome

            if missing_completion:
                run.warning("At least one of your genomes does not have completion and/or redundancy scores, which makes "
                            "it impossible to use Qscore to pick best representatives from each cluster. One of these "
                            "genomes is %s. Anvi'o switched you to the 'centrality' method for picking representatives." % (genome))

        self.genomes_info_dict = full_dict


    def are_redundant(self, genome1, genome2):
        cluster1 = self.genome_name_to_cluster_name[genome1]
        cluster2 = self.genome_name_to_cluster_name[genome2]
        return True if cluster1 == cluster2 else False


    def gen_cluster_report(self):
        if not any([self.clusters, self.cluster_to_representative]):
            raise ConfigError("gen_cluster_report :: Run dereplicate() before trying to generate a cluster report")

        self.cluster_report = {
            cluster_name: {
                'representative': self.cluster_to_representative[cluster_name],
                'genomes': self.clusters[cluster_name],
                'size': len(self.clusters[cluster_name]),
            } for cluster_name in self.clusters
        }


    def gen_cluster_report_output(self):
        cluster_report_output = pd.DataFrame(self.cluster_report).T.reset_index().rename(columns={'index': 'cluster'})
        cluster_report_output['genomes'] = cluster_report_output['genomes'].apply(lambda s: ','.join(s))
        return cluster_report_output[['cluster', 'size', 'representative', 'genomes']]


    def gen_fasta_report_output(self):
        fasta_report_output = {
            'name': [],
            'source': [],
            'cluster': [],
            'cluster_rep': [],
            'path': [],
        }

        for name, path in self.output_fasta_paths.items():
            cluster = self.genome_name_to_cluster_name[name]

            fasta_report_output['name'].append(name)
            fasta_report_output['source'].append(self.similarity.genome_sources[name])
            fasta_report_output['cluster'].append(cluster)
            fasta_report_output['cluster_rep'].append(self.cluster_report[cluster]['representative'])
            fasta_report_output['path'].append(path)

        return pd.DataFrame(fasta_report_output).sort_values(by=['cluster', 'name'])


    def rename_clusters(self):
        conversion = {key: 'cluster_%06d' % count for count, key in enumerate(self.clusters.keys(), start=1)}
        self.clusters = {conversion[key]: self.clusters[key] for key in self.clusters}
        self.genome_name_to_cluster_name = {k: conversion[v] for k, v in self.genome_name_to_cluster_name.items()}


    def dereplicate(self):
        num_genome_pairs = sum(1 for _ in combinations(self.genome_names, 2))
        progress.new('Dereplication', progress_total_items=num_genome_pairs)

        counter = 0
        genome_pairs = combinations(self.genome_names, 2)
        for genome1, genome2 in genome_pairs:
            if counter % 10000 == 0:
                progress.increment(increment_to=counter)
                progress.update('%d / %d pairwise comparisons made' % (counter, num_genome_pairs))

            if genome1 == genome2 or self.are_redundant(genome1, genome2):
                counter += 1
                continue

            similarity = float(self.similarity_matrix[genome1][genome2])
            if similarity >= self.similarity_threshold:
                self.update_clusters(genome1, genome2)

            counter += 1

        progress.increment(increment_to=num_genome_pairs)
        progress.update('All %d pairwise comparisons have been made' % num_genome_pairs)

        # remove empty clusters and rename so that names are sequential
        self.clusters = {cluster: genomes for cluster, genomes in self.clusters.items() if genomes}
        self.rename_clusters()

        self.cluster_to_representative = self.get_representative_for_each_cluster()
        self.gen_cluster_report()

        progress.end()

        run.info('Number of redundant genomes', len(self.genome_names) - len(self.cluster_report))
        run.info('Final number of dereplicated genomes', len(self.cluster_report))


    def update_clusters(self, genome1, genome2):
        cluster1 = self.genome_name_to_cluster_name[genome1]
        cluster2 = self.genome_name_to_cluster_name[genome2]

        if cluster1 == cluster2:
            return

        # empty the smaller cluster
        from_cluster, to_cluster = (cluster1, cluster2) \
                                   if len(self.clusters[cluster1]) < len(self.clusters[cluster2]) \
                                   else (cluster2, cluster1)

        self.genome_name_to_cluster_name.update({g: to_cluster for g in self.clusters[from_cluster]})
        self.clusters[to_cluster]   |= self.clusters[from_cluster]
        self.clusters[from_cluster] -= self.clusters[from_cluster]


    def pick_representative_with_largest_Qscore(self, cluster):
        """This function will return the a genome in a cluster with highest substantive completion esitmate.

        If there are multiple genomes with the same substantive completion estimate, then it will return the
        longest one. If there are multiple that has the same substantive completion estimate AND length, then
        it will return the first one.
        """

        if not cluster:
            return None

        if len(cluster) == 1:
            return list(cluster)[0]

        # get all substantive completion and lenght values for genomes within the cluster
        substantive_completion_and_length_values = [(g, self.genomes_info_dict[g]['percent_completion'] - self.genomes_info_dict[g]['percent_redundancy'], self.genomes_info_dict[g]['total_length']) for g in cluster]

        # calculate the maximum substantive completion value found
        max_substantive_completion = max([e[1] for e in substantive_completion_and_length_values])
        
        # get all the genomes with that exact completion value, sorted by their length (in case there are more than one genomes with identical completion estimates)
        genomes_with_max_substantive_completion = sorted([e for e in substantive_completion_and_length_values if e[1] == max_substantive_completion], key=lambda x: x[2], reverse=True)

        # return the first one (which is supposed to be the longest one)
        return genomes_with_max_substantive_completion[0][0]


    def pick_representative_with_largest_genome(self, cluster):
        lengths_dict = {name: self.genomes_info_dict[name]['total_length'] for name in cluster}
        return max(lengths_dict, key = lambda name: lengths_dict.get(name, 0))


    def pick_representative_with_largest_centrality(self, cluster):
        d = {}
        for name in cluster:
            d[name] = 0

            for target in cluster:
                d[name] += float(self.similarity_matrix[name][target])

        new_dict = {}
        for name, val in d.items():
            new_dict[val] = name

        max_val = max(new_dict.keys())

        return new_dict[max_val]


    def get_representative_for_each_cluster(self):
        cluster_to_representative = {}
        for cluster_name in self.clusters.keys():
            cluster = self.clusters[cluster_name]

            if cluster == []:
                continue

            if self.representative_method == 'Qscore':
                representative_name = self.pick_representative_with_largest_Qscore(cluster)
            elif self.representative_method == 'length':
                representative_name = self.pick_representative_with_largest_genome(cluster)
            elif self.representative_method == 'centrality':
                representative_name = self.pick_representative_with_largest_centrality(cluster)

            cluster_to_representative[cluster_name] = representative_name

        return cluster_to_representative



class GenomeSimilarity:
    def __init__(self, args):
        self.args = args

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.fasta_txt = A('fasta_text_file', null)
        self.internal_genomes = A('internal_genomes', null)
        self.external_genomes = A('external_genomes', null)
        self.clustering_distance = A('centrality', null) or constants.distance_metric_default
        self.clustering_linkage = A('linkage', null) or constants.linkage_method_default
        self.output_dir = A('output_dir', null)
        self.pan_db = A('pan_db', null)
        self.output_dir_exists = A('output_dir_exists', null)
        self.just_do_it = A('just_do_it', null)
        self.representative_method = A('representative_method', null)

        if (self.internal_genomes or self.external_genomes):
            self.genome_desc = genomedescriptions.GenomeDescriptions(args, run = terminal.Run(verbose=False))
        else:
            self.genome_desc = None

        self.clusterings = {}
        self.hash_to_name = {}
        self.name_to_temp_path = {}

        self.sanity_check(ok_if_exists=self.output_dir_exists)
        self.genome_names, self.genome_sources = self.get_genome_names_list_and_genome_source_dict()


    def sanity_check(self, ok_if_exists):
        clustering.is_distance_and_linkage_compatible(self.clustering_distance, self.clustering_linkage)
        filesnpaths.check_output_directory(self.output_dir, ok_if_exists=(ok_if_exists or self.just_do_it))

        if self.pan_db:
            utils.is_pan_db(self.pan_db)


    def cluster(self):
        for report_name in self.results:
            try:
                self.clusterings[report_name] = clustering.get_newick_tree_data_for_dict(self.results[report_name],
                                                                                         linkage=self.clustering_linkage,
                                                                                         distance=self.clustering_distance)
            except Exception as e:
                if anvio.DEBUG:
                    print(e)
                run.warning("Bad news :/ Something went wrong while anvi'o was processing the clustering for "
                            "'%s'. You can find the offending file if you search for the output file in "
                            "the temporary output directory '%s'. You potentially waited so long, so anvi'o "
                            "will continue as if nothing happened. If you want to diagnose this problem, rerun with "
                            "--debug" % (report_name, J(self.temp_dir, 'output')))


    def add_to_pan_db(self):
        if self.pan_db:
            utils.is_pan_db(self.pan_db)
            pan = db.DB(self.pan_db, anvio.__pan__version__)

            db_genome_names = set([])
            G = lambda g: pan.get_meta_value(g).strip().split(',')
            for genome_name in G('external_genome_names') + G('internal_genome_names'):
                db_genome_names.add(genome_name) if len(genome_name) else None

            found_only_in_db = db_genome_names.difference(self.genome_names)
            found_only_in_results = self.genome_names.difference(db_genome_names)

            if len(found_only_in_results) > 0:
                raise ConfigError("Some genome names found in your results do not seem to be exist in the pan database. "
                   "Here are the list of them: " + ", ".join(list(found_only_in_results)))

            if len(found_only_in_db) > 0:
                run.warning("Some genome names found in pan database do not seem to be exist in the similarity report. "
                   "anvi'o will add the ones that are found in the database anyway, but here is the list of missing ones: "
                   "" + ", ".join(list(found_only_in_db)))

            run.warning(None, header="MISC DATA MAGIC FOR YOUR PAN DB", lc='green')
            for report_name in self.results:
                target_data_group = '%s_%s' % (self.similarity_type, report_name)
                l_args = argparse.Namespace(pan_db=self.pan_db, just_do_it=True, target_data_group=target_data_group)
                TableForLayerAdditionalData(l_args, r=terminal.Run(verbose=False)).add(self.results[report_name], list(self.results[report_name].keys()))

                clustering_present = True if report_name in self.clusterings else False
                if clustering_present:
                    TableForLayerOrders(self.args, r=terminal.Run(verbose=False)).add({self.method + '_' + report_name: {'data_type': 'newick',
                                                                                       'data_value': self.clusterings[report_name]}})

                extra_msg = ' and order' if clustering_present else ''
                run.info_single("Additional data%s for %s are now in pan db" % (extra_msg, target_data_group.replace('_', ' ')), mc='green')


    def report(self, output_dir = None, ok_if_exists = False):
        output_dir = output_dir if output_dir else self.output_dir

        filesnpaths.check_output_directory(output_dir, ok_if_exists=(ok_if_exists or self.just_do_it))
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        run.warning(None, header="%s RESULTS" % self.similarity_type.upper(), lc='green')
        for report_name in self.results:
            output_path_for_report = J(output_dir, self.method + '_' + report_name)

            utils.store_dict_as_TAB_delimited_file(self.results[report_name], output_path_for_report + '.txt')

            if report_name in self.clusterings:
                with open(output_path_for_report + '.newick', 'w') as f:
                    f.write(self.clusterings[report_name])

            run.info_single('Matrix and clustering of \'%s\' written to output directory' % report_name.replace('_',' '), mc='green')


    def get_genome_names_list_and_genome_source_dict(self):
        def get_names(f):
            d = utils.get_TAB_delimited_file_as_dictionary(f, expected_fields=['name'], indexing_field=-1) if f else {}
            return [line['name'] for line in d.values()]

        names = {
            'fasta': get_names(self.fasta_txt),
            'internal': get_names(self.internal_genomes),
            'external': get_names(self.external_genomes),
        }

        for source1, source2 in combinations(names, 2):
            names_in_both = [n for n in names[source1] if n in names[source2]]
            if len(names_in_both):
                raise ConfigError("Ok, so you provided `%s` and `%s` as sequence sources, but some names from these sources are shared "
                                  "so anvi'o doesn't know how these names should be treated. Here is the list of names that are shared "
                                  "by both: [%s]" % (source1, source2, ', '.join([str(n) for n in names_in_both])))

        self.genome_sources = {}
        for source, names in names.items():
            self.genome_sources.update({name: source for name in names})

        self.genome_names = set(self.genome_sources.keys())

        return self.genome_names, self.genome_sources


    def get_fasta_sequences_dir(self):
        if self.genome_desc:
            if self.representative_method == "Qscore":
                # if representative method is Qscrore, we need SCG information to be included
                self.genome_desc.load_genomes_descriptions(skip_functions=True, init=True)
            else:
                # init=False so that no complaint if there are no SCGs
                self.genome_desc.load_genomes_descriptions(skip_functions=True, init=False)

        self.temp_dir,\
        self.hash_to_name,\
        self.genome_names,\
        self.name_to_temp_path = utils.create_fasta_dir_from_sequence_sources(self.genome_desc, self.fasta_txt)

        return self.temp_dir



class FastANI(GenomeSimilarity):
    def __init__(self, args):
        self.args = args
        self.results = {}

        GenomeSimilarity.__init__(self, args)

        self.similarity_type = 'ANI'
        self.method = 'fastANI'

        self.program = fastani.ManyToMany(args=self.args)

        self.fastANI_sanity_check()


    def fastANI_sanity_check(self):
        pass


    def create_fasta_path_file(self):
        """ Creates a text file with the paths of all the fasta files

        It is the file created here that is passed to fastANI as its --rl and --ql parameters
        """
        file_with_fasta_paths = J(self.temp_dir, 'fasta_paths.txt')
        with open(file_with_fasta_paths, 'w') as f:
            f.write('\n'.join(self.name_to_temp_path.values()) + '\n')

        return file_with_fasta_paths


    def process(self, directory=None):
        self.temp_dir = directory if directory else self.get_fasta_sequences_dir()
        fastANI_input_file = self.create_fasta_path_file()

        self.results = self.program.run_command(
            query_targets=fastANI_input_file,
            reference_targets=fastANI_input_file,
            output_path=J(self.temp_dir, 'output'),
            run_dir=self.temp_dir,
            name_conversion_dict={v: k for k, v in self.name_to_temp_path.items()}
        )

        self.cluster()

        if directory is None:
            shutil.rmtree(self.temp_dir)



class ANI(GenomeSimilarity):
    """This class handles specifically pyANI. See FastANI class for fastani handle"""

    def __init__(self, args):
        self.args = args
        self.results = {}

        GenomeSimilarity.__init__(self, args)

        self.similarity_type = 'ANI'

        self.args.quiet = True
        self.program = pyani.PyANI(self.args)

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.min_alignment_fraction = A('min_alignment_fraction', null)
        self.min_full_percent_identity = A('min_full_percent_identity', null)
        self.significant_alignment_length = A('significant_alignment_length', null)
        self.method = A('method', null)

        self.ANI_sanity_check()


    def ANI_sanity_check(self):
        if self.min_alignment_fraction < 0 or self.min_alignment_fraction > 1:
            raise ConfigError("The minimum alignment fraction must be a value between 0.0 and 1.0. %.2f does "
                              "not work for anvi'o :/" % self.min_alignment_fraction)

        if self.significant_alignment_length and self.significant_alignment_length < 0:
            raise ConfigError("You missed concept :/ Alignment length can't be smaller than 0.")

        if self.significant_alignment_length and not self.min_alignment_fraction:
            raise ConfigError("Using the --significant-alignment-length parameter Without the --min-alignment-fraction "
                              "parameter does not make any sense. But how could you know that unless you read the help "
                              "menu? Yep. Anvi'o is calling the bioinformatics police on you :(")

        if len(self.genome_names) > 50:
            run.warning("You are comparing %d genomes. That's no small task. For context, 10 Prochlorococcus genomes "
                        "using 4 threads took 2m20s on a 2016 Macbook Pro. And 50 genomes took 49m53s. Consider "
                        "instead using fastANI (--program fastANI). It's orders of magnitudes faster." % len(self.genome_names))


    def restore_names_in_dict(self, input_dict):
        """
        Takes dictionary that contains hashes as keys
        and replaces it back to genome names using conversion_dict.

        If value is dict, it calls itself.
        """
        new_dict = {}
        for key, value in input_dict.items():
            if isinstance(value, dict):
                value = self.restore_names_in_dict(value)

            if key in self.hash_to_name:
                new_dict[self.hash_to_name[key]] = value
            else:
                new_dict[key] = value

        return new_dict


    def decouple_weak_associations(self):
        """
        potentially modifies results dict using any of:
            {self.min_alignment_fraction, self.significant_alignment_length, self.min_full_percent_identity}
        """
        # in this list we will keep the tuples of genome-genome associations
        # that need to be set to zero in all result dicts:
        genome_hits_to_zero = []
        num_anvio_will_remove_via_full_percent_identity = 0
        num_anvio_wants_to_remove_via_alignment_fraction = 0
        num_saved_by_significant_length_param = 0

        if self.min_full_percent_identity:
            p = self.results.get('full_percentage_identity')
            if not p:
                raise ConfigError("You asked anvi'o to remove weak hits through the --min-full-percent-identity "
                                  "parameter, but the results dictionary does not contain any information about "
                                  "full percentage identity :/ These are the items anvi'o found instead: '%s'. Please let a "
                                  "developer know about this if this doesn't make any sense." % (', '.join(self.results.keys())))

            for g1 in p:
                for g2 in p:
                    if g1 == g2:
                        continue

                    if float(p[g1][g2]) < self.min_full_percent_identity/100 or float(p[g2][g1]) < self.min_full_percent_identity/100:
                        num_anvio_will_remove_via_full_percent_identity += 1
                        genome_hits_to_zero.append((g1, g2), )

            if len(genome_hits_to_zero):
                g1, g2 = genome_hits_to_zero[0]

                run.warning("THIS IS VERY IMPORTANT! You asked anvi'o to remove any hits between "
                            "two genomes if they had a full percent identity less than '%.2f'. Anvi'o found %d "
                            "such instances between the pairwise comparisons of your %d genomes, and is about "
                            "to set all ANI scores between these instances to 0. For instance, one of your "
                            "genomes, '%s', had a full percentage identity of %.3f relative to '%s', "
                            "another one of your genomes, which is below your threshold, and so the "
                            "ANI scores will be ignored (set to 0) for all downstream "
                            "reports you will find in anvi'o tables and visualizations. Anvi'o "
                            "kindly invites you to carefully think about potential implications of "
                            "discarding hits based on an arbitrary alignment fraction, but does not "
                            "judge you because it is not perfect either." % (self.min_full_percent_identity/100,
                                                                             num_anvio_will_remove_via_full_percent_identity,
                                                                             len(p), g1, float(p[g1][g2]), g2))

        if self.min_alignment_fraction:
            if 'alignment_coverage' not in self.results:
                raise ConfigError("You asked anvi'o to remove weak hits through the --min-alignment-fraction "
                                  "parameter, but the results dictionary does not contain any information about "
                                  "alignment fractions :/ These are the items anvi'o found instead: '%s'. Please let a "
                                  "developer know about this if this doesn't make any sense." % (', '.join(self.results.keys())))

            if self.significant_alignment_length is not None and 'alignment_lengths' not in self.results:
                raise ConfigError("The pyANI results do not contain any alignment lengths data. Perhaps the method you "
                                  "used for pyANI does not produce such data. Well. That's OK. But then you can't use the "
                                  "--significant-alignment-length parameter :/")

            d = self.results['alignment_coverage']
            l = self.results['alignment_lengths']

            for g1 in d:
                for g2 in d:
                    if g1 == g2:
                        continue

                    if float(d[g1][g2]) < self.min_alignment_fraction or float(d[g2][g1]) < self.min_alignment_fraction:
                        num_anvio_wants_to_remove_via_alignment_fraction += 1

                        if self.significant_alignment_length and min(float(l[g1][g2]), float(l[g2][g1])) > self.significant_alignment_length:
                            num_saved_by_significant_length_param += 1
                            continue
                        else:
                            genome_hits_to_zero.append((g1, g2), )

            if num_anvio_wants_to_remove_via_alignment_fraction - num_saved_by_significant_length_param > 0:
                g1, g2 = genome_hits_to_zero[num_anvio_will_remove_via_full_percent_identity]

                if num_saved_by_significant_length_param:
                    additional_msg = "By the way, anvi'o saved %d weak hits because they were longer than the length of %d nts you\
                                      specified using the --significant-alignment-length parameter. " % \
                                            (num_saved_by_significant_length_param, self.significant_alignment_length)
                else:
                    additional_msg = ""

                run.warning("THIS IS VERY IMPORTANT! You asked anvi'o to remove any hits between two genomes if the hit "
                            "was produced by a weak alignment (which you defined as alignment fraction less than '%.2f'). Anvi'o "
                            "found %d such instances between the pairwise comparisons of your %d genomes, and is about "
                            "to set all ANI scores between these instances to 0. For instance, one of your genomes, '%s', "
                            "was %.3f identical to '%s', another one of your genomes, but the aligned fraction of %s to %s was only %.3f "
                            "and was below your threshold, and so the ANI scores will be ignored (set to 0) for all downstream "
                            "reports you will find in anvi'o tables and visualizations. %sAnvi'o kindly invites you "
                            "to carefully think about potential implications of discarding hits based on an arbitrary alignment "
                            "fraction, but does not judge you because it is not perfect either." % \
                                                    (self.min_alignment_fraction,
                                                     num_anvio_wants_to_remove_via_alignment_fraction,
                                                     len(d), g1, float(self.results['percentage_identity'][g1][g2]), g2, g1, g2,
                                                     float(self.results['alignment_coverage'][g1][g2]), additional_msg))

            elif num_saved_by_significant_length_param:
                 run.warning("THIS IS VERY IMPORTANT! You asked anvi'o to remove any hits between two genomes if the hit "
                             "was produced by a weak alignment (which you defined as an alignment fraction less "
                             "than '%.2f'). Anvi'o found %d such instances between the pairwise "
                             "comparisons of your %d genomes, but the --significant-alignment-length parameter "
                             "saved them all, because each one of them were longer than %d nts. So your filters kinda cancelled "
                             "each other out. Just so you know." % \
                                                    (self.min_alignment_fraction,
                                                     num_anvio_wants_to_remove_via_alignment_fraction, len(d),
                                                     self.significant_alignment_length))

        # time to zero those values out:
        genome_hits_to_zero = set(genome_hits_to_zero)
        for report_name in self.results:
            for g1, g2 in genome_hits_to_zero:
                self.results[report_name][g1][g2] = 0
                self.results[report_name][g2][g1] = 0


    def process(self, directory=None):
        self.temp_dir = directory if directory else self.get_fasta_sequences_dir()

        self.results = self.program.run_command(self.temp_dir)
        self.results = self.restore_names_in_dict(self.results)
        self.results = self.compute_additonal_matrices(self.results)
        self.decouple_weak_associations()

        self.cluster()

        if directory is None:
            shutil.rmtree(self.temp_dir)


    def compute_additonal_matrices(self, results):
        # full percentage identity
        try:
            df = lambda matrix_name: pd.DataFrame(results[matrix_name]).astype(float)
            results['full_percentage_identity'] = (df('percentage_identity') * df('alignment_coverage')).to_dict()
        except KeyError:
            # method did not produce percentage_identity score--that's okay, no full percentage
            # identity for you
            pass

        return results


class SourMash(GenomeSimilarity):
    def __init__(self, args):
        GenomeSimilarity.__init__(self, args)

        self.results = {}
        self.clusterings = {}

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.min_similarity = A('min_mash_similarity', null) or 0
        self.kmer_size = A('kmer_size', null) or 13 # lucky number 13
        self.method = 'sourmash'
        self.similarity_type = 'SourMash'

        self.program = sourmash.Sourmash(args)


    def reformat_results(self, results):
        file_to_name = {}
        lines = list(results.keys())
        files = list(results[lines[0]].keys())

        for genome_fasta_path in files:
            genome_fasta_hash = genome_fasta_path[::-1].split(".")[1].split("/")[0][::-1]
            file_to_name[genome_fasta_path] = self.hash_to_name[genome_fasta_hash]

        reformatted_results = {}
        for num, file1 in enumerate(files):
            genome_name_1 = file_to_name[file1]
            reformatted_results[genome_name_1] = {}
            key = lines[num]

            for file2 in files:
                genome_name_2 = file_to_name[file2]
                val = float(results[key][file2])
                reformatted_results[genome_name_1][genome_name_2] = val if val >= self.min_similarity else 0

        return reformatted_results


    def process(self, directory=None):
        self.temp_dir = directory if directory else self.get_fasta_sequences_dir()

        self.results = self.program.process(self.temp_dir, self.name_to_temp_path.values())
        for report in self.results:
            self.results[report] = self.reformat_results(self.results[report])

        self.cluster()

        if directory is None:
            shutil.rmtree(self.temp_dir)


program_class_dictionary = {
    'pyANI': ANI,
    'fastANI': FastANI,
    'sourmash': SourMash,
}
