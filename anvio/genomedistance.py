#!/usr/bin/env python
# -*- coding: utf-8
"""Code for genome distance calculation"""

import shutil
import pandas as pd

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.genomedescriptions as genomedescriptions

from itertools import combinations

from anvio.errors import ConfigError
from anvio.drivers import pyani, sourmash
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.tables.miscdata import TableForLayerOrders

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Mahmoud Yousef"
__email__ = "mahmoudyousef@uchicago.edu"

run = terminal.Run()
progress = terminal.Progress()


class GenomeDictionary:
    def __init__(self, args, genome_names, genome_desc, fasta_txt, distance_matrix):
        self.args = args
        self.distance_threshold = self.args.distance_threshold
        self.method = self.args.representative_method

        self.genome_names = genome_names
        self.distance_matrix = distance_matrix

        self.clusters = {}
        self.cluster_report = {}
        self.genomes_info_dict = {}
        self.representative_names = {}
        self.genome_name_to_cluster_name = {}

        self.init_clusters()
        self.populate_genomes_info_dict(genome_desc, fasta_txt)


    def init_clusters(self):
        """Each genome starts in its own cluster, i.e. as many clusters as genomes"""
        for count, genome_name in enumerate(self.genome_names, start=1):
            cluster_name = 'cluster_%05d' % count
            self.genome_name_to_cluster_name[genome_name] = cluster_name
            self.clusters[cluster_name] = [genome_name]


    def populate_genomes_info_dict(self, genome_desc, fasta_txt):
        full_dict = {}

        if genome_desc:
            full_dict = genome_desc.genomes

        if fasta_txt:
            fastas = utils.get_TAB_delimited_file_as_dictionary(fasta_txt, expected_fields=['name', 'path'], only_expected_fields=True)

            for name in fastas:
                if name in full_dict.keys(): #redundant name
                    continue

                full_dict[name] = {}
                full_dict[name]['percent_completion'] = 0
                full_dict[name]['percent_redundancy'] = 0
                full_dict[name]['total_length'] = sum(utils.get_read_lengths_from_fasta(fastas[name]['path']).values())

        self.genomes_info_dict = full_dict


    def update_clusters(self, genome1, genome2):
        from_cluster = self.genome_name_to_cluster_name[genome1]
        to_cluster = self.genome_name_to_cluster_name[genome2]

        if from_cluster == to_cluster:
            return

        for genome in self.clusters[from_cluster]:
            self.clusters[to_cluster].append(genome)
            self.genome_name_to_cluster_name[genome] = to_cluster
            self.clusters[from_cluster].remove(genome)


    def are_redundant(self, genome1, genome2):
        cluster1 = self.genome_name_to_cluster_name[genome1]
        cluster2 = self.genome_name_to_cluster_name[genome2]
        return True if cluster1 == cluster2 else False


    def gen_cluster_report(self):
        if not any([self.clusters, self.representative_names]):
            raise ConfigError("gen_cluster_report :: Run dereplicate() before trying to generate a cluster report")

        self.cluster_report = {
            cluster_name: {
                'representative': self.representative_names[cluster_name],
                'genomes': self.clusters[cluster_name]
            } for cluster_name in self.clusters
        }


    def rename_clusters(self):
        new_cluster_names = ['cluster_%05d' % count for count in range(1, len(self.clusters) + 1)]
        self.clusters = {new_cluster_name: genomes for new_cluster_name, genomes in zip(new_cluster_names, self.clusters.values())}


    def dereplicate(self):
        genome_pairs = combinations(self.genome_names, 2)
        for genome1, genome2 in genome_pairs:
            if genome1 == genome2 or self.are_redundant(genome1, genome2):
                continue

            distance = float(self.distance_matrix[genome1][genome2])
            if distance > self.distance_threshold:
                self.update_clusters(genome1, genome2)

        # remove empty clusters and rename so that names are sequential
        self.clusters = {cluster: genomes for cluster, genomes in self.clusters.items() if genomes}
        self.rename_clusters()

        self.representative_names = self.get_representative_for_each_cluster()
        self.gen_cluster_report()

        return self.cluster_report


    def pick_best_of_two(self, one, two):
        if not one and not two:
            return None
        elif not one and len(two) == 1:
            return two[0]
        elif not two and len(one) == 1:
            return one[0]

        best_one = self.pick_representative(one)
        best_two = self.pick_representative(two)

        if not best_one and best_two:
            return best_two
        elif not best_two and best_one:
            return best_one

        try:
            score1 = self.genomes_info_dict[best_one]['percent_completion'] - self.genomes_info_dict[best_one]['percent_redundancy']
        except:
            raise ConfigError("At least one of your genomes does not contain completion or redundancy estimates. Here is an example: %s." % best_one)
        try:
            score2 = self.genomes_info_dict[best_two]['percent_completion'] - self.genomes_info_dict[best_two]['percent_redundancy']
        except:
            raise ConfigError("At least one of your genomes does not contain completion or redundancy estimates. Here is an example: %s." % best_two)

        if score1 > score2:
            return best_one
        elif score2 > score1:
            return best_two
        else:
            len1 = self.genomes_info_dict[best_one]['total_length']
            len2 = self.genomes_info_dict[best_two]['total_length']

            if len2 > len1:
                return best_two
            else:
                return best_one


    def pick_representative(self, cluster):
        if not cluster:
            return None
        elif len(cluster) == 1:
            return cluster[0]

        medium = int(len(cluster) / 2)
        best = self.pick_best_of_two(cluster[:medium], cluster[medium:])

        return best


    def pick_longest_representative(self, cluster):
        max_name = cluster[0]
        max_val = self.genomes_info_dict[max_name]['total_length']

        for name in cluster[1:]:
            val = self.genomes_info_dict[name]['total_length']

            if val > max_val:
                max_name = name
                max_val = val

        return max_name


    def pick_closest_distance(self, cluster):
        d = {}
        for name in cluster:
            d[name] = 0

            for target in cluster:
                d[name] += float(self.distance_matrix[name][target])

        new_dict = {}
        for name, val in d.items():
            new_dict[val] = name

        max_val = max(new_dict.keys())

        return new_dict[max_val]


    def get_representative_for_each_cluster(self):
        representative_names = {}
        for h in self.clusters.keys():
            cluster = self.clusters[h]

            if cluster == []:
                continue

            if self.method == "Qscore":
                name = self.pick_representative(cluster)
            elif self.method == "length":
                name = self.pick_longest_representative(cluster)
            elif self.method == "distance":
                name = self.pick_closest_distance(cluster)

            representative_names[h] = name

        return representative_names



class GenomeDistance:
    def __init__(self, args):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.fasta_txt = A('fasta_text_file', null)
        self.internal_genomes = A('internal_genomes', null)
        self.external_genomes = A('external_genomes', null)

        if (self.internal_genomes or self.external_genomes):
            self.genome_desc = genomedescriptions.GenomeDescriptions(args, run = terminal.Run(verbose=False))
        else:
            self.genome_desc = None

        self.hash_to_name = {}
        self.name_to_temp_path = {}
        self.genome_names = set([])


    def get_fasta_sequences_dir(self):
        if self.genome_desc is not None:
            self.genome_desc.load_genomes_descriptions(skip_functions=True)

        temp_dir, hash_to_name, genome_names, name_to_temp_path = utils.create_fasta_dir_from_sequence_sources(self.genome_desc, fasta_txt=self.fasta_txt)

        self.hash_to_name = hash_to_name
        self.genome_names = genome_names
        self.name_to_temp_path = name_to_temp_path
        return temp_dir


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


    def remove_redundant_genomes(self, distance_matrix):
        self.progress.new('Dereplication')
        self.progress.update('Generating cluster info report')

        genomes_dict = GenomeDictionary(self.args, self.genome_names, self.genome_desc, self.fasta_txt, distance_matrix)
        cluster_report = genomes_dict.dereplicate()

        self.progress.end()

        return cluster_report



class ANI(GenomeDistance):
    def __init__(self, args):
        self.args = args
        self.results = {}

        GenomeDistance.__init__(self, args)

        self.args.quiet = True
        self.program = pyani.PyANI(self.args)

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.min_alignment_fraction = A('min_alignment_fraction', null)
        self.min_full_percent_identity = A('min_full_percent_identity', null)
        self.significant_alignment_length = A('significant_alignment_length', null)


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
                raise ConfigError("You asked anvi'o to remove weak hits through the --min-full-percent-identity\
                                   parameter, but the results dictionary does not contain any information about\
                                   full percentage identity :/ These are the items anvi'o found instead: '%s'. Please let a\
                                   developer know about this if this doesn't make any sense." % (', '.join(self.results.keys())))

            for g1 in p:
                for g2 in p:
                    if g1 == g2:
                        continue

                    if float(p[g1][g2]) < self.min_full_percent_identity/100 or float(p[g2][g1]) < self.min_full_percent_identity/100:
                        num_anvio_will_remove_via_full_percent_identity += 1
                        genome_hits_to_zero.append((g1, g2), )

            if len(genome_hits_to_zero):
                g1, g2 = genome_hits_to_zero[0]

                run.warning("THIS IS VERY IMPORTANT! You asked anvi'o to remove any hits between\
                             two genomes if they had a full percent identity less than '%.2f'. Anvi'o found %d\
                             such instances between the pairwise comparisons of your %d genomes, and is about\
                             to set all ANI scores between these instances to 0. For instance, one of your \
                             genomes, '%s', had a full percentage identity of %.3f relative to '%s',\
                             another one of your genomes, which is below your threshold, and so the\
                             ANI scores will be ignored (set to 0) for all downstream\
                             reports you will find in anvi'o tables and visualizations. Anvi'o\
                             kindly invites you to carefully think about potential implications of\
                             discarding hits based on an arbitrary alignment fraction, but does not\
                             judge you because it is not perfect either." % (self.min_full_percent_identity/100,
                                                                             num_anvio_will_remove_via_full_percent_identity,
                                                                             len(p), g1, float(p[g1][g2]), g2))

        if self.min_alignment_fraction:
            if 'alignment_coverage' not in self.results:
                raise ConfigError("You asked anvi'o to remove weak hits through the --min-alignment-fraction\
                                   parameter, but the results dictionary does not contain any information about\
                                   alignment fractions :/ These are the items anvi'o found instead: '%s'. Please let a\
                                   developer know about this if this doesn't make any sense." % (', '.join(self.results.keys())))

            if self.significant_alignment_length is not None and 'alignment_lengths' not in self.results:
                raise ConfigError("The pyANI results do not contain any alignment lengths data. Perhaps the method you\
                                   used for pyANI does not produce such data. Well. That's OK. But then you can't use the\
                                   --significant-alignment-length parameter :/")

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
                    additional_msg = "By the way, anvi'o saved %d weak hits becasue they were longer than the length of %d nts you\
                                      specified using the --significant-alignment-length parameter. " % \
                                            (num_saved_by_significant_length_param, self.significant_alignment_length)
                else:
                    additional_msg = ""

                run.warning("THIS IS VERY IMPORTANT! You asked anvi'o to remove any hits between two genomes if the hit\
                             was produced by a weak alignment (which you defined as alignment fraction less than '%.2f'). Anvi'o\
                             found %d such instances between the pairwise comparisons of your %d genomes, and is about\
                             to set all ANI scores between these instances to 0. For instance, one of your genomes, '%s',\
                             was %.3f identical to '%s', another one of your genomes, but the aligned fraction of %s to %s was only %.3f\
                             and was below your threshold, and so the ANI scores will be ignored (set to 0) for all downstream\
                             reports you will find in anvi'o tables and visualizations. %sAnvi'o kindly invites you\
                             to carefully think about potential implications of discarding hits based on an arbitrary alignment\
                             fraction, but does not judge you because it is not perfect either." % \
                                                    (self.min_alignment_fraction,
                                                     num_anvio_wants_to_remove_via_alignment_fraction,
                                                     len(d), g1, float(self.results['percentage_identity'][g1][g2]), g2, g1, g2,
                                                     float(self.results['alignment_coverage'][g1][g2]), additional_msg))

            elif num_saved_by_significant_length_param:
                 run.warning("THIS IS VERY IMPORTANT! You asked anvi'o to remove any hits between two genomes if the hit\
                              was produced by a weak alignment (which you defined as an alignment fraction less\
                              than '%.2f'). Anvi'o found %d such instances between the pairwise\
                              comparisons of your %d genomes, but the --significant-alignment-length parameter\
                              saved them all, because each one of them were longer than %d nts. So your filters kinda cancelled\
                              each other out. Just so you know." % \
                                                    (self.min_alignment_fraction,
                                                     num_anvio_wants_to_remove_via_alignment_fraction, len(d),
                                                     self.significant_alignment_length))

        # time to zero those values out:
        genome_hits_to_zero = set(genome_hits_to_zero)
        for report_name in self.results:
            for g1, g2 in genome_hits_to_zero:
                self.results[report_name][g1][g2] = 0
                self.results[report_name][g2][g1] = 0


    def process(self, temp=None):
        temp_dir = temp if temp else self.get_fasta_sequences_dir()

        self.results = self.program.run_command(temp_dir)
        self.results = self.restore_names_in_dict(self.results)
        self.results = self.compute_additonal_matrices(self.results)
        self.decouple_weak_associations()

        if temp is None:
            shutil.rmtree(temp_dir)


    def compute_additonal_matrices(self, results):
        # full percentage identity
        df = lambda matrix_name: pd.DataFrame(results[matrix_name]).astype(float)
        results['full_percentage_identity'] = (df('percentage_identity') * df('alignment_coverage')).to_dict()

        return results



class SourMash(GenomeDistance):
    def __init__(self, args):
        GenomeDistance.__init__(self, args)

        self.results = {}
        self.program = sourmash.Sourmash(args)
        self.min_distance = args.min_mash_distance if 'min_mash_distance' in vars(args) else 0


    def reformat_results(self, results):
        file_to_name = {}
        lines = list(results.keys())
        files = list(results[lines[0]].keys())

        for genome_fasta_path in files:
            genome_fasta_hash = genome_fasta_path[::-1].split(".")[1].split("/")[0]
            genome_fasta_hash = genome_fasta_hash[::-1]
            file_to_name[genome_fasta_path] = self.hash_to_name[genome_fasta_hash]

        reformatted_results = {}
        for num, file1 in enumerate(files):
            genome_name_1 = file_to_name[file1]
            reformatted_results[genome_name_1] = {}
            key = lines[num]

            for file2 in files:
                genome_name_2 = file_to_name[file2]
                val = float(results[key][file2])
                reformatted_results[genome_name_1][genome_name_2] = val if val >= self.min_distance else 0

        return reformatted_results


    def process(self, temp=None):
        temp_dir = temp if temp else self.get_fasta_sequences_dir()

        self.results['mash_distance'] = self.program.process(temp_dir, self.name_to_temp_path.values())
        self.results['mash_distance'] = self.reformat_results(self.results['mash_distance'])

        if temp is None:
            shutil.rmtree(temp_dir)
