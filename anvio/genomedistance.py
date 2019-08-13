#!/usr/bin/env python
# -*- coding: utf-8
"""Code for genome distance calculation"""

import shutil


import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.genomedescriptions as genomedescriptions


from anvio.drivers import pyani
from anvio.errors import ConfigError
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
    def __init__(self, args, genome_names, genome_desc, data):
        self.args = args
        self.percent_alignment_threshold = args.percent_alignment_threshold
        self.correlation_threshold = args.correlation_threshold
        self.average_identity_threshold = args.average_identity_threshold
        self.genome_names = genome_names
        self.genomes_dict = genome_desc.genomes
        self.data = data
        self.hash = {}
        self.groups = {}
        self.init_groups()

    def init_groups(self):
        hash = 0
        for name in self.genome_names:
            self.hash[name] = hash
            self.groups[hash] = [name]
            hash = hash + 1

    def group_genomes(self, genome1, genome2):
        from_hash = self.hash[genome1]
        to_hash = self.hash[genome2]
        if from_hash == to_hash:
            return
        for genome in self.groups[from_hash]:
            self.groups[to_hash].append(genome)
            self.hash[genome] = to_hash
            self.groups[from_hash].remove(genome)

    def are_redundant(self, genome1, genome2):
        hash1 = self.hash[genome1]
        hash2 = self.hash[genome2]
        if hash1 == hash2:
            return True
        return False

    def dereplicate(self):
        for genome1 in self.genome_names:
            for genome2 in self.genome_names:
                if genome1 == genome2 or self.are_redundant(genome1, genome2):
                    continue

                if float(self.data['correlations'][genome1][genome2]) >= self.correlation_threshold:
                    if float(self.data['percent_alignment'][genome1][genome2]) >= self.percent_alignment_threshold:
                        if float(self.data['percentage_identity'][genome1][genome2]) >= self.average_identity_threshold:
                            self.group_genomes(genome1, genome2)
        return

    def pick_best_of_two(self, one, two):
        if (one is None or one == []) and (two is None or two == []):
            return None
        elif (one is None or one == []) and len(two) == 1:
            return two[0]
        elif (two is None or two == []) and len(one) == 1:
            return one[0]
        
        best_one = self.pick_representative(one)
        best_two = self.pick_representative(two)
        if (best_one is None or best_one == []) and best_two != [] and best_two is not None:
            return best_two
        elif best_two is None or best_two == [] and best_one != [] and best_one is not None:
            return best_one

        #do some calc
        try:
            score1 = self.genomes_dict[best_one]['percent_completion'] - self.genomes_dict[best_one]['percent_redundancy']
        except:
            raise ConfigError("At least one of your genomes does not contain completion or redundancy estimates. Here is an example: %s." % best_one)
        try:
            score2 = self.genomes_dict[best_two]['percent_completion'] - self.genomes_dict[best_two]['percent_redundancy']
        except:
            raise ConfigError("At least one of your genomes does not contain completion or redundancy estimates. Here is an example: %s." % best_two)

        if score1 > score2:
            return best_one
        elif score2 > score1:
            return best_two
        else:
            len1 = self.genomes_dict[best_one]['total_length']
            len2 = self.genomes_dict[best_two]['total_length']
            if len2 > len1:
                return best_two
            else:
                return best_one

    def pick_representative(self, group):
        if group is None or group == []:
            return None
        elif len(group)== 1:
            return group[0]
        medium = int(len(group) / 2)
        best = self.pick_best_of_two(group[:medium], group[medium:])
        return best

    def get_dereplicated_genome_names(self):
        names = []
        for hash in self.groups.keys():
            group = self.groups[hash]
            if group == []:
                continue
            names.append(self.pick_representative(group))
        return names


class GenomeDistance:
    def __init__(self, args):
        self.args = args
        self.run = run
        self.progress = progress
        if args.internal_genomes is None and args.external_genomes is None:
            self.genome_desc = None
        else:
            self.genome_desc = genomedescriptions.GenomeDescriptions(args, run = terminal.Run(verbose=False))
        #self.fasta_txt = args.fasta_txt or None
        self.fasta_txt = None
        self.hash_to_name = {}
        self.genome_names = set([])

    def get_fasta_sequences_dir(self):
        if self.genome_desc is not None:
            self.genome_desc.load_genomes_descriptions(skip_functions=True) #we want a full init
        temp_dir, hash_to_name, genome_names = utils.create_fasta_dir_from_sequence_sources(self.genome_desc, self.fasta_txt)
        self.hash_to_name = hash_to_name
        self.genome_names = genome_names
        return temp_dir

    def restore_names_in_dict(self, input_dict):
        new_dict = {}
        for key, value in input_dict.items():
            if isinstance(value, dict):
                value = self.restore_names_in_dict(value)
            
            if key in self.hash_to_name:
                new_dict[self.hash_to_name[key]] = value
            else:
                new_dict[key] = value
        return new_dict
    
    def rehash_names(self, names):
        new_dict = dict((a, b) for b, a in self.hash_to_name.items())
        hashes = {}
        for name in names:
            hashes[name] = new_dict[name]
        return hashes

    def retrieve_genome_names(self):
        return self.genome_names

    def remove_redundant_genomes(self, data):
        self.progress.new('Dereplication')
        self.progress.update('Identifying redundant genomes...')
        dict = GenomeDictionary(self.args, self.genome_names, self.genome_desc, data)
        dict.dereplicate()

        self.progress.update('Removing redundant genomes...')
        names = dict.get_dereplicated_genome_names()
        self.progress.end()

        return names

class ANI(GenomeDistance):
    def __init__(self, args):
        GenomeDistance.__init__(self, args)
        self.program = pyani.PyANI(args)
    def process(self, temp=None):
        temp_dir=temp
        if temp is None:
            temp_dir = self.get_fasta_sequences_dir()
        results = self.program.run_command(temp_dir)
        results = self.restore_names_in_dict(results)
        if temp is None:
            shutil.rmtree(temp_dir)
        return results

    def add_default_correlation(self):
        dict = {}
        for name in self.genome_names:
            dict[name] = {}
            for target in self.genome_names:
                dict[name][target] = 1
        return dict
    
    def gen_additional_stats(self, results):
        self.progress.new('Analysis Extension')
        self.progress.update('Generating additional pairwise genome statistics...')
        if 'correlations' not in results.keys():
            self.run.warning("Correlation values were not calculated for your genomes. Please\
                             refer to the log file for more information. Meanwhile, anvi'o\
                             will add default values of 1")
            results['correlations'] = self.add_default_correlation()
        results['percent_alignment'] = {}
        for name in self.genome_names:
            results['percent_alignment'][name] = {}
            for target in self.genome_names:
                alignment = min(float(results['alignment_lengths'][name][target]),
                                float(results['alignment_lengths'][target][name]))
                len1 = int(self.genome_desc.genomes[name]['total_length'])
                len2 = int(self.genome_desc.genomes[target]['total_length'])
                results['percent_alignment'][name][target] = 100 * alignment / min(len1, len2)
        self.progress.end()

        return results
