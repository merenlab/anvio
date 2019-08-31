#!/usr/bin/env python
# -*- coding: utf-8
"""Code for genome distance calculation"""

import shutil

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
    def __init__(self, args, genome_names, genome_desc, fasta_txt, data):
        self.args = args
        self.distance_threshold = args.distance_threshold
        self.method = args.representative_method
        self.genome_names = genome_names
        self.genomes_dict = {}
        self.data = data
        self.hash = {}
        self.groups = {}
        self.init_groups()
        self.init_full_genome_dict(genome_desc, fasta_txt)

    def init_groups(self):
        hash = 0
        for name in self.genome_names:
            self.hash[name] = hash
            self.groups[hash] = [name]
            hash = hash + 1

    def init_full_genome_dict(self, genome_desc, fasta_txt):
        full_dict = {}
        if genome_desc is not None:
            full_dict = genome_desc.genomes
        if fasta_txt is not None:
            fastas = utils.get_TAB_delimited_file_as_dictionary(fasta_txt, expected_fields=['name', 'path'], only_expected_fields=True)
            for name in fastas:
                if name in full_dict.keys(): #redundant name
                    continue
                full_dict[name] = {}
                full_dict[name]['percent_completion'] = 0
                full_dict[name]['percent_redundancy'] = 0
                full_dict[name]['total_length'] = 0
        self.genomes_dict = full_dict



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
        genome_pairs = combinations(self.genome_names, 2)
        for pair in genome_pairs:
            genome1 = pair[0]
            genome2 = pair[1]
            if genome1 == genome2 or self.are_redundant(genome1, genome2):
                continue

            distance = float(self.data[genome1][genome2])
            if distance > self.distance_threshold:
                self.group_genomes(genome1, genome2)

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

    def pick_longest_rep(self, group):
        max_name = group[0]
        max_val = self.genomes_dict[max_name]['total_length']

        for name in group[1:]:
            val = self.genomes_dict[name]['total_length']
            if val > max_val:
                max_name = name
                max_val = val

        return max_name

    def pick_closest_distance(self, group):
        dict = {}
        for name in group:
            dict[name] = 0
            for target in group:
                dict[name] += float(self.data[name][target])
        new_dict = {}
        for name, val in dict.items():
            new_dict[val] = name
        max_val = max(new_dict.keys())
        return new_dict[max_val]

    def get_dereplicated_genome_names(self):
        names = []
        for hash in self.groups.keys():
            group = self.groups[hash]
            if group == []:
                continue
            if self.method == "Qscore":
                name = self.pick_representative(group)
            elif self.method == "length":
                name = self.pick_longest_rep(group)
            elif self.method == "distance":
                name = self.pick_closest_distance(group)
            names.append(name)
        return names

    def get_all_redundant_groups(self, names):
        dict = {}
        for name in names:
            hash = self.hash[genome]
            dict[name] = self.groups[hash]
        return dict

class GenomeDistance:
    def __init__(self, args):
        self.args = args
        self.run = run
        self.progress = progress
        if args.internal_genomes is None and args.external_genomes is None:
            self.genome_desc = None
        else:
            self.genome_desc = genomedescriptions.GenomeDescriptions(args, run = terminal.Run(verbose=False))
        self.fasta_txt = args.fasta_text_file if 'fasta_text_file' in vars(args) else None
        self.hash_to_name = {}
        self.genome_names = set([])
        self.temp_paths = set([])
        self.dict = {}

    def get_fasta_sequences_dir(self):
        if self.genome_desc is not None:
            self.genome_desc.load_genomes_descriptions(skip_functions=True)
        temp_dir, hash_to_name, genomes = utils.create_fasta_dir_from_sequence_sources(self.genome_desc, fasta_txt=self.fasta_txt)
        self.hash_to_name = hash_to_name
        self.genome_names = genomes[0]
        self.temp_paths = genomes[1]
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
        self.dict = GenomeDictionary(self.args, self.genome_names, self.genome_desc, self.fasta_txt, data)
        self.dict.dereplicate()

        self.progress.update('Removing redundant genomes...')
        names = self.dict.get_dereplicated_genome_names()
        self.progress.end()

        return names

    def retrieve_genome_groups(self, names):
        return self.dict.get_all_redundant_groups(names)

class ANI(GenomeDistance):
    def __init__(self, args):
        GenomeDistance.__init__(self, args)
        self.program = pyani.PyANI(args)
        self.min_coverage = args.min_alignment_coverage if 'min_alignment_coverage' in vars(args) else 0
    
    def collapse_results(self, results):
        self.progress.new('PyANI')
        self.progress.update('Collapsing results into distance matrix...')
        matrix = {}
        for name in self.genome_names:
            matrix[name] = {}
            for target in self.genome_names:
                if float(results['alignment_coverage'][name][target]) < self.min_coverage:
                    matrix[name][target] = 0
                else:
                    matrix[name][target] = float(results['percentage_identity'][name][target])

        self.progress.end()
        return matrix


    def process(self, temp=None):
        temp_dir=temp
        if temp is None:
            temp_dir = self.get_fasta_sequences_dir()
        results = self.program.run_command(temp_dir)
        results = self.restore_names_in_dict(results)
        results = self.collapse_results(results)
        if temp is None:
            shutil.rmtree(temp_dir)
        return results


class SourMash(GenomeDistance):
    def __init__(self, args):
        GenomeDistance.__init__(self, args)
        self.program = sourmash.Sourmash(args)
        self.min_distance = args.min_mash_distance if 'min_mash_distance' in vars(args) else 0
    
    def reformat_results(self, dict):
        file_hash = {}
        lines = list(dict.keys())
        files = list(dict[lines[0]].keys())
        for file in files:
            hash = file[::-1].split(".")[1].split("/")[0]
            hash = hash[::-1]
            file_hash[file] = self.hash_to_name[hash]
        results = {}
        for num, file1 in enumerate(files):
            name1 = file_hash[file1]
            key = lines[num]
            results[name1] = {}
            for file2 in files:
                name2 = file_hash[file2]
                val = float(dict[key][file2])
                if val < self.min_distance:
                    results[name1][name2] = 0
                else:
                    results[name1][name2] = val
        return results

    def process(self, temp=None):
        temp_dir=temp
        if temp is None:
            temp_dir = self.get_fasta_sequences_dir()
        results = self.program.process(temp_dir, list(self.temp_paths))
        results = self.reformat_results(results)
        if temp is None:
            shutil.rmtree(temp_dir)
        return results
