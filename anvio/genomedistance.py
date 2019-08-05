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


class GenomeDistance:
    def __init__(self, args):
        self.args = args
        if args.internal_genomes is None and args.external_genomes is None:
            self.genome_desc = None
        else:
            self.genome_desc = genomedescriptions.GenomeDescriptions(args, run = terminal.Run(verbose=False))
        self.fasta_txt = args.fasta_txt or None
        self.hash_to_name = {}
        self. genome_names = set([])

    def get_fasta_sequences_dir(self):
        if self.genome_desc is not None:
            self.genome_desc.load_genomes_descriptions(skip_functions=True, init=False)
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
