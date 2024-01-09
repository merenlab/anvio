# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity Generating Retroelements"""

import subprocess
import json
import re
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
from anvio.errors import ConfigError 

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Katy Lambert-Slosarska"
__email__ = "klambertslosarska@gmail.com"
__status__ = "Development"

class DGR_Finder:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_path = A('input_file')
        self.step = A('step')
        self.word_size = A('word_size')
        self.skip_Ns = A('skip_Ns')
        self.skip_dashes = A('skip_dashes')

        filesnpaths.is_file_fasta_formatted(self.input_path)
        if self.step < 0 or self.word_size < 0:
            raise ConfigError('The step value and/or word size value you are trying to input should be positive.')
        
        self.run.info('Input FASTA file', self.input_path)
        self.run.info('Step size', self.step)
        self.run.info('BLASTn word size', self.word_size)
        self.run.info('Skip "N" characters', self.skip_Ns)
        self.run.info('Skip "-" characters', self.skip_dashes)
    
    def split_sequences(self, start=0):
        """
        This function splits the sequence given into sections of the step value length.
        
        Parameters 
        ==========
        start : integer 
            Start index of the first split (Default: 0)
        
        Returns 
        =======
        section_sequences : list of strings
            A list of the split sequences
        """
        section_sequences = []
        for sequence in SeqIO.parse(self.input_path, "fasta"):
            for i in range(start, len(sequence.seq) - self.step + 1, self.step):
                section = sequence.seq[i:i + self.step]
                section_record = SeqRecord(section, id=f"{sequence.id}_part{i//self.step}_start_bp{i}_end_bp{i + self.step}", description="")
                section_sequences.append(section_record)
                if i + self.step > len(sequence.seq):
                    print(sequence.seq)
        return section_sequences
    
    def run_blastn(self):
        """
        This function runs the BLASTn search of the split sequences against the original input FASTA to find regions of matching nucleotides.
        
        Running the BLASTn generates an xml file of results

        Returns 
        =======
        """
        # split sequences output file name 
        output_file = f"shredded_sequences_step_{self.step}_wordsize_{self.word_size}.fasta"
        #blast output file name
        blast_output = f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml"
        # Start at half the step size of the output file
        overlap_start = self.step // 2
        first_sequences = self.split_sequences(0)
        second_sequences = self.split_sequences(overlap_start)

        all_sequences = first_sequences + second_sequences

        # Write combined sequences to output file
        with open(output_file, "w") as output_handle:
            SeqIO.write(all_sequences, output_handle, "fasta")
        blast_command = ["blastn", "-query", output_file, "-subject", self.input_path, "-out", blast_output, 
                         "-word_size", str(self.word_size), "-dust", "no", "-outfmt", "5"]
        subprocess.run(blast_command)
        tree = ET.parse(blast_output)
        root = tree.getroot()