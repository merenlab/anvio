# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity Generating Retroelements"""

import json
import re
import xml.etree.ElementTree as ET
import csv

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST

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
        
        Running the BLASTn generates an xml file of results.

        Returns
        root : xml file
            An xml of BLASTn results
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
        
        blast = BLAST(output_file, target_fasta =self.input_path, search_program = 'blastn', output_file=blast_output, additional_params = '-dust no')
        blast.evalue = 10 #set Evalue to be same as blastn default
        blast.makedb(dbtype = 'nucl')

        blast.blast(outputfmt = '5', word_size = self.word_size)

        #blast_command = ["blastn", "-query", output_file, "-subject", self.input_path, "-out", blast_output, 
                         #"-word_size", str(self.word_size), "-dust", "no", "-outfmt", "5"]
        #subprocess.run(blast_command)
        tree = ET.parse(blast_output)
        root = tree.getroot()
        return root
        
    #then have new function that takes root as an argument to filter for hits - filter blast hit, with self and root as param. copy in code. 
    def filter_blastn_for_none_identical(self, root):
        """
        This function takes the BLASTn xml output and refines the results to those with less than 100% identity.
        
        Takes the xml file and filters for hits with less than 100% identity, then gives every hit a name
        with its original position in the sequence, counts the bases that are mismatching and on which strand they occur.
        Finally initialises all of these within a dictionary.

        Parameters 
        ==========
        root : xml file 
            An xml of BLASTn results
        
        Returns
        mismatch_hits : dict
            A dictionary of all of the BLASTn hits that are less than 100%
        =======
        
        """
        max_percent_identity = 100
        mismatch_hits = {}

        for iteration in root.findall(".//Iteration"):
            for hit in iteration.findall(".//Hit"):
                for hsp in hit.findall(".//Hsp"):
                    # Get the number of identical positions and their alignment length
                    identical_positions = int(hsp.find('Hsp_identity').text)
                    alignment_length = int(hsp.find('Hsp_align-len').text)

                    percentage_identity = (identical_positions / alignment_length) * 100
                    
                    # Check if the percentage identity is within the threshold (under 100%)
                    if percentage_identity < max_percent_identity:
                        #need to write in the objects for the list
                        section_id = iteration.find('Iteration_query-def').text
                        hsp_num = hsp.find('Hsp_num').text

                        hit_identity = '_'.join([section_id, f'_BLAST_hsp_is_{hsp_num}'])
                        pattern = r"start_bp(\d+)_end_bp(\d+)"

                        # Use re.search to find the pattern in the input string
                        match = re.search(pattern, section_id)

                        # Extract the start and end values from the matched groups
                        query_start_position = int(match.group(1))
                        query_end_position = int(match.group(2))
                        
                        mismatch_hits[hit_identity] = {}

                        qseq = str(hsp.find('Hsp_qseq').text)
                        hseq = str(hsp.find('Hsp_hseq').text)
                        midline = str(hsp.find('Hsp_midline').text)
                        subject_genome_start_position = int(hsp.find('Hsp_hit-from').text)
                        subject_genome_end_position = int(hsp.find('Hsp_hit-to').text)
                        alignment_length = int(hsp.find('Hsp_align-len').text)
                        query_genome_start_position = query_start_position + int(hsp.find('Hsp_query-from').text)
                        query_genome_end_position = query_end_position + int(hsp.find('Hsp_query-to').text)
                        query_frame = int(hsp.find('Hsp_query-frame').text)
                        subject_frame = int(hsp.find('Hsp_hit-frame').text)

                        query_mismatch_positions = []

                        #query_mismatch_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                        #subject_mismatch_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}

                        # Unique characters that may appear in qseq and hseq
                        all_possible_characters = set(qseq + hseq)

                        # Initialize counts with all possible characters
                        query_mismatch_counts = {char: 0 for char in all_possible_characters}
                        subject_mismatch_counts = {char: 0 for char in all_possible_characters}

                        chars_to_skip = [self.skip_dashes]

                        if self.skip_Ns:
                            chars_to_skip.append('N')
                        for idx in range(len(qseq)):
                            if qseq[idx] in chars_to_skip:
                                continue
                            if hseq[idx] in chars_to_skip:
                                continue
                            if qseq[idx] != hseq[idx]:
                                query_mismatch_counts[qseq[idx]]+=1
                                query_mismatch_positions.append(idx)
                                subject_mismatch_counts[hseq[idx]]+=1

                        
                        mismatch_hits[hit_identity] = {
                            'query_seq': qseq,
                            'hit_seq': hseq,
                            'midline': midline,
                            'subject_genome_start_position': subject_genome_start_position,
                            'subject_genome_end_position': subject_genome_end_position,
                            'query_mismatch_counts': query_mismatch_counts,
                            'subject_mismatch_counts': subject_mismatch_counts,
                            'position': query_mismatch_positions,
                            'alignment_length': alignment_length,
                            'query_genome_start_position': query_genome_start_position,
                            'query_genome_end_position': query_genome_end_position,
                            'query_frame': query_frame,
                            'subject_frame': subject_frame
                            }
        return mismatch_hits