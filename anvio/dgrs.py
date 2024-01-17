# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity Generating Retroelements"""

import re
import xml.etree.ElementTree as ET
import csv
import os
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import anvio
import anvio.terminal as terminal
import anvio.utils as utils
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
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path= None
        self.fasta_file_path = A('input_file')
        self.step = A('step')
        self.word_size = A('word_size')
        self.skip_Ns = A('skip_Ns')
        self.skip_dashes = A('skip_dashes')
        self.number_of_mismatches = A('number_of_mismatches')
        self.percentage_mismatch = A('percentage_mismatch')
        self.temp_dir = A('temp_dir') or filesnpaths.get_temp_directory_path()

        self.sanity_check()
        if self.step < 0 or self.word_size < 0:
            raise ConfigError('The step value and/or word size value you are trying to input should be positive.')

        self.run.info('Input FASTA file', self.fasta_file_path)
        self.run.info('Step size', self.step)
        self.run.info('BLASTn word size', self.word_size)
        self.run.info('Skip "N" characters', self.skip_Ns)
        self.run.info('Skip "-" characters', self.skip_dashes)
        

    def sanity_check(self):
        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                              "class, not multiple :/")
        if self.fasta_file_path:
            # check  fasta input
            filesnpaths.is_file_fasta_formatted(self.fasta_file_path)

    def get_blast_results(self):
        """
        This function runs the BLASTn search, depending on the input file type, running this against the .
        
        Running the BLASTn generates an xml file of results.

        Returns
        blast_output : xml file
            An xml of BLASTn results
        =======
        """
        #initialise temporary dictionary
        tmp_directory_path = self.temp_dir
        target_file_path = os.path.join(tmp_directory_path,f"input_file.fasta")
        self.run.info('temporary input for blast', target_file_path)

        if self.fasta_file_path or (self.contigs_db_path and not self.profile_db_path):
            shredded_sequence_file = os.path.join(tmp_directory_path,f"shredded_sequences_step_{self.step}_wordsize_{self.word_size}.fasta")
            blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")      
            if self.fasta_file_path:
                os.system(f"cp {self.input_path} {target_file_path}")
            elif self.contigs_db_path:
                utils.export_sequence_from_contigs_db(self.contigs_db_path, target_file_path)
            # Start at half the step size of the output file
            overlap_start = self.step // 2
            first_sequences = self.split_sequences()
            second_sequences = self.split_sequences(overlap_start)

            all_sequences = first_sequences + second_sequences

            # Write combined sequences to output file
            with open(shredded_sequence_file, "w") as output_handle:
                SeqIO.write(all_sequences, output_handle, "fasta") 
            
            blast = BLAST(shredded_sequence_file, target_fasta =target_file_path, search_program = 'blastn', output_file=blast_output, additional_params = '-dust no')
            blast.evalue = 10 #set Evalue to be same as blastn default
            blast.makedb(dbtype = 'nucl')
            blast.blast(outputfmt = '5', word_size = self.word_size)

        elif self.contigs_db_path and self.profile_db_path:
            #contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
                #self.splits_basic_info = contigs_db.db.smart_get(t.splits_info_table_name, column='split', data=split_names)
                #self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
                #self.genes_are_called_in_contigs_db = contigs_db.meta['genes_are_called']
                #self.genes_annotated_with_functions_in_contigs_db = contigs_db.meta['gene_function_sources'] is not None and len(contigs_db.meta['gene_function_sources']) > 0
                #contigs_db.disconnect()
            #NEED SNV WINDOW MAKER func here that takes SNV info from .db using tables? or tableops? 
            #then take in these windows outputs nd blast to give blast_output in xaml file. but need to original sequence to comapre the SNV windows to?
            pass

        return(blast_output)
                


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
        #initialise temporary dictionary
        tmp_directory_path = filesnpaths.get_temp_directory_path()
        # split sequences output file name 
        shredded_sequence_file = os.path.join(tmp_directory_path,f"shredded_sequences_step_{self.step}_wordsize_{self.word_size}.fasta")
        #blast output file name
        blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")
        target_file_path = os.path.join(tmp_directory_path,f"input_file.fasta")
        print(f"cp {self.input_path} {target_file_path}")
        os.system(f"cp {self.input_path} {target_file_path}")
        self.run.info('temporary input for blast', target_file_path) 

        # Start at half the step size of the output file
        overlap_start = self.step // 2
        first_sequences = self.split_sequences()
        second_sequences = self.split_sequences(overlap_start)

        all_sequences = first_sequences + second_sequences

        # Write combined sequences to output file
        with open(shredded_sequence_file, "w") as output_handle:
            SeqIO.write(all_sequences, output_handle, "fasta")
        #need a temporary directory where intermediate files are written, to call on them. 
        
        blast = BLAST(shredded_sequence_file, target_fasta =target_file_path, search_program = 'blastn', output_file=blast_output, additional_params = '-dust no')
        blast.evalue = 10 #set Evalue to be same as blastn default
        blast.makedb(dbtype = 'nucl')
        blast.blast(outputfmt = '5', word_size = self.word_size)

        #blast_command = ["blastn", "-query", output_file, "-subject", self.input_path, "-out", blast_output, 
                         #"-word_size", str(self.word_size), "-dust", "no", "-outfmt", "5"]
        #subprocess.run(blast_command)
        return blast_output
    #need to return blast_output which is a path so next function calls the path. 
        
    #then have new function that takes root as an argument to filter for hits - filter blast hit, with self and root as param. copy in code. 
    def filter_blastn_for_none_identical(self, blast_output):
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
        tree = ET.parse(blast_output)
        root = tree.getroot()
        
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
    
    def filter_for_TR_VR(self, mismatch_hits):
        """
        This function takes the none identical hits of the BLASTn and filters for template and variable regions.

        This works by filtering for sequences that have an overrepresentation of one base that is mismatching and a certain number 
        one type of base mismatching within the sequence, defined by the number of mismatches argument. 

        Parameters 
        ==========
        mismatch_hits : dict
            A dictionary of all of the BLASTn hits that are less than 100%
        
        Returns
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions
        
        """
        num_DGR = 0

        #possible DGR dictionary 
        DGRs_found_dict = {}

        for sequence_component, hit_data in mismatch_hits.items():
            query_mismatch_counts = hit_data['query_mismatch_counts']
            subject_mismatch_counts = hit_data['subject_mismatch_counts']
            position = hit_data['position']
            subject_genome_start_position = hit_data['subject_genome_start_position']
            subject_genome_end_position = hit_data['subject_genome_end_position']
            alignment_length = hit_data['alignment_length']
            subject_sequence = Seq(hit_data['hit_seq'])
            midline = hit_data['midline']
            query_sequence = Seq(hit_data['query_seq'])
            shredded_sequence_name = sequence_component
            query_genome_start_position = hit_data['query_genome_start_position']
            query_genome_end_position = hit_data['query_genome_end_position']
            query_frame = hit_data['query_frame']
            subject_frame = hit_data['subject_frame']
            TR_sequence_found = None
            VR_sequence_found = None

        # get number of mismatches
            mismatch_length_bp = len(position)

            # if num of mismatches = 0, skip DGR search sanity check 
            if mismatch_length_bp == 0:
                continue
                #old code mismatch_dict[hit_id]['is_DGR'] = False
            else:
                # Calculate the percentage identity of each alignment
                is_TR = False
                for letter, count in query_mismatch_counts.items():
                    percentage_of_mismatches = (count / mismatch_length_bp)
                    if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches): 
                        #make the nums changeable params w/ sanity check that percentage_of_mismatches > 0.5 and mismatch_length > 0
                        is_TR = True
                        #need to check if the new TR youre looping through exsists in the DGR_found_dict, compare start stop position (likely not equal) 
                        #take longest one, bit like the FIlter code, replace sequence with longest TR. Check if VR already exsists, 
                        num_DGR += 1
                        #creates an empty dict, that has itself empty dicts frothe VRs so you can fill it and create a new key 
                        DGRs_found_dict[f'DGR_{num_DGR:03d}'] = {'VRs':{'VR1':{}}}
                        if letter == 'T':
                            #this section needs work, doesnt change T to A or reverse midline and reverse complement the sequences :(
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence'] = str(query_sequence.reverse_complement())
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence'] = str(subject_sequence.reverse_complement())
                            #overwrite midline string 
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['midline'] =  ''.join(reversed(midline))
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['base'] = letter.replace('T', 'A')
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_reverse_complement'] = True
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence_found'] = 'query'
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence_found'] = 'subject'
                        else:
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence'] = str(query_sequence)
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence'] = str(subject_sequence)
                            #overwrite midline string 
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['midline'] = midline
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['base'] = letter
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_reverse_complement'] = False
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence_found'] = 'query'
                            DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence_found'] = 'subject'               
                        
                        DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_start_position'] = query_genome_start_position
                        DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_end_position'] = query_genome_end_position
                        DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_start_position'] = subject_genome_start_position
                        DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_end_position'] = subject_genome_end_position
                        DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['percentage_of_mismatches'] = percentage_of_mismatches

                if not is_TR:
                    # Calculate the percentage identity of each alignment
                    for letter, count in subject_mismatch_counts.items():
                            percentage_of_mismatches = (count / mismatch_length_bp)
                            if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches): 
                                #make the nums changeable params w/ sanity check that percentage_of_mismatches > 0.5 and mismatch_length > 0
                                is_TR = True
                                num_DGR += 1
                                #creates an empty dict, that has itself empty dicts frothe VRs so you can fill it and create a new key 
                                DGRs_found_dict[f'DGR_{num_DGR:03d}'] = {'VRs':{'VR1':{}}}
                                if letter == 'T':
                                    #this section needs work, doesnt change T to A or reverse midline and reverse complement the sequences :(
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence'] = str(subject_sequence.reverse_complement())
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence'] = str(query_sequence.reverse_complement())
                                    #overwrite midline string 
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['midline'] =  ''.join(reversed(midline))
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['base'] = letter.replace('T', 'A')
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_reverse_complement'] = True
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence_found'] = 'subject'
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence_found'] = 'query'
                                else:
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence'] = str(subject_sequence)
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence'] = str(query_sequence)
                                    #overwrite midline string 
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['midline'] = midline
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['base'] = letter
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_reverse_complement'] = False
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_sequence_found'] = 'subject'
                                    DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_sequence_found'] = 'query'

                                DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_start_position'] = subject_genome_start_position
                                DGRs_found_dict[f'DGR_{num_DGR:03d}']['TR_end_position'] = subject_genome_end_position
                                DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_start_position'] = query_genome_start_position
                                DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['VR_end_position'] = query_genome_end_position
                                DGRs_found_dict[f'DGR_{num_DGR:03d}']['VRs']['VR1']['percentage_of_mismatches'] = percentage_of_mismatches
                            
        print(f'number of DGRs is {num_DGR}')
        if anvio.DEBUG:
            self.run.warning(f"The temp directory, '{self.temp_dir}', is kept. Don't forget to clean it up later!", header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (use `--debug` to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(self.temp_dir)
        return DGRs_found_dict
    
    def create_found_tr_vr_csv(self, DGRs_found_dict):
        """
        This function creates a csv tabular format of the template and variable regions that are found from this tool.
        Parameters 
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions
        
        Returns
        =======
        : csv
            A csv tabular file containing the template and variable regions
        
        """
        base_input_name = None
        if self.fasta_file_path:
             base_input_name = os.path.basename(self.fasta_file_path)
        elif self.contigs_db_path:
             base_input_name = os.path.basename(self.contigs_db_path)

        csv_file_path = f'DGRs_found_from_{base_input_name}_percentage_{self.percentage_mismatch}_number_mismatches_{self.number_of_mismatches}.csv'
        with open(csv_file_path, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            
            # Write header
            csv_writer.writerow(["DGR", "VR_sequence", "Midline","VR_sequence_found", "VR_start_position", "VR_end_position", "Mismatch %",
                                "TR_sequence", "Base","TR_sequence_found", "Reverse Complement", "TR_start_position", "TR_end_position"])
            
            # Write data
            for dgr, info in DGRs_found_dict.items():
                vr_data = info['VRs']['VR1']
                #Print vr_data for debugging
                print(f'DGR: {dgr}, vr_data: {vr_data}')
                csv_writer.writerow([dgr, vr_data['VR_sequence'], vr_data['midline'], vr_data['VR_sequence_found'], vr_data['VR_start_position'], vr_data['VR_end_position'],
                                    vr_data['percentage_of_mismatches'], info['TR_sequence'], info['base'], info['TR_sequence_found'], info['TR_reverse_complement'],
                                    info['TR_start_position'], info['TR_end_position']])
                return csv_file_path
            