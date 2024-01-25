# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity Generating Retroelements"""

import re
import xml.etree.ElementTree as ET
import csv
import os
import shutil
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths
import anvio.tables as t

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.variabilityops import NucleotidesEngine

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Katy Lambert-Slosarska"
__email__ = "klambertslosarska@gmail.com"
__status__ = "Development"

run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)

class DGR_Finder:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path= A('profile_db')
        self.fasta_file_path = A('input_file')
        self.step = A('step')
        self.word_size = A('word_size')
        self.skip_Ns = A('skip_Ns')
        self.skip_dashes = A('skip_dashes')
        self.number_of_mismatches = A('number_of_mismatches')
        self.percentage_mismatch = A('percentage_mismatch')
        self.temp_dir = A('temp_dir') or filesnpaths.get_temp_directory_path()
        self.min_dist_bw_snvs = A('distance_between_snv')
        self.variable_buffer_length = A('variable_buffer_length')
        self.departure_from_reference_percentage = A('departure_from_reference_percentage')

        self.sanity_check()

        if self.fasta_file_path:
            self.run.info('Input FASTA file', self.fasta_file_path)
        if self.contigs_db_path:
            self.run.info('Contigs.db', self.contigs_db_path)
        if self.profile_db_path:
            self.run.info('Profile.db', self.profile_db_path)
        if self.fasta_file_path or self.contigs_db_path and not self.profile_db_path:
            self.run.info('Step size', self.step)
        self.run.info('BLASTn word size', self.word_size)
        self.run.info('Skip "N" characters', self.skip_Ns)
        self.run.info('Skip "-" characters', self.skip_dashes)
        if self.profile_db_path and self.contigs_db_path:
            self.run.info('Minimum distance between SNVs', self.min_dist_bw_snvs)
            self.run.info('Variable buffer length', self.variable_buffer_length)
            self.run.info('Departure from reference percentage', self.departure_from_reference_percentage)
    def sanity_check(self):
        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                              "class, not multiple :/")
        if self.fasta_file_path:
            # check fasta input
            filesnpaths.is_file_fasta_formatted(self.fasta_file_path)

        if self.step < 0 or self.word_size < 0:
            raise ConfigError('The step value and/or word size value you are trying to input should be positive integer.')

        if self.variable_buffer_length < 0:
            raise ConfigError('The variable buffer length value you are trying to input should be positive integer.')

        if self.departure_from_reference_percentage < 0:
            raise ConfigError('The departure from reference percentage value you are trying to input should be a positive decimal number.')

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
        self.target_file_path = os.path.join(tmp_directory_path,f"input_file.fasta")
        self.run.info('Temporary input for blast', self.target_file_path)

        if self.fasta_file_path or (self.contigs_db_path and not self.profile_db_path):
            shredded_sequence_file = os.path.join(tmp_directory_path,f"shredded_sequences_step_{self.step}_wordsize_{self.word_size}.fasta")
            blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")      
            if self.fasta_file_path:
                os.system(f"cp {self.fasta_file_path} {self.target_file_path}")
            elif self.contigs_db_path:
                utils.export_sequences_from_contigs_db(self.contigs_db_path, self.target_file_path)
            # Start at half the step size of the output file
            overlap_start = self.step // 2
            first_sequences = self.split_sequences()
            second_sequences = self.split_sequences(overlap_start)

            all_sequences = first_sequences + second_sequences

            # Write combined sequences to output file
            with open(shredded_sequence_file, "w") as output_handle:
                SeqIO.write(all_sequences, output_handle, "fasta")

            # Reading and printing the contents of the file
            with open(shredded_sequence_file, "r") as handle:
                fasta_content = handle.read()
                print("Contents of input_original.fasta:")
                print(fasta_content)

            blast = BLAST(shredded_sequence_file, target_fasta =self.target_file_path, search_program = 'blastn', output_file=blast_output, additional_params = '-dust no')
            blast.evalue = 10 #set Evalue to be same as blastn default
            blast.makedb(dbtype = 'nucl')
            blast.blast(outputfmt = '5', word_size = self.word_size)

        elif self.contigs_db_path and self.profile_db_path:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            #self.splits_basic_info = contigs_db.db.smart_get(t.splits_info_table_name, column = 'split')
            #self.splits_of_interest = contigs_db.db.smart_get(t.splits_info_table_name, column='split')
            self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
            #self.splits_of_interest = list(self.splits_basic_info.keys())
            #self.splits_of_interest = ', '.join(map(str, self.splits_basic_info.keys()))
            #print(self.splits_of_interest)
            #split_path = 'split_of_interest_file.txt

            contigs_db.disconnect()

            #open merged profile-db and get the variable nucleotide table as a dictionary then acess the split names as a list to use in get_snvs
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            self.variable_nucleotides_dict = profile_db.db.get_table_as_dict(t.variable_nts_table_name)


            # Use a list comprehension to extract the values associated with the target key
            split_names = [self.variable_nucleotides_dict[key]['split_name'] for key in self.variable_nucleotides_dict if 'split_name' in self.variable_nucleotides_dict[key]]
            self.split_names_unique = list(dict.fromkeys(split_names))
            sample_id_list = [self.variable_nucleotides_dict[key]['sample_id'] for key in self.variable_nucleotides_dict if 'sample_id' in self.variable_nucleotides_dict[key]]
            sample_id_list = list(set(sample_id_list))
            departure_from_reference = [self.variable_nucleotides_dict[key]['departure_from_reference'] for key in self.variable_nucleotides_dict if 'departure_from_reference' in self.variable_nucleotides_dict[key]]

            profile_db.disconnect()

            #Sort pandas dataframe of SNVs by contig name and then by position of SNV within contig
            self.snv_panda = self.get_snvs().sort_values(by=['contig_name', 'pos_in_contig'])

            self.snv_panda['departure_from_reference'] = self.snv_panda.apply(lambda row: self.variable_nucleotides_dict.get(row.name, {}).get('departure_from_reference', None), axis=1)

            self.all_possible_windows = {} # we will keep this as a dictionary that matches contig name to list of window tuples within that contig
            # structure of self.all_possible_windows: {'contig_0' : [(start0, stop0), (start1, stop1), (start2, stop2), ....... ], 
            #                                           'contig_1': [(start0, stop0), (start1, stop1), (start2, stop2), ....... ], 
            #                                           ....}
            # the windows will not necessarily be sorted within each inner list (yet) because we add windows from one sample at a time

            for split in self.split_names_unique:
                for sample in sample_id_list:
                    split_subset = self.snv_panda.loc[(self.snv_panda.split_name==split)&
                                                      (self.snv_panda.sample_id==sample)&
                                                      (self.snv_panda.departure_from_reference>=self.departure_from_reference_percentage)]
                    if split_subset.shape[0] == 0:
                        continue
                    contig_name = split_subset.contig_name.unique()[0]
                    pos_list = split_subset.pos_in_contig.to_list()

                    if contig_name not in self.all_possible_windows:
                        # If not, initialize it with an empty dictionary
                        self.all_possible_windows[contig_name] = []
                        # subset pandas df with split name

                    #get list of pos within that split
                    for i in range(len(pos_list) - 1):
                        current_pos = pos_list[i]
                        next_pos = pos_list[i + 1]
                        distance = next_pos - current_pos
                        range_start = current_pos
                        range_end = current_pos

                        while i + 1 < len(pos_list) and distance <= self.min_dist_bw_snvs:
                            i += 1
                            current_pos = pos_list[i]
                            if i + 1 < len(pos_list):
                                next_pos = pos_list[i + 1]
                                distance = next_pos - current_pos
                                range_end = current_pos
                        if distance <= self.min_dist_bw_snvs:
                            range_end = next_pos

                        if range_end > range_start:
                            window_start = range_start - self.variable_buffer_length
                            window_end = range_end + self.variable_buffer_length

                        contig_len = len(self.contig_sequences[contig_name]['sequence'])

                        if window_start <0:
                            window_start = 0
                        if window_end > contig_len:
                            window_end = contig_len

                        # Add the window to the contig's list
                        self.all_possible_windows[contig_name].append((window_start, window_end))

            all_merged_snv_windows = {} # this dictionary will be filled up with the merged window list for each contig
            # loop to merge overlaps within a given contig
            for contig_name, window_list in self.all_possible_windows.items():
                # before we check overlaps, we need to sort the list of windows within each contig by the 'start' position (at index 0)
                sorted_windows_in_contig = sorted(window_list, key=lambda x: x[0]) # this list is like the old variable 'all_entries'
                # at this point, sorted_windows_in_contig contains (start, stop) tuples in order of 'start'

                merged_windows_in_contig = []
                while 1:
                    if not len(sorted_windows_in_contig):
                        break

                    entry = sorted_windows_in_contig.pop(0)
                    overlapping_entries = [entry]
                    start, end = entry
                    matching_entries_indices = []

                    for i in range(0, len(sorted_windows_in_contig)):
                        n_start, n_end = sorted_windows_in_contig[i]
                        if self.range_overlapping(start, end, n_start, n_end):
                            matching_entries_indices.append(i)
                            start = min(start, n_start)
                            end = max(end, n_end)

                    # remove each overlapping window from the list and simultaneously add to list of overlapping entries
                    # we do this in backwards order so that pop() doesn't change the indices we need to remove
                    for i in sorted(matching_entries_indices, reverse=True):
                        overlapping_entries.append(sorted_windows_in_contig.pop(i))

                    merged_ranges = self.combine_ranges(overlapping_entries)
                    merged_windows_in_contig.append(merged_ranges)

                all_merged_snv_windows[contig_name] = merged_windows_in_contig

                for i in range(len(merged_windows_in_contig)):
                    for j in range(i, len(merged_windows_in_contig)):
                        if i != j:
                            if self.range_overlapping(merged_windows_in_contig[i][0], merged_windows_in_contig[i][1], merged_windows_in_contig[j][0], merged_windows_in_contig[j][1]):
                                print(f"overlapping at indices {i} and {j} for contig {contig_name}:\n{merged_windows_in_contig[i]}\n{merged_windows_in_contig[j]}")
        #print(all_merged_snv_windows)
        # below is the old version, I've moved it above with minor changes to variable names
            # # Initialize an empty list for unique overlapping sequences
            # all_entries = []

            # for contig_name, window_dict in self.all_possible_windows.items():
            #     for window_name, window_values in window_dict.items():
            #         all_entries.append((contig_name, window_name, window_values['start_position'], window_values['end_position']))

            # # now it is time to identify clusters. the following state
            # clusters = []
            # while 1:
            #     if not len(all_entries):
            #         break

            #     entry = all_entries.pop(0)
            #     cluster = [entry]
            #     contig_name, window_number, start, end = entry
            #     matching_entries = []

            #     for i in range(0, len(all_entries)):
            #         contig_name, n_window_number, n_start, n_end = all_entries[i]
            #         if self.range_overlapping(start, end, n_start, n_end):
            #             matching_entries.append(i)
            #             start = min(start, n_start)
            #             end = max(end, n_end)
            #     #for i in range(0, len(all_entries)):
            #         #contig_name, n_window_number, n_start, n_end = all_entries[i]
            #         #if self.range_overlapping(start, end, n_start, n_end):
            #             #matching_entries.append(i)

            #     # add all matching entries
            #     for i in sorted(matching_entries, reverse=True):
            #         cluster.append(all_entries.pop(i))

            #     # combine ranges of the cluster from entries and then add the combined lsit to the final clusters
            #     combined_result = self.combine_ranges(cluster)
            #     clusters.append(combined_result)
        #print(self.contig_sequences)

        #dumb katy code for turning dict to fasta see utils.func below
        #fasta_file_path = os.path.join(self.temp_dir, "input_original.fasta")

        # Writing to the file
        #with open(fasta_file_path, "w") as handle:
            #for contig_name, contig_data in self.contig_sequences.items():
                #sequence = contig_data.get('sequence', '')  # Get the sequence from the inner dictionary
                #id = contig_name
                #SeqIO.write(SeqRecord(Seq(sequence), id=id), handle, "fasta")

            #export contigs_db to fasta file
            utils.export_sequences_from_contigs_db(self.contigs_db_path, self.target_file_path)

            #get short sequences from all_merged_snv_window and create new fasta from self.target_file_path
            contig_records = []
            for record in SeqIO.parse(self.target_file_path, "fasta"):
                contig_name = record.id
                if contig_name in all_merged_snv_windows:
                    positions = all_merged_snv_windows[contig_name]
                    for i, (start, end) in enumerate(positions, start):
                        section_sequence = record.seq[start:end]
                        section_id = f"{contig_name}_section_{i}_start_bp{start}_end_bp{end}"
                        contig_records.append(SeqRecord(section_sequence, id=section_id, description=""))

            # Write SeqRecord objects to a new FASTA file
            output_fasta_path = os.path.join(self.temp_dir,"output.fasta")
            with open(output_fasta_path, "w") as output_handle:
                SeqIO.write(contig_records, output_handle, "fasta")

            print(f"FASTA file written to {output_fasta_path}")

            # Reading and printing the contents of the file
            with open(output_fasta_path, "r") as handle:
                fasta_content = handle.read()
                print("Contents of input_original.fasta:")
                print(fasta_content)

            blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")

            blast = BLAST(output_fasta_path, target_fasta = self.target_file_path, search_program = 'blastn', output_file=blast_output, additional_params = '-dust no')
            blast.evalue = 10 #set Evalue to be same as blastn default
            blast.makedb(dbtype = 'nucl')
            blast.blast(outputfmt = '5', word_size = self.word_size)

        return(blast_output)

    def split_sequence_at_given_pos(sequence, positions):
        sections = []
        for start, end in positions:
            section_sequence = sequence[start - 1:end]  # Adjust positions for 0-based indexing
            sections.append(section_sequence)
        return sections

    def combine_ranges(self, entries):
        """
        This function takes a list of (contig_name, key, start, end) tuples and takes the longest sequence possible - the smallest start and largest end.

        Returns
        =======
        a tuple containing (contig_name, 'combined', combined_start, combined_end) where the variables are the following:
        contig_name : str
            header of the contig sequence
        combined_start, combined_end : integers
            A new start and end position for a contig sequence, to get the longest possible string. 
        """

        #extract all starts and stops
        all_start = []
        all_end = []
        contig_name = None
        for start, end in entries:
            #contig_name = contig # these should all be the same so it doesn't matter that we overwrite it every iteration of the loop
            all_start.append(start)
            all_end.append(end)
        # do le math
        combined_start = min(all_start)
        combined_end = max(all_end)

        return (combined_start, combined_end)

    def range_overlapping(self, start1, end1, n_start, n_end):
        """
        This function checks if the sections of sequences overlap based on the start and end positions.
        Parameters
        ==========
        start1, end1, n_start, n_end : integer
            Start and end of windows containing SNVs with 20 bp buffer on either side

        Returns
        =======
            :boolean

        """
        return (n_start >= start1 and n_start <= end1) or (n_end >= start1 and n_end <= end1)

    def check_overlap(window1, window2):
                        contig_name_1, start_position_1, end_position_1 = window1[0][1], window1[1][1], window1[2][1]
                        contig_name_2, start_position_2, end_position_2 = window2[0][1], window2[1][1], window2[2][1]

                        return (
                            contig_name_1 == contig_name_2
                            and start_position_1 <= end_position_2
                            and end_position_1 >= start_position_2
                        )

    def get_snvs(self):
        """
        This function takes the contigs.db and the profile.db and finds the SNVs.

        Returns
        =======
        n.data : panada df
            A dataframe with all of the information in the nucleotide variability table of a profile,
            including the splits, parent contigs, positions of the SNVs, the bases of the SNVs and the coverage.
        """
        args = argparse.Namespace(contigs_db=self.contigs_db_path,
                                profile_db=self.profile_db_path,
                                splits_of_interest_set= set(self.split_names_unique),
                                compute_gene_coverage_stats=True)

        n = NucleotidesEngine(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        n.process()

        return n.data

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
        for sequence in SeqIO.parse(self.target_file_path, "fasta"):
            for i in range(start, len(sequence.seq) - self.step + 1, self.step):
                section = sequence.seq[i:i + self.step]
                section_record = SeqRecord(section, id=f"{sequence.id}_part{i//self.step}_start_bp{i}_end_bp{i + self.step}", description="")
                section_sequences.append(section_record)
                if i + self.step > len(sequence.seq):
                    print(sequence.seq)
        return section_sequences

    #def run_blastn(self):
        #"""
        #This function runs the BLASTn search of the split sequences against the original input FASTA to find regions of matching nucleotides.

        #Running the BLASTn generates an xml file of results.

        #Returns
        #blast_output : xml file
            #BLASTn results
        #=======
        #"""
        #blast output file name
        #blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")
        #self.target_file_path = os.path.join(tmp_directory_path,f"input_file.fasta")
        #print(f"cp {self.fasta_file_path} {self.target_file_path}")
        #os.system(f"cp {self.fasta_file_path} {self.target_file_path}")
        #self.run.info('temporary input for blast', self.target_file_path)

        # Start at half the step size of the output file
        #overlap_start = self.step // 2
        #first_sequences = self.split_sequences()
        #second_sequences = self.split_sequences(overlap_start)

        #all_sequences = first_sequences + second_sequences

        # Write combined sequences to output file
        #with open(shredded_sequence_file, "w") as output_handle:
            #SeqIO.write(all_sequences, output_handle, "fasta")
        #need a temporary directory where intermediate files are written, to call on them.

        #blast = BLAST(shredded_sequence_file, target_fasta =self.target_file_path, search_program = 'blastn', output_file=blast_output, additional_params = '-dust no')
        #blast.evalue = 10 #set Evalue to be same as blastn default
        #blast.makedb(dbtype = 'nucl')
        #blast.blast(outputfmt = '5', word_size = self.word_size)

        #blast_command = ["blastn", "-query", output_file, "-subject", self.fasta_file_path, "-out", blast_output,
                         #"-word_size", str(self.word_size), "-dust", "no", "-outfmt", "5"]
        #subprocess.run(blast_command)
        #return blast_output

     #def find_SNV_window(self, profile.db)
        #if SNV:
            #for row()

    def filter_blastn_for_none_identical(self, blast_output):
        """
        This function takes the BLASTn xml output and refines the results to those with less than 100% identity.

        Takes the xml file and filters for hits with less than 100% identity, then gives every hit a name
        with its original position in the sequence, counts the bases that are mismatching and on which strand they occur.
        Finally initialises all of these within a dictionary.

        Parameters
        ==========
        blast_output : xml file
            BLASTn results

        Returns
        mismatch_hits : dict
            A dictionary of all of the BLASTn hits that are less than 100%
        =======

        """
        print(blast_output)
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