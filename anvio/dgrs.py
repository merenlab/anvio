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
        self.num_threads = A('num-threads')
        self.number_of_mismatches = A('number_of_mismatches')
        self.percentage_mismatch = A('percentage_mismatch')
        self.min_mismatching_base_types_vr = A('min_mismatching_base_types_vr')
        self.temp_dir = A('temp_dir') or filesnpaths.get_temp_directory_path()
        self.min_dist_bw_snvs = A('distance_between_snv')
        self.variable_buffer_length = A('variable_buffer_length')
        self.departure_from_reference_percentage = A('departure_from_reference_percentage')
        self.gene_caller_to_consider_in_context = A('gene_caller') or 'prodigal'
        self.min_range_size = A('minimum_range_size')
        self.hmm = A('hmm_usage')
        self.discovery_mode = A('discovery_mode')
        self.output_directory = A('output_dir') or 'DGR-OUTPUT'
        self.parameter_outputs = A('parameter-output')

        # performance
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        self.sanity_check()

        self.run.info('BLASTn word size', self.word_size)
        self.run.info('Skip "N" characters', self.skip_Ns)
        self.run.info('Skip "-" characters', self.skip_dashes)
        self.run.info('Discovery mode', self.discovery_mode)
        self.run.info('Number of Mismatches', self.number_of_mismatches)
        self.run.info('Percentage of Mismatching Bases', self.percentage_mismatch)
        self.run.info('Minimum Mismatching Base Types in VR', self.min_mismatching_base_types_vr)
        self.run.info('Output Directory', self.output_directory)
        self.run.info('Gene Caller Provided', self.gene_caller_to_consider_in_context)
        if self.fasta_file_path:
            self.run.info('Input FASTA file', self.fasta_file_path)
        if self.contigs_db_path:
            self.run.info('Contigs.db', self.contigs_db_path)
        if self.profile_db_path:
            self.run.info('Profile.db', self.profile_db_path)
        if self.fasta_file_path or self.contigs_db_path and not self.profile_db_path:
            self.run.info('Step size', self.step)
        if self.profile_db_path and self.contigs_db_path:
            self.run.info('Minimum distance between SNVs', self.min_dist_bw_snvs)
            self.run.info('Minimum length of SNV window', self.min_range_size)
            self.run.info('Variable buffer length', self.variable_buffer_length)
            self.run.info('Departure from reference percentage', self.departure_from_reference_percentage)
        if self.contigs_db_path:
            self.run.info('HMM(s) Provided', ", ".join(self.hmm))

    def sanity_check(self):
        """Basic checks for a smooth operation"""
        #filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=False)
        #filesnpaths.check_output_directory(self.output_directory)
        #filesnpaths.is_output_dir_writable(self.output_directory)
        if os.path.exists(self.output_directory) and not anvio.FORCE_OVERWRITE:
            raise ConfigError(f"Hold up your directory ({self.output_directory}) exists. To avoid overwriting data we are stopping you here. You have three options "
                            "if you want to continue: rename your directory, delete your existing directory, or rerun with flag `--force-overwrite`.")
        else:
            filesnpaths.gen_output_directory(self.output_directory)

        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                            "class, not multiple :/")

        if ((self.contigs_db_path and self.profile_db_path) or (self.contigs_db_path)) and not self.hmm:
            raise ConfigError("Just so you know in order to report a DGR you need to have Reverse Transcriptases "
                            "annotated or else you are not finding DGRs. If you would like to do this please run "
                            "your 'contigs.db' with `anvi-run-hmms -I Reverse_Transcriptase` (type: 6 clades of DGR "
                            "Retroelements from doi.org/10.1038/s41467-021-23402-7 including other known reverse transcriptases). "
                            "then re-run `anvi-report-dgrs` with the `--hmm-usage` flag. This is one option, but you can use other HMM "
                            "profiles if you like.")

        if self.fasta_file_path:
            # check fasta input
            filesnpaths.is_file_fasta_formatted(self.fasta_file_path)

        if self.step < 0 or self.word_size < 0:
            raise ConfigError('The step value and/or word size value you are trying to input should be positive integer.')

        if self.variable_buffer_length < 0:
            raise ConfigError('The variable buffer length value you are trying to input should be positive integer.')

        if self.departure_from_reference_percentage < 0:
            raise ConfigError('The departure from reference percentage value you are trying to input should be a positive decimal number.')

        if self.contigs_db_path and self.hmm:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            hmm_hits_info_dict = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
            hmm_hits_info_dict = contigs_db.db.smart_get(t.hmm_hits_info_table_name, column = 'source')
            self.hmms_provided = list(hmm_hits_info_dict.keys())

            bad_hmm = []
            for i in self.hmm:
                if i not in self.hmms_provided:
                    bad_hmm.append(i)
                if bad_hmm == ['Reverse_Transcriptase']:
                    raise ConfigError("The classic reverse transcriptase HMM has not been run with these contigs. "
                                        "Please run 'anvi-run-hmms' with the '-I' flag and the 'Reverse_Transcriptase' HMM, "
                                        f"so that anvi'o knows where there are Reverse Transcriptases in your {self.contigs_db_path}.")
                elif len(bad_hmm)>0:
                    if 'Reverse_Transcriptase' in bad_hmm:
                        raise ConfigError("The classic reverse transcriptase HMM has not been run with these contigs. "
                                        "Please run 'anvi-run-hmms' with the '-I' flag and the 'Reverse_Transcriptase' HMM, "
                                        f"so that anvi'o knows where there are Reverse Transcriptases in your {self.contigs_db_path}. "
                                        f"Also you don't have these HMMs you requested: {bad_hmm} in {self.hmms_provided}")
                    else:
                        raise ConfigError(f"You requested these HMMs to be searched through: {bad_hmm} in these HMMs {self.hmms_provided} "
                                        f"that are in your {self.contigs_db_path}. The HMMs you give 'anvi-report-dgrs' need to be in your "
                                        "contigs.db.")

            if self.gene_caller_to_consider_in_context:
                genes_in_contigs_dict = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
                unique_sources = set()
                for gene_id, gene_info in genes_in_contigs_dict.items():
                    unique_sources.add(gene_info['source'])
                unique_sources_list = list(unique_sources)
                if self.gene_caller_to_consider_in_context not in unique_sources_list:
                    raise ConfigError(f"Anvi'o can't find {self.gene_caller_to_consider_in_context} in your {self.contigs_db_path}. "
                                    f"Here are the sources of your genes: {unique_sources_list}.")
            contigs_db.disconnect()

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
        self.run.info('Temporary (contig) reference input for blast', self.target_file_path)

        if self.fasta_file_path or (self.contigs_db_path and not self.profile_db_path):
            shredded_sequence_file = os.path.join(tmp_directory_path,f"shredded_sequences_step_{self.step}_wordsize_{self.word_size}.fasta")
            self.blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")
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

            blast = BLAST(shredded_sequence_file, target_fasta = self.target_file_path, search_program = 'blastn', output_file = self.blast_output, additional_params = '-dust no')
            blast.evalue = 10 #set Evalue to be same as blastn default
            blast.makedb(dbtype = 'nucl')
            blast.blast(outputfmt = '5', word_size = self.word_size)

        elif self.contigs_db_path and self.profile_db_path:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)

            contigs_db.disconnect()

            #open merged profile-db and get the variable nucleotide table as a dictionary then acess the split names as a list to use in get_snvs
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            #Sort pandas dataframe of SNVs by contig name and then by position of SNV within contig
            self.snv_panda = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
            self.snv_panda['contig_name'] = self.snv_panda.split_name.str.split('_split_').str[0]
            self.split_names_unique = utils.get_all_item_names_from_the_database(self.profile_db_path)

            profile_db.disconnect()

            sample_id_list = list(set(self.snv_panda.sample_id.unique()))

            self.all_possible_windows = {} # we will keep this as a dictionary that matches contig name to list of window tuples within that contig
            # structure of self.all_possible_windows: {'contig_0' : [(start0, stop0), (start1, stop1), (start2, stop2), ....... ],
            #                                           'contig_1': [(start0, stop0), (start1, stop1), (start2, stop2), ....... ],
            #                                           ....}
            # the windows will not necessarily be sorted within each inner list (yet) because we add windows from one sample at a time

            # filter snv_panda for min departure from ref and codon position
            if self.discovery_mode:
                self.snv_panda = self.snv_panda.loc[self.snv_panda.departure_from_reference>=self.departure_from_reference_percentage]
            else:
                self.snv_panda = self.snv_panda.loc[(self.snv_panda.departure_from_reference>=self.departure_from_reference_percentage) &
                                                    ((self.snv_panda.base_pos_in_codon == 1) | (self.snv_panda.base_pos_in_codon == 2))]

            for split in self.split_names_unique:
                split_subset = self.snv_panda.loc[self.snv_panda.split_name==split]
                for sample in sample_id_list:
                    split_subset_sample = split_subset.loc[split_subset.sample_id==sample]
                    if split_subset_sample.shape[0] == 0:
                        continue
                    contig_name = split_subset_sample.contig_name.unique()[0]
                    pos_list = split_subset_sample.pos_in_contig.to_list()

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

                        if (range_end - range_start) < self.min_range_size:
                            continue
                        else:
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

            #export contigs_db to fasta file
            utils.export_sequences_from_contigs_db(self.contigs_db_path, self.target_file_path)

            #get short sequences from all_merged_snv_window and create new fasta from self.target_file_path
            contig_records = []
            for record in SeqIO.parse(self.target_file_path, "fasta"):
                contig_name = record.id
                if contig_name in all_merged_snv_windows:
                    self.positions= all_merged_snv_windows[contig_name]
                    for i, (start, end) in enumerate(self.positions, start):
                        section_sequence = record.seq[start:end]
                        section_id = f"{contig_name}_section_{i}_start_bp{start}_end_bp{end}"
                        contig_records.append(SeqRecord(section_sequence, id=section_id, description=""))

            # Write SeqRecord objects to a new FASTA file
            output_fasta_path = os.path.join(self.temp_dir,"output.fasta")
            with open(output_fasta_path, "w") as output_handle:
                SeqIO.write(contig_records, output_handle, "fasta")

            self.run.info('Temporary (SNV window) query input for blast', output_fasta_path)

            self.blast_output = os.path.join(tmp_directory_path,f"blast_output_step_{self.step}_wordsize_{self.word_size}.xml")
            blast = BLAST(output_fasta_path, target_fasta = self.target_file_path, search_program = 'blastn',
                        output_file=self.blast_output, additional_params = '-dust no', num_threads=self.num_threads)
            blast.evalue = 10 #set Evalue to be same as blastn default
            blast.makedb(dbtype = 'nucl')
            blast.blast(outputfmt = '5', word_size = self.word_size)

        return

    def split_sequence_at_given_pos(self, sequence):
        sections = []
        for start, end in self.positions:
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
            all_start.append(start)
            all_end.append(end)
        # do le maths
        combined_start = min(all_start)
        combined_end = max(all_end)

        return (combined_start, combined_end)

    def range_overlapping(self, start1, end1, n_start, n_end):
        """
        This function checks if the sections of sequences overlap based on the start and end positions.
        Parameters
        ==========
        start1, end1, n_start, n_end : integer
            Start and end of windows containing SNVs with 20 bp buffer on either side.

        Returns
        =======
        A boolean indicating whether the ranges overlap.
        """
        return (n_start >= start1 and n_start <= end1) or (n_end >= start1 and n_end <= end1) or (start1 >= n_start and start1 <= n_end and end1 >= n_start and end1 <= n_end)

    def check_overlap(window1, window2):
        """
        WRITE ME PLEASE!!!
        This function checks if the sections of sequences overlap based on the start and end positions.
        Parameters
        ==========
        start1, end1, n_start, n_end : integer
            Start and end of windows containing SNVs with 20 bp buffer on either side.

        Returns
        =======
        A boolean indicating whether the ranges overlap.
        """
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
        This function is not currently in use but here for Merens sake.
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

        n = NucleotidesEngine(args, r=terminal.Run(verbose=True), p=terminal.Progress(verbose=True))
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

    def filter_blastn_for_none_identical(self):
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
        tree = ET.parse(self.blast_output)
        root = tree.getroot()

        max_percent_identity = 100
        self.mismatch_hits = {}

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

                        self.mismatch_hits[hit_identity] = {}

                        qseq = str(hsp.find('Hsp_qseq').text)
                        hseq = str(hsp.find('Hsp_hseq').text)
                        midline = str(hsp.find('Hsp_midline').text)
                        subject_genome_start_position = min([int(hsp.find('Hsp_hit-from').text), int(hsp.find('Hsp_hit-to').text)])
                        subject_genome_end_position = max([int(hsp.find('Hsp_hit-from').text), int(hsp.find('Hsp_hit-to').text)])
                        alignment_length = int(hsp.find('Hsp_align-len').text)
                        query_genome_start_position = query_start_position + min([int(hsp.find('Hsp_query-from').text), int(hsp.find('Hsp_query-to').text)])
                        query_genome_end_position = query_start_position + max([int(hsp.find('Hsp_query-from').text), int(hsp.find('Hsp_query-to').text)])
                        query_frame = int(hsp.find('Hsp_query-frame').text)
                        subject_frame = int(hsp.find('Hsp_hit-frame').text)

                        query_mismatch_positions = []

                        # Unique characters that may appear in qseq and hseq
                        all_possible_characters = set(qseq + hseq)

                        # Initialize counts with all possible characters
                        query_mismatch_counts = {char: 0 for char in all_possible_characters}
                        subject_mismatch_counts = {char: 0 for char in all_possible_characters}

                        chars_to_skip = [self.skip_dashes, self.skip_Ns]

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

                        # get contigs name
                        query_contig = section_id.split('_section', 1)[0]
                        subject_contig = hit.find('Hit_def').text

                        self.mismatch_hits[hit_identity] = {
                            'query_seq': qseq,
                            'hit_seq': hseq,
                            'midline': midline,
                            'query_contig': query_contig,
                            'subject_contig': subject_contig,
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
        return self.mismatch_hits

    def add_new_DGR(self, DGR_number, TR_is_query, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig, base,
                    is_reverse_complement, VR_sequence, subject_genome_start_position, subject_genome_end_position, subject_contig, midline,
                    percentage_of_mismatches):
        """
        This function is adding the DGRs that are found initially to the DGRs_found_dict.

        This is to remove redundancy from the previous function when everything was in 'filter_for_TR_VR', now it should be easier to read and edit this code.

        Parameters
        ==========
        (a lot)
        all of the input arguments from the BLASTn output : dict keys (different types listed below)
            DGR_number : integer argument
            TR_is_query, is_reverse_complement : boolean
            TR_sequence, query_contig, base, VR_sequence, subject_contig : string
            query_genome_start_position, query_genome_end_position, subject_genome_start_position, subject_genome_end_position : integer (bp position)
            midline : charaacter with defined spacing
            percentage_of_mismatches : float
                Keys of a dictionary containing all of the BLASTn hits that are less than 100%

        Creates
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions and the corrosponding info for those regions

        """

        DGR_key = f'DGR_{DGR_number:03d}'
        #print(f"adding new DGR, {DGR_key}. TR_is_query: {TR_is_query}")
        self.DGRs_found_dict[DGR_key] = {}
        # TR stuff
        self.DGRs_found_dict[DGR_key]['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['base'] = base
        self.DGRs_found_dict[DGR_key]['TR_reverse_complement'] = is_reverse_complement
        # VR stuff
        self.DGRs_found_dict[DGR_key]['VRs'] = {'VR_001':{}}
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_sequence'] = VR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['midline'] = midline
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['percentage_of_mismatches'] = percentage_of_mismatches

        # query-subject dependent stuff
        if TR_is_query:
            self.DGRs_found_dict[DGR_key]['TR_start_position'] = query_genome_start_position
            self.DGRs_found_dict[DGR_key]['TR_end_position'] = query_genome_end_position
            self.DGRs_found_dict[DGR_key]['TR_contig'] = query_contig
            self.DGRs_found_dict[DGR_key]['TR_sequence_found'] = 'query'

            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_start_position'] = subject_genome_start_position
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_end_position'] = subject_genome_end_position
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_contig'] = subject_contig
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_sequence_found'] = 'subject'

            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_start_position'] = query_genome_start_position
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_end_position'] = query_genome_end_position

        else:
            self.DGRs_found_dict[DGR_key]['TR_start_position'] = subject_genome_start_position
            self.DGRs_found_dict[DGR_key]['TR_end_position'] = subject_genome_end_position
            self.DGRs_found_dict[DGR_key]['TR_contig'] = subject_contig
            self.DGRs_found_dict[DGR_key]['TR_sequence_found'] = 'subject'

            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_start_position'] = query_genome_start_position
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_end_position'] =   query_genome_end_position
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_contig'] = query_contig
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_sequence_found'] = 'query'

            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_start_position'] = subject_genome_start_position
            self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_end_position'] = subject_genome_end_position

    def update_existing_DGR(self, existing_DGR_key, VR_sequence, TR_sequence, midline, percentage_of_mismatches, query_genome_start_position,
                            query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position, subject_contig):
        """
        This function is updating the DGRs in the DGRs_found_dict with those DGRs that overlap each other.

        This is to remove redundancy from the previous function when everything was in 'filter_for_TR_VR', now it should be easier to read and edit this code.

        Parameters
        ==========
        (a lot)
        all of the input arguments from the BLASTn output : dict keys (different types listed below)
            existing_DGR_key : integer argument
            TR_is_query, is_reverse_complement : boolean
            TR_sequence, query_contig, base, VR_sequence, subject_contig : string
            query_genome_start_position, query_genome_end_position, subject_genome_start_position, subject_genome_end_position : integer (bp position)
            midline : charaacter with defined spacing
            percentage_of_mismatches : float
                Keys of a dictionary containing all of the BLASTn hits that are less than 100%

        Creates
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions and the corrosponding info for those regions

        """
        # VR info
        num_VR = len(self.DGRs_found_dict[existing_DGR_key]['VRs']) + 1
        new_VR_key = f'VR_{num_VR:03d}'
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key] = {}
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_sequence'] = VR_sequence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['midline'] = midline
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['percentage_of_mismatches'] = percentage_of_mismatches

        if self.DGRs_found_dict[existing_DGR_key]['TR_sequence_found'] == 'query':
            # update TR start and end
            self.DGRs_found_dict[existing_DGR_key]['TR_start_position'] = min(query_genome_start_position, self.DGRs_found_dict[existing_DGR_key]['TR_start_position'])
            self.DGRs_found_dict[existing_DGR_key]['TR_end_position'] = max(query_genome_end_position, self.DGRs_found_dict[existing_DGR_key]['TR_end_position'])

            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_start_position'] = subject_genome_start_position
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_end_position'] = subject_genome_end_position
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_contig'] = subject_contig
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_sequence_found'] = 'subject'
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_start_position'] = query_genome_start_position
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_end_position'] = query_genome_end_position
        else:
            self.DGRs_found_dict[existing_DGR_key]['TR_start_position'] = min(subject_genome_start_position, self.DGRs_found_dict[existing_DGR_key]['TR_start_position'])
            self.DGRs_found_dict[existing_DGR_key]['TR_end_position'] = max(subject_genome_end_position, self.DGRs_found_dict[existing_DGR_key]['TR_end_position'])

            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_start_position'] = query_genome_start_position
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_end_position'] = query_genome_end_position
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_contig'] = query_contig
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_sequence_found'] = 'query'
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_start_position'] = subject_genome_start_position
            self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_end_position'] = subject_genome_end_position

    def filter_for_TR_VR(self):
        #NOTE: need to check this code to clean up and maybe put into one function to remove redundancy - Iva offered to help :)
        # We did this on 21.05.2024, created functions add_new_DGR and update_existing_DGR
        """
        This function takes the none identical hits of the BLASTn and filters for template and variable regions.

        This works by filtering for sequences that have an overrepresentation of one base that is mismatching and a certain number
        one type of base mismatching within the sequence, defined by the number of mismatches argument.

        Parameters
        ==========
        mismatch_hits : dict
            A dictionary of all of the BLASTn hits that are less than 100%

        Creates
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        """
        num_DGR = 0

        #possible DGR dictionary
        self.DGRs_found_dict = {}

        for sequence_component, hit_data in self.mismatch_hits.items():
            query_mismatch_counts = hit_data['query_mismatch_counts']
            subject_mismatch_counts = hit_data['subject_mismatch_counts']
            position = hit_data['position']
            subject_genome_start_position = hit_data['subject_genome_start_position']
            subject_genome_end_position = hit_data['subject_genome_end_position']
            alignment_length = hit_data['alignment_length']
            subject_sequence = Seq(hit_data['hit_seq'])
            original_midline = hit_data['midline']
            query_sequence = Seq(hit_data['query_seq'])
            shredded_sequence_name = sequence_component
            query_genome_start_position = hit_data['query_genome_start_position']
            query_genome_end_position = hit_data['query_genome_end_position']
            query_frame = hit_data['query_frame']
            subject_frame = hit_data['subject_frame']
            query_contig = hit_data['query_contig']
            subject_contig = hit_data['subject_contig']

            # get number of mismatches
            mismatch_length_bp = len(position)

            # if num of mismatches = 0, skip DGR search sanity check
            if mismatch_length_bp == 0:
                continue
            else:
                is_VR = False
                for letter, count in subject_mismatch_counts.items():
                    non_zero_bases = sum(1 for count in subject_mismatch_counts.values() if count > 0)
                    if non_zero_bases >= self.min_mismatching_base_types_vr:
                        is_VR = True

                # Calculate the percentage identity of each alignment
                is_TR = False
                for letter, count in query_mismatch_counts.items():
                    percentage_of_mismatches = (count / mismatch_length_bp)
                    if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches):
                        is_TR = True
                        # if the letter is T, then we assume that it is an A base and we reverse completment EVERYTHING
                        if letter == 'T':
                            TR_sequence = str(query_sequence.reverse_complement())
                            VR_sequence = str(subject_sequence.reverse_complement())
                            midline = ''.join(reversed(original_midline))
                            base = 'A'
                            is_reverse_complement = True
                        else:
                            TR_sequence = str(query_sequence)
                            VR_sequence = str(subject_sequence)
                            midline = original_midline
                            base = letter
                            is_reverse_complement = False

                        # Only add the DGR if both is_TR and is_VR are True
                        if is_TR and is_VR:
                            #need to check if the new TR and VR you're looping through exsists in the DGR_found_dict, see if position overlap
                            if not self.DGRs_found_dict:
                                # add first DGR
                                num_DGR += 1
                                self.add_new_DGR(num_DGR, True, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                            base, is_reverse_complement, VR_sequence, subject_genome_start_position, subject_genome_end_position,
                                            subject_contig, midline, percentage_of_mismatches)

                            else:
                                was_added = False
                                for dgr in self.DGRs_found_dict:
                                    if self.DGRs_found_dict[dgr]['TR_contig'] == query_contig and self.range_overlapping(query_genome_start_position,
                                                                                                                    query_genome_end_position,
                                                                                                                    self.DGRs_found_dict[dgr]['TR_start_position'],
                                                                                                                    self.DGRs_found_dict[dgr]['TR_end_position']):
                                        was_added = True
                                        self.update_existing_DGR(dgr, VR_sequence, TR_sequence, midline, percentage_of_mismatches, query_genome_start_position,
                                                        query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                                                        subject_contig)
                                        break
                                if not was_added:
                                    # add new TR and its first VR
                                    num_DGR += 1
                                    self.add_new_DGR(num_DGR, True, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                            base, is_reverse_complement, VR_sequence, subject_genome_start_position, subject_genome_end_position,
                                            subject_contig, midline, percentage_of_mismatches)

                if not is_TR:
                    is_VR = False
                    for letter, count in query_mismatch_counts.items():
                        non_zero_bases = sum(1 for count in query_mismatch_counts.values() if count > 0)
                        if non_zero_bases >= self.min_mismatching_base_types_vr:
                            is_VR = True
                    # Calculate the percentage identity of each alignment
                    for letter, count in subject_mismatch_counts.items():
                            percentage_of_mismatches = (count / mismatch_length_bp)
                            if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches):
                                is_TR = True
                                # if the letter is T, then we assume that it is an A base and we reverse completment EVERYTHING
                                if letter == 'T':
                                    TR_sequence = str(subject_sequence.reverse_complement())
                                    VR_sequence = str(query_sequence.reverse_complement())
                                    midline = ''.join(reversed(original_midline))
                                    base = 'A'
                                    is_reverse_complement = True
                                else:
                                    TR_sequence = str(subject_sequence)
                                    VR_sequence = str(query_sequence)
                                    midline = original_midline
                                    base = letter
                                    is_reverse_complement = False

                                # Only add the DGR if both is_TR and is_VR are True
                                if is_TR and is_VR:
                                    #need to check if the new TR youre looping through exsists in the DGR_found_dict, see if position overlap
                                    if not self.DGRs_found_dict:
                                        # add first DGR
                                        num_DGR += 1
                                        self.add_new_DGR(num_DGR, False, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                            base, is_reverse_complement, VR_sequence, subject_genome_start_position, subject_genome_end_position,
                                            subject_contig, midline, percentage_of_mismatches)

                                    else:
                                        was_added = False
                                        for dgr in self.DGRs_found_dict:
                                            if self.DGRs_found_dict[dgr]['TR_contig'] == subject_contig and self.range_overlapping(subject_genome_start_position,
                                                                                                                            subject_genome_end_position,
                                                                                                                            self.DGRs_found_dict[dgr]['TR_start_position'],
                                                                                                                            self.DGRs_found_dict[dgr]['TR_end_position']):
                                                was_added = True
                                                #TODO can rename concensus_TR
                                                self.update_existing_DGR(dgr, VR_sequence, TR_sequence, midline, percentage_of_mismatches, query_genome_start_position,
                                                        query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                                                        subject_contig)
                                                break

                                        if not was_added:
                                            # add new TR and its first VR
                                            num_DGR += 1
                                            self.add_new_DGR(num_DGR, False, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                            base, is_reverse_complement, VR_sequence, subject_genome_start_position, subject_genome_end_position,
                                            subject_contig, midline, percentage_of_mismatches)

        if anvio.DEBUG:
            self.run.warning(f"The temp directory, '{self.temp_dir}', is kept. Don't forget to clean it up later!", header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (use `--debug` to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(self.temp_dir)

        return


    def get_gene_info(self):
        """
        This function collects information genes in which the variable regions act in and prints all the information to a CSV file.

        Parameters
        ==========

        Creates
        =======
        DGR_genes_found.csv : csv
            A cav file containing the variable regions and the genes they act in, inclduing basic information.

        """
        # initiate a dictionary for the gene where we find VR
        self.vr_gene_info = {}
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)

        # are there genes?
        if not contigs_db.meta['genes_are_called']:
            self.run.warning("There are no gene calls in your contigs database, therefore there is context to "
                            "learn about :/ Your reports will not include a file to study the genomic context "
                            "that surrounds the DGRs")

            contigs_db.disconnect()

        # are there functions?
        function_sources_found = contigs_db.meta['gene_function_sources'] or []
        if not len(function_sources_found):
            self.run.warning("There are no functions for genes in your contigs database :/ Your reports on DGRs "
                            "will not include the function of genes with VRs. SHAME.")

        # now we will go through each VRs to populate `self.vr_gene_info`
        # with gene calls and functions
        gene_calls_per_contig = {}

        for dgr, tr in self.DGRs_found_dict.items():
            for vr, vr_data in tr['VRs'].items():
                contig_name = vr_data['VR_contig']
                vr_start = vr_data['VR_start_position']
                vr_end = vr_data['VR_end_position']
                if dgr not in self.vr_gene_info:
                    self.vr_gene_info[dgr] = {vr:{}}
                else:
                    self.vr_gene_info[dgr][vr] = {}

                # get any genes overlapping the start and end position of the VR
                if contig_name not in gene_calls_per_contig:
                    where_clause = f'''contig="{contig_name}" and source="{self.gene_caller_to_consider_in_context}"'''
                    gene_calls_per_contig[contig_name] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

                gene_calls_in_contig = gene_calls_per_contig[contig_name]

                # find gene caller overlapping with VR
                for gene_callers_id, gene_call in gene_calls_in_contig.items():
                    if gene_call['start'] <= vr_start and gene_call['stop'] >= vr_end:
                        gene_call['gene_callers_id'] = gene_callers_id

                        # if there are funtion sources, let's recover them for our gene of interest
                        if function_sources_found:
                            where_clause = f'''gene_callers_id="{gene_callers_id}"'''
                            hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
                            if hits:
                                gene_call['functions'] = [h['function'] for h in hits]
                                gene_call['sources'] = [h['source'] for h in hits]
                                gene_call['accessions'] = [h['accession'] for h in hits]
                            else:
                                gene_call['functions'] = []
                                gene_call['sources'] = []
                                gene_call['accessions'] = []
                        else:
                            gene_call['functions'] = []
                            gene_call['sources'] = []
                            gene_call['accessions'] = []

                        # While we are here, let's add more info about the gene
                        # DNA sequence:
                        dna_sequence = self.contig_sequences[contig_name]['sequence'][gene_call['start']:gene_call['stop']]
                        rev_compd = None
                        if gene_call['direction'] == 'f':
                            gene_call['DNA_sequence'] = dna_sequence
                            rev_compd = False
                        else:
                            gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)
                            rev_compd = True

                        # add AA sequence
                        where_clause = f'''gene_callers_id == "{gene_callers_id}"'''
                        aa_sequence = contigs_db.db.get_some_rows_from_table_as_dict(t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                        gene_call['AA_sequence'] = aa_sequence[gene_callers_id]['sequence']

                        # gene length
                        gene_call['length'] = gene_call['stop'] - gene_call['start']

                        # add fasta header
                        header = '|'.join([f"contig:{contig_name}",
                                    f"start:{gene_call['start']}",
                                    f"stop:{gene_call['stop']}",
                                    f"direction:{gene_call['direction']}",
                                    f"rev_compd:{rev_compd}",
                                    f"length:{gene_call['length']}"])
                        #gene_call['header'] = ' '.join([str(gene_callers_id), header])

                        self.vr_gene_info[dgr][vr] = gene_call
                        break

        #MAKE into WRITE DGR_genes_found write function
        #define output path
        output_directory_path = self.output_directory
        output_path_for_genes_found = os.path.join(output_directory_path, "DGR_genes_found.csv")

        # Define the header for the CSV file
        csv_header = ['DGR_ID', 'VR_ID', 'Contig', 'Start', 'Stop', 'Direction', 'Partial', 'Call_Type', 'Gene_Caller_Source', 'Version', 'Gene_Caller_ID', 'DNA_Sequence', 'AA_Sequence', 'Length', 'Gene_Functions', 'Gene_Function_Source', 'Gene_Function_Accession']

        # Open the CSV file in write mode
        with open(output_path_for_genes_found, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(csv_header)  # Write the header row
            # Iterate through the dictionary and write each gene's information to the CSV file
            for dgr_id, vr_data in self.vr_gene_info.items():
                for vr_id, gene_info in vr_data.items():
                    if not gene_info:
                        continue

                    # Convert list of functions to string
                    gene_functions = ''.join(gene_info['functions'])
                    gene_annotation_source = ''.join(gene_info['sources'])
                    gene_annotation_accession = ''.join(gene_info['accessions'])

                    writer.writerow([
                        dgr_id,
                        vr_id,
                        gene_info['contig'],
                        gene_info['start'],
                        gene_info['stop'],
                        gene_info['direction'],
                        gene_info['partial'],
                        gene_info['call_type'],
                        gene_info['source'],
                        gene_info['version'],
                        gene_info['gene_callers_id'],
                        gene_info['DNA_sequence'],
                        gene_info['AA_sequence'],
                        gene_info['length'],
                        gene_functions,
                        gene_annotation_source,
                        gene_annotation_accession
                    ])
        return


    def get_hmm_info(self):
        """
        This function creates a dictionary of the HMMs provided that are the closest to the template
        regions, this is done by finding the middle of the template region and the middle of the HMM gene
        and comparing those pairs to each other to find the shortest distance.
        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======
        : csv
            A csv tabular file containing the template and variable regions

        """

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.hmm_hits_in_splits_dict = contigs_db.db.get_table_as_dict(t.hmm_hits_splits_table_name)
        genes_in_contigs = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        self.hmm_hits_dict = contigs_db.db.get_table_as_dict(t.hmm_hits_table_name)

        found_HMMS_dict = {}
        # go to hmm_hits, for every unique gene_caller_ids pick lowest e-value. now have gene callers ID
        for index, entry in self.hmm_hits_dict.items():
            if entry['source'] in self.hmm:
                gene_callers_id = str(entry['gene_callers_id'])
                gene_name = entry['gene_name']
                e_value = entry['e_value']
                HMM_source = entry['source']
                if not gene_callers_id in found_HMMS_dict.keys():
                    found_HMMS_dict[gene_callers_id] = {'gene_name':gene_name,
                                                            'e_value':e_value,
                                                            'HMM_source':HMM_source}
                elif e_value < found_HMMS_dict[gene_callers_id]['e_value']:
                    found_HMMS_dict[gene_callers_id]['gene_name'] = gene_name
                    found_HMMS_dict[gene_callers_id]['e_value'] = e_value
                    found_HMMS_dict[gene_callers_id]['HMM_source'] = HMM_source

        # Check if the gene_caller_id exists in genes_in_contigs
        for gene_callers_id, hmm_dict in found_HMMS_dict.items():
            # Retrieve the 'start' and 'stop' values
            int_gene_callers_id = int(gene_callers_id)
            start_value = genes_in_contigs[int_gene_callers_id]['start']
            stop_value = genes_in_contigs[int_gene_callers_id]['stop']
            contig = genes_in_contigs[int_gene_callers_id]['contig']
            gene_annotation_source = genes_in_contigs[int_gene_callers_id]['source']
            direction = genes_in_contigs[int_gene_callers_id]['direction']

            # add new info about each hmm
            hmm_dict['HMM_start'] = start_value
            hmm_dict['HMM_stop'] = stop_value
            hmm_dict['HMM_Midpoint'] = int((start_value + stop_value)/2)
            hmm_dict['contig'] = contig
            hmm_dict['HMM_direction'] = direction
            hmm_dict['Gene_annotation_source'] = gene_annotation_source

            found_HMMS_dict[gene_callers_id] = hmm_dict

        #look at general concensus TR in the level up so all the TRs have the same HMM if in the same DGR.
        for DGR_id, DGR_info in self.DGRs_found_dict.items():
            TR_start_position = DGR_info['TR_start_position']
            TR_end_position = DGR_info['TR_end_position']
            TR_middle_position = (TR_start_position + TR_end_position) / 2

            # Initialize closest_distances dictionary inside the loop
            closest_distances = {}  # Initialize an empty dictionary to store closest distances
            HMM_found = False

            for gene_callers_id, hmm_dict in found_HMMS_dict.items():
                if DGR_info['TR_contig'] == hmm_dict['contig']:
                    HMM_found = True
                    HMM_midpoint = hmm_dict['HMM_Midpoint']
                    distance = abs(TR_middle_position - HMM_midpoint)

                    if not closest_distances:
                        closest_distances = {'gene_callers_id': gene_callers_id, 'distance': distance}
                    elif distance < closest_distances['distance']:
                        closest_distances['gene_callers_id'] = gene_callers_id
                        closest_distances['distance'] = distance
            if not HMM_found:
                HMM_gene_callers_id = ''
                DGR_info['HMM_gene_callers_id'] = ''
                DGR_info['distance_to_HMM'] = ''
                DGR_info['HMM_start'] = ''
                DGR_info['HMM_stop'] = ''
                DGR_info['HMM_direction'] = ''
                DGR_info['HMM_source'] = ''
                DGR_info['HMM_gene_name'] = ''
                DGR_info['HMM_gene_source'] = ''
            else:
                HMM_gene_callers_id = closest_distances['gene_callers_id']
                DGR_info['HMM_gene_callers_id'] = HMM_gene_callers_id
                DGR_info['distance_to_HMM'] = closest_distances['distance']
                DGR_info['HMM_start'] = found_HMMS_dict[HMM_gene_callers_id]['HMM_start']
                DGR_info['HMM_stop'] = found_HMMS_dict[HMM_gene_callers_id]['HMM_stop']
                DGR_info['HMM_direction'] = found_HMMS_dict[HMM_gene_callers_id]['HMM_direction']
                DGR_info['HMM_source'] = found_HMMS_dict[HMM_gene_callers_id]['HMM_source']
                DGR_info['HMM_gene_name'] = found_HMMS_dict[HMM_gene_callers_id]['gene_name']
                DGR_info['HMM_gene_source'] = found_HMMS_dict[HMM_gene_callers_id]['Gene_annotation_source']

        return

    def create_found_tr_vr_csv(self):
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
        if not self.DGRs_found_dict:
            raise ConfigError("No DGRS were found so no output file will be written")

        output_directory_path = self.output_directory
        output_path_dgrs = os.path.join(output_directory_path, "DGRs_found.csv")
        with open(output_path_dgrs, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)

            # Write header
            csv_writer.writerow(["DGR", "VR", "VR_contig", "VR_sequence", "Midline", "VR_start_position", "VR_end_position", "Mismatch %",
                                "TR_contig", "TR_sequence", "Base", "Reverse Complement", "TR_start_position", "TR_end_position","HMM_source",
                                "distance_to_HMM", "HMM_gene_name", "HMM_direction", "HMM_start", "HMM_stop",  "HMM_gene_callers_id"])

            for dgr, tr in self.DGRs_found_dict.items():
                for vr, vr_data in tr['VRs'].items():
                    # Access the closest_HMM dictionary within each VR
                    closest_hmm_dict = tr.get('closest_HMM', {})

                    # Write VR data to CSV
                    csv_row = [dgr, vr, vr_data['VR_contig'], vr_data['VR_sequence'], vr_data['midline'],
                            vr_data['VR_start_position'], vr_data['VR_end_position'], vr_data['percentage_of_mismatches'],
                            tr['TR_contig'], vr_data['TR_sequence'], tr['base'], tr['TR_reverse_complement'],
                            vr_data['TR_start_position'], vr_data['TR_end_position'], tr['HMM_source'], tr["distance_to_HMM"],
                            tr["HMM_gene_name"], tr["HMM_direction"], tr["HMM_start"], tr["HMM_stop"], tr["HMM_gene_callers_id"]]

                    # Write the CSV row
                    csv_writer.writerow(csv_row)
            return

    def parameter_output_sheet(self):
        """
        This function creates a csv tabular format of all the parameters the user input in the current run.

        Returns
        =======
        : csv
            A csv tabular file containing the template and variable regions

        """
        output_path_parameters = os.path.join(self.output_directory, "Parameters_used_in_DGRs_found.csv")
        with open(output_path_parameters, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            headers = ["Parameter", "Value"]
            csv_writer.writerow(headers)

            parameters = [
                ("Contig.db", self.contigs_db_path if self.contigs_db_path else None),
                ("Profile.db", self.profile_db_path if self.profile_db_path else None),
                ("Word Size of BLAST", self.word_size if self.word_size else "8"),
                ("Number of Threads for BLAST", self.num_threads if self.num_threads else None),
                ("Skip 'Ns'", self.skip_Ns if self.skip_Ns else "FALSE"),
                ("Skip '-'", self.skip_dashes if self.skip_dashes else "FLASE"),
                ("Number of Mismatches", self.number_of_mismatches if self.number_of_mismatches else "7"),
                ("Percentage of Mismatches", self.percentage_mismatch if self.percentage_mismatch else "0.8"),
                ("Minimum Mismatching Base Types in VR", self.min_mismatching_base_types_vr if self.min_mismatching_base_types_vr else "3"),
                ("Temporary Directory", self.temp_dir if self.temp_dir else None),
                ("Distance between SNVs", self.min_dist_bw_snvs if self.min_dist_bw_snvs else "5"),
                ("Variable Buffer Length", self.variable_buffer_length if self.variable_buffer_length else "20"),
                ("Departure from Reference (Percentage)", self.departure_from_reference_percentage if self.departure_from_reference_percentage else "0.1"),
                ("Minimum Range size of High Density SNVs", self.min_range_size if self.min_range_size else "5"),
                ("Gene caller", self.gene_caller_to_consider_in_context if self.gene_caller_to_consider_in_context else "prodigal"),
                ("HMMs Provided to Search through", self.hmm if self.hmm else "Reverse_Transcriptase"),
                ("Discovery mode", self.discovery_mode if self.discovery_mode else "FALSE"),
                ("Output Directoy", self.output_directory if self.output_directory else "default")
            ]

            csv_writer.writerows(parameters)

        return

    def process(self,args):
        """Here we process all of the functions in our class and call upon different functions depending on the inputs used"""
        self.sanity_check()
        self.get_blast_results()
        self.filter_blastn_for_none_identical()
        self.filter_for_TR_VR()
        if args.parameter_output:
            self.run.info_single("Writing to Parameters used file.")
            self.run.info_single('\n')
            self.parameter_output_sheet()

        if self.fasta_file_path:
            return

        else:
            self.run.info_single("Computing the Genes the Variable Regions occur in and creating a 'DGR_genes_found.csv'.")
            self.run.info_single('\n')
            self.get_gene_info()
            self.run.info_single("Computing the closest HMMs to the Template Regions and printing them in your output csv.")
            self.get_hmm_info()
            self.create_found_tr_vr_csv()

        return