# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity Generating Retroelements"""

import re
import xml.etree.ElementTree as ET
import csv
import os
import shutil
import argparse
import copy
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths
import anvio.tables as t
import multiprocess as multiprocessing

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.variabilityops import NucleotidesEngine
from anvio.summaryhtml import SummaryHTMLOutput
from anvio.sequencefeatures import PrimerSearch
from anvio.constants import nucleotides

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Katy Lambert-Slosarska"
__email__ = "klambertslosarska@gmail.com"
__status__ = "Development"

PL = terminal.pluralize
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
        self.parameter_outputs = A('parameter_output')
        self.just_do_it = A('just_do_it')
        self.metagenomics_contigs_mode = A('metagenomics_contigs_mode')
        self.collections_given = A('collection_name')
        self.skip_recovering_genomic_context = A('skip_recovering_genomic_context')
        self.num_genes_to_consider_in_context = A('num_genes_to_consider_in_context') or 3
        self.samples_txt = A('samples_txt')
        self.variable_region_primer_base_length = A('variable_region_primer_base_length') or 12 #TODO test different values for this. If Illumina reads are 250 bases then depends on length of VR
        self.skip_compute_DGR_variability_profiling = A('skip_compute_DGR_variability_profiling')
        self.skip_primer_variability = A('skip_primer_variability')


        # performance
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        if self.num_threads:
            self.run.info('Threads Used', self.num_threads)
        self.run.info('BLASTn word size', self.word_size)
        self.run.info('Skip "N" characters', self.skip_Ns)
        self.run.info('Skip "-" characters', self.skip_dashes)
        self.run.info('Discovery mode', self.discovery_mode)
        self.run.info('Number of Mismatches', self.number_of_mismatches)
        self.run.info('Percentage of Mismatching Bases', self.percentage_mismatch)
        self.run.info('Minimum Mismatching Base Types in VR', self.min_mismatching_base_types_vr)
        self.run.info('Metagenomics Contigs Mode', self.metagenomics_contigs_mode)
        if self.metagenomics_contigs_mode:
            self.run.info('Collection(s) Provided', (self.collections_given))
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
        if not self.skip_recovering_genomic_context:
            self.run.info('Number of genes to consider in context', self.num_genes_to_consider_in_context)
        #computing variability profiling for every VR in every DGR by searching through raw reads?
        if not self.skip_compute_DGR_variability_profiling:
            self.run.info('Samples.txt', self.samples_txt)
            #self.run.info("R1/R2 for raw reads present?", "True" if self.raw_r1_r2_reads_are_present else "False")
            self.run.info('Variable Region Primer Base Length', self.variable_region_primer_base_length)



    def sanity_check(self):
        """Basic checks for a smooth operation"""

        try:
            output_dir = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=False or self.just_do_it)

            if output_dir:
                filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=False)
        except:
            raise ConfigError(f"Hold up your directory ({self.output_directory}) already exists. To avoid overwriting data we are stopping "
                            "you here. You have three options if you want to continue: rename your directory, delete your existing directory, "
                            "or rerun with the flag `--just-do-it` (which will overwrite your directory. You have been warned).")

        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                            "class, not multiple :/")

        if ((self.contigs_db_path and self.profile_db_path) or (self.contigs_db_path)) and not self.hmm:
            raise ConfigError("Just so you know in order to report a DGR you need to have Reverse Transcriptases "
                            "annotated or else you are not finding DGRs. If you would like to do this please run "
                            "your 'contigs.db' with `anvi-run-hmms -I Reverse_Transcriptase` (type: 6 clades of DGR "
                            "Retroelements from doi.org/10.1038/s41467-021-23402-7 including other known reverse "
                            "transcriptases). Then re-run `anvi-report-dgrs` with the `--hmm-usage` flag. This is "
                            "one option, but you can use other HMM profiles if you like.")

        if self.fasta_file_path:
            # check fasta input
            filesnpaths.is_file_fasta_formatted(self.fasta_file_path)

        if self.step < 0 or self.word_size < 0:
            raise ConfigError('The step value and/or word size value you are trying to input should be positive integer.')

        if self.variable_buffer_length < 0:
            raise ConfigError('The variable buffer length value you are trying to input should be positive integer.')

        if self.departure_from_reference_percentage < 0:
            raise ConfigError('The departure from reference percentage value you are trying to input should be a positive decimal number.')

        #if  self.metagenomics_contigs_mode:
            #if self.collections_provided not in self.collections_provided_in_profile:
                #raise ConfigError(f"You requested this collection was searched through: {self.collections_provided} in these collections {self.collections_provided_in_profile} "
                                        #f"that are in your {self.profile_db_path}. The collections you give 'anvi-report-dgrs' need to be in your "
                                        #"profile.db.")

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

        html_files_exist = any(file.endswith('.html') for file in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, file)))
        if html_files_exist:
            raise ConfigError("Files with .html suffix exist in the directory. Please delete before carrying on. (later this will delete them for you)")

        if not self.skip_compute_DGR_variability_profiling:
            self.samples_txt_dict = utils.get_samples_txt_file_as_dict(self.samples_txt)
            self.raw_r1_r2_reads_are_present = False
            self.raw_r1_r2_reads_are_present = all('r1' in v and 'r2' in v and isinstance(v['r1'], str) and isinstance(v['r2'], str) and v['r1'] and v['r2'] for v in self.samples_txt_dict.values())

            if not self.raw_r1_r2_reads_are_present and not self.skip_compute_DGR_variability_profiling:
                    raise ConfigError("You asked anvi'o to calculate DGR profiling variability across samples, but your samples-txt "
                                    "does not include raw R1/R2 reads :(")

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

            #open merged profile-db and get the variable nucleotide table as a dictionary then access the split names as a list to use in get_snvs
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            #Sort pandas data-frame of SNVs by contig name and then by position of SNV within contig
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
        This function is not currently in use but here for Merens' sake.
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
                        query_frame = str(hsp.find('Hsp_query-frame').text)
                        subject_frame = str(hsp.find('Hsp_hit-frame').text)

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
                    is_reverse_complement, TR_frame, VR_sequence, VR_frame, subject_genome_start_position, subject_genome_end_position, subject_contig, midline,
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
            midline : character with defined spacing
            percentage_of_mismatches : float
                Keys of a dictionary containing all of the BLASTn hits that are less than 100%

        Creates
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions and the corresponding info for those regions

        """

        DGR_key = f'DGR_{DGR_number:03d}'
        #print(f"adding new DGR, {DGR_key}. TR_is_query: {TR_is_query}")
        self.DGRs_found_dict[DGR_key] = {}
        # TR stuff
        self.DGRs_found_dict[DGR_key]['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['base'] = base

        # VR stuff
        self.DGRs_found_dict[DGR_key]['VRs'] = {'VR_001':{}}
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_sequence'] = VR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['midline'] = midline
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_reverse_complement'] = is_reverse_complement
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_frame'] = TR_frame
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['percentage_of_mismatches'] = percentage_of_mismatches
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_frame'] = VR_frame

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



    def update_existing_DGR(self, existing_DGR_key, TR_is_query, TR_frame, VR_sequence, VR_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement,
                            query_genome_start_position, query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position, subject_contig):
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
            midline : character with defined spacing
            percentage_of_mismatches : float
                Keys of a dictionary containing all of the BLASTn hits that are less than 100%

        Creates
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions and the corresponding info for those regions

        """
        # VR info
        num_VR = len(self.DGRs_found_dict[existing_DGR_key]['VRs']) + 1
        new_VR_key = f'VR_{num_VR:03d}'
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key] = {}
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_sequence'] = VR_sequence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['midline'] = midline
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['percentage_of_mismatches'] = percentage_of_mismatches
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_frame'] = VR_frame
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_reverse_complement'] = is_reverse_complement
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_frame'] = TR_frame

        # query-subject dependent stuff
        if TR_is_query:
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
        #NOTE: the arguments for add_new_DGR and update_existing_DGR that are based on the query and subject genome are positional, so they need
        # to be in the same place as those function deal with whether these are TR or VR. BUT the frame argument is not positional and you need
        # to change that based on whether it is in the query or subject. YES it is gross, but if it works, it works.

        """
        This function takes the none identical hits of the BLASTn and filters for template and variable regions.

        This works by filtering for sequences that have an over representation of one base that is mismatching and a certain number
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
            query_frame = int(hit_data['query_frame'])
            subject_frame = int(hit_data['subject_frame'])
            query_contig = hit_data['query_contig']
            subject_contig = hit_data['subject_contig']

            # get number of mismatches
            mismatch_length_bp = len(position)

            # if num of mismatches = 0, skip DGR search sanity check
            if mismatch_length_bp == 0:
                continue
            else:
                # Calculate the percentage identity of each alignment
                is_TR = False
                for letter, count in query_mismatch_counts.items():
                    percentage_of_mismatches = (count / mismatch_length_bp)
                    if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches):
                        is_TR = True
                        # if the letter is T, then we assume that it is an A base and we reverse complement EVERYTHING
                        if letter == 'T':
                            TR_sequence = str(query_sequence.reverse_complement())
                            VR_sequence = str(subject_sequence.reverse_complement())
                            midline = ''.join(reversed(original_midline))
                            base = 'A'
                            is_reverse_complement = True
                            query_frame = query_frame * -1
                            subject_frame = subject_frame * -1
                        else:
                            TR_sequence = str(query_sequence)
                            VR_sequence = str(subject_sequence)
                            midline = original_midline
                            base = letter
                            is_reverse_complement = False

                        #to test for VR diversity of base types in the protein sequence
                        for letter, count in subject_mismatch_counts.items():
                            non_zero_bases = sum(1 for count in subject_mismatch_counts.values() if count > 0)
                        if not non_zero_bases >= self.min_mismatching_base_types_vr:
                            continue

                        #need to check if the new TR you're looping through exists in the DGR_found_dict, see if position overlap
                        if not self.DGRs_found_dict:
                            # add first DGR
                            num_DGR += 1
                            self.add_new_DGR(num_DGR, True, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                        base, is_reverse_complement, query_frame, VR_sequence, subject_frame, subject_genome_start_position, subject_genome_end_position,
                                        subject_contig, midline, percentage_of_mismatches)
                        else:
                            was_added = False
                            for dgr in self.DGRs_found_dict:
                                if self.DGRs_found_dict[dgr]['TR_contig'] == query_contig and self.range_overlapping(query_genome_start_position,
                                                                                                                query_genome_end_position,
                                                                                                                self.DGRs_found_dict[dgr]['TR_start_position'],
                                                                                                                self.DGRs_found_dict[dgr]['TR_end_position']):
                                    was_added = True
                                    self.update_existing_DGR(dgr, True, query_frame, VR_sequence, subject_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement, query_genome_start_position,
                                                    query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                                                    subject_contig)
                                    break
                            if not was_added:
                                # add new TR and its first VR
                                num_DGR += 1
                                self.add_new_DGR(num_DGR, True, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                        base, is_reverse_complement, query_frame, VR_sequence, subject_frame, subject_genome_start_position, subject_genome_end_position,
                                        subject_contig, midline, percentage_of_mismatches)

                if not is_TR:
                    # Calculate the percentage identity of each alignment
                    for letter, count in subject_mismatch_counts.items():
                            percentage_of_mismatches = (count / mismatch_length_bp)
                            if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches):
                                is_TR = True
                                # if the letter is T, then we assume that it is an A base and we reverse complement EVERYTHING
                                if letter == 'T':
                                    TR_sequence = str(subject_sequence.reverse_complement())
                                    VR_sequence = str(query_sequence.reverse_complement())
                                    midline = ''.join(reversed(original_midline))
                                    base = 'A'
                                    is_reverse_complement = True
                                    query_frame = query_frame * -1
                                    subject_frame = subject_frame * -1
                                    #TODO: reverse the strand orientation of the TR and VR.
                                else:
                                    TR_sequence = str(subject_sequence)
                                    VR_sequence = str(query_sequence)
                                    midline = original_midline
                                    base = letter
                                    is_reverse_complement = False

                                #to test for VR diversity of base types in the protein sequence
                                for letter, count in query_mismatch_counts.items():
                                    non_zero_bases = sum(1 for count in query_mismatch_counts.values() if count > 0)
                                if not non_zero_bases >= self.min_mismatching_base_types_vr:
                                    continue
                                #need to check if the new TR you're looping through exists in the DGR_found_dict, see if position overlap
                                if not self.DGRs_found_dict:
                                    # add first DGR
                                    num_DGR += 1
                                    self.add_new_DGR(num_DGR, False, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                        base, is_reverse_complement, subject_frame, VR_sequence, query_frame, subject_genome_start_position, subject_genome_end_position,
                                        subject_contig, midline, percentage_of_mismatches)
                                else:
                                    was_added = False
                                    for dgr in self.DGRs_found_dict:
                                        if self.DGRs_found_dict[dgr]['TR_contig'] == subject_contig and self.range_overlapping(subject_genome_start_position,
                                                                                                                        subject_genome_end_position,
                                                                                                                        self.DGRs_found_dict[dgr]['TR_start_position'],
                                                                                                                        self.DGRs_found_dict[dgr]['TR_end_position']):
                                            was_added = True
                                            #TODO can rename consensus_TR
                                            self.update_existing_DGR(dgr, False, subject_frame, VR_sequence, query_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement, query_genome_start_position,
                                                    query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                                                    subject_contig)
                                            break
                                    if not was_added:
                                        # add new TR and its first VR
                                        num_DGR += 1
                                        self.add_new_DGR(num_DGR, False, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                        base, is_reverse_complement, subject_frame, VR_sequence, query_frame, subject_genome_start_position, subject_genome_end_position,
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
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Creates
        =======

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

                        # if there are function sources, let's recover them for our gene of interest
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

    # Function to add bin info if positions fall within bin ranges
    #TODO: have warning for DGR over two splits that are only in part of a collection. But ok if across 2 splits but both in same bin
    def add_bin_info(self, start_pos, end_pos, dgr_name, vr_name, bin_ranges_dict):
        tr_or_vr = 'vr' if vr_name else 'tr'
        bin_info = None
        for bin_name, ranges in bin_ranges_dict.items():
            for (bin_start, bin_end) in ranges:
                if start_pos >= bin_start and end_pos <= bin_end:
                    bin_info = bin_name
                    break
            if bin_info:
                break

        return bin_info


    def collections_mode(self):
            """
            This function is if the user only wants to search for DGRs that are in the same collection. This is known as the metagenomics mode.
            It filters through the dgrs found dictionary so that only the DGRs within the same collection are outputted.

            Parameters
            ==========
            DGRs_found_dict : dict
                A dictionary containing the template and variable regions

            Creates
            =======
            DGRs_found_dict : dict
                A dictionary containing the template and variable regions

            """

            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            contig_db = dbops.ContigsDatabase(self.contigs_db_path)

            #check that the collections provided are in the collections available in the profile.
            where_clause = f'''collection_name == "{self.collections_given}"'''
            self.split_collections_dict = profile_db.db.get_some_rows_from_table_as_dict(t.collections_splits_table_name, where_clause=where_clause,
                error_if_no_data=True)

            splits_list = []
            for key, value in self.split_collections_dict.items():
                splits_list.append(value['split'])

            #Construct the where_clause for the new SQL query
            splits_str = ', '.join(f'"{split}"' for split in splits_list)
            where_clause_splits = f'''split IN ({splits_str})'''

            # Retrieve rows from the splits basic info database using the new where_clause, to get the splits in the given collection
            self.splits_in_collections_dict = contig_db.db.get_some_rows_from_table_as_dict(t.splits_info_table_name, where_clause=where_clause_splits,
                error_if_no_data=True)

            # Add bin_name to the final dictionary
            for split, info in self.splits_in_collections_dict.items():
                for key, value in self.split_collections_dict.items():
                    if value['split'] == split:
                        info['bin_name'] = value['bin_name']
                        #need to also add in if start or end of split is inside the bin range
                        break

            from collections import defaultdict
            #create a dictionary of ranges of all the start and stop positions in a bin.
            bin_ranges_dict = defaultdict(list)

            for key, value in self.splits_in_collections_dict.items():
                bin_name = value['bin_name']
                start = value['start']
                end = value['end']
                bin_ranges_dict[bin_name].append((start, end))

            # Sort the ranges within each bin_name (if needed)
            for bin_name in bin_ranges_dict:
                bin_ranges_dict[bin_name].sort()

            for dgr_id, dgr_data in self.DGRs_found_dict.items():
                # Add bin info to TR
                tr_start, tr_end = dgr_data['TR_start_position'], dgr_data['TR_end_position']
                dgr_data['TR_bin'] = self.add_bin_info(tr_start, tr_end, dgr_id, None, bin_ranges_dict)

                for vr_id, vr_data in dgr_data['VRs'].items():
                    # Add bin info to VR
                    vr_start, vr_end = vr_data['VR_start_position'], vr_data['VR_end_position']
                    vr_data['VR_bin'] = self.add_bin_info(vr_start, vr_end, dgr_id, vr_id, bin_ranges_dict)

            #Initialize a new dictionary to store filtered DGRs
            self.dgrs_in_collections = {}

            # Iterate over DGRs_found_dict and filter TR and VR pairs in the same bin
            for dgr_name, dgr_data in self.DGRs_found_dict.items():
                tr_bin_info = dgr_data.get('TR_bin')
                print(f"Processing {dgr_name}: tr_bin_info = {tr_bin_info}")

                if tr_bin_info is None:
                    continue  # Skip if TR does not have bin info

                filtered_vrs = {}
                for vr_name, vr_data in dgr_data.get('VRs', {}).items():
                    vr_bin_info = vr_data.get('VR_bin')
                    print(f"Checking VR {vr_name} of {dgr_name}: vr_bin_info = {vr_bin_info}")
                    if vr_bin_info == tr_bin_info:
                        filtered_vrs[vr_name] = vr_data
                        print(f"Including VR {vr_name} in {dgr_name} because vr_bin {vr_bin_info} == tr_bin {tr_bin_info}")

                # Add to filtered_dgrs only if there are matching VRs
                if filtered_vrs:
                    self.dgrs_in_collections[dgr_name] = {**dgr_data, 'VRs': filtered_vrs}

                self.get_hmm_info(self.dgrs_in_collections)

                self.create_found_tr_vr_csv(self.dgrs_in_collections)

                self.recover_genomic_context_surrounding_dgrs(self.dgrs_in_collections)

                self.report_genomic_context_surrounding_dgrs(self.dgrs_in_collections)

                self.compute_dgr_variability_profiling(self.dgrs_in_collections)


            profile_db.disconnect()
            contig_db.disconnect()

            return

    def get_hmm_info(self, dgrs_dict):
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

        #look at general consensus TR in the level up so all the TRs have the same HMM if in the same DGR.
        for DGR_id, DGR_info in dgrs_dict.items():
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



    def create_found_tr_vr_csv(self, dgrs_dict):
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
        if not self.DGRs_found_dict and not self.dgrs_in_collections:
            raise ConfigError("No DGRS were found so no output file will be written :(")

        output_directory_path = self.output_directory

        if dgrs_dict == self.DGRs_found_dict:
            output_path_dgrs = os.path.join(output_directory_path, "DGRs_found.csv")
            headers = ["DGR", "VR", "VR_contig", "VR_strand_direction", "VR_sequence", "Midline", "VR_start_position", "VR_end_position","VR_bin", "Mismatch %",
                    "TR_contig", "TR_strand_direction", "TR_sequence", "Base", "Reverse Complement", "TR_start_position", "TR_end_position", "TR_bin", "HMM_source",
                    "distance_to_HMM", "HMM_gene_name", "HMM_direction", "HMM_start", "HMM_stop", "HMM_gene_callers_id"]
        elif dgrs_dict == self.dgrs_in_collections:
            # Create new directory for DGRs_found_in_collections
            self.collections_dir = os.path.join(output_directory_path, "DGRs_found_in_collections")
            if not os.path.exists(self.collections_dir):
                os.makedirs(self.collections_dir)
            output_path_dgrs = os.path.join(self.collections_dir, "DGRs_found_with_collections_mode.csv")
            headers = ["DGR", "VR", "VR_contig","VR_strand_direction", "VR_sequence", "Midline", "VR_start_position", "VR_end_position", "VR_bin",  "VR_frame","VR_reverse_comp_for_primer", "Mismatch %", "TR_contig", "TR_strand_direction", "TR_sequence", "Base", "Reverse Complement", "TR_start_position", "TR_end_position", "TR_bin", "TR_reverse_comp_for_primer", "HMM_source", "distance_to_HMM", "HMM_gene_name", "HMM_direction", "HMM_start", "HMM_stop", "HMM_gene_callers_id"]

        with open(output_path_dgrs, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)

            for dgr, tr in dgrs_dict.items():
                for vr, vr_data in tr['VRs'].items():
                    csv_row = [dgr, vr, vr_data['VR_contig'],vr_data['VR_frame'], vr_data['VR_sequence'], vr_data['midline'],
                            vr_data['VR_start_position'], vr_data['VR_end_position'], vr_data.get('VR_bin', 'N/A'), vr_data.get('VR_reverse_comp_for_primer', 'FALSE'), vr_data['percentage_of_mismatches'],
                            tr['TR_contig'],vr_data['TR_frame'], vr_data['TR_sequence'], tr['base'], vr_data['TR_reverse_complement'],
                            vr_data['TR_start_position'], vr_data['TR_end_position'], tr.get('TR_bin', 'N/A'), vr_data.get('TR_reverse_comp_for_primer', 'FALSE'), tr['HMM_source'], tr["distance_to_HMM"],
                            tr["HMM_gene_name"], tr["HMM_direction"], tr["HMM_start"], tr["HMM_stop"], tr["HMM_gene_callers_id"]]
                    csv_writer.writerow(csv_row)
            return



    def recover_genomic_context_surrounding_dgrs(self, dgrs_dict):
        """Learn about what surrounds the variable region sites of each found DGR"""

        # in which we will store the genomic context that surrounds dgrs for downstream fun
        self.genomic_context_surrounding_dgrs = {}
        self.genomic_context_surrounding_rts = {}

        # we know when we are not wanted
        if self.skip_recovering_genomic_context:
            self.run.info_single('Skipping genomic context recovery due to self.skip_recovering_genomic_context being True')
            return

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)

        # are there genes?
        if not contigs_db.meta['genes_are_called']:
            self.run.warning("There are no gene calls in your contigs database, therefore there is context to "
                            "learn about :/ Your reports will not include a file to study the genomic context "
                            "that surrounds the variable regions associated with the diversity-generating retroelements.")

            contigs_db.disconnect()
            return

        # are there functions?
        function_sources_found = contigs_db.meta['gene_function_sources'] or []
        if not len(function_sources_found):
            self.run.warning("There are no functions for genes in your contigs database :/ Your reports on the "
                            "genomic context that surrounds the variable regions associated with the diversity-generating retroelements will not have any functions "
                            "for genes. PITY.")
            return

        self.progress.new('Recovering genomic context surrounding the DGRs', progress_total_items=len(self.DGRs_found_dict))
        self.progress.update('...')

        # now we will go through each dgr to populate `self.genomic_context_surrounding_dgrs`
        # with gene calls and functions
        gene_calls_per_TR_contig = {}
        gene_calls_per_VR_contig = {}
        gene_calls_per_RT_contig = {}
        trs_with_no_gene_calls_around = set([])
        vrs_with_no_gene_calls_around = set([])
        rts_with_no_gene_calls_around = set([])


        for dgr_key, dgr_data in dgrs_dict.items():
            dgr_id = dgr_key
            self.progress.update(f"{dgr_id}", increment=True)

            TR_contig_name = dgr_data.get('TR_contig')
            TR_start = dgr_data.get('TR_start_position')
            TR_end = dgr_data.get('TR_end_position')

            RT_start = int(dgr_data.get('HMM_start'))
            RT_end = int(dgr_data.get('HMM_stop'))

            # Initialize TR context
            tr_context_genes = []

            # Initialize RT context
            rt_context_genes = []

            if TR_contig_name not in gene_calls_per_TR_contig:
                where_clause = f'''contig="{TR_contig_name}" and source="{self.gene_caller_to_consider_in_context}"'''
                gene_calls_per_TR_contig[TR_contig_name] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

            gene_calls_in_TR_contig = gene_calls_per_TR_contig[TR_contig_name]

            if not len(gene_calls_in_TR_contig):
                trs_with_no_gene_calls_around.add(dgr_id)
                print(f'No gene calls found around TR {dgr_id}')
                continue

            # Process TR gene calls
            min_distance_to_TR_start, min_distance_to_TR_end = float('inf'), float('inf')
            closest_gene_call_to_TR_start, closest_gene_call_to_TR_end = None, None
            for gene_callers_id, gene_call in gene_calls_in_TR_contig.items():
                if abs(gene_call['start'] - TR_start) < min_distance_to_TR_start:
                    closest_gene_call_to_TR_start = gene_callers_id
                    min_distance_to_TR_start = abs(gene_call['start'] - TR_start)

                if abs(gene_call['start'] - TR_end) < min_distance_to_TR_end:
                    closest_gene_call_to_TR_end = gene_callers_id
                    min_distance_to_TR_end = abs(gene_call['start'] - TR_end)

            TR_range = range(closest_gene_call_to_TR_start - self.num_genes_to_consider_in_context,
                            closest_gene_call_to_TR_end + self.num_genes_to_consider_in_context)
            tr_gene_caller_ids_of_interest = [c for c in TR_range if c in gene_calls_in_TR_contig]

            for gene_callers_id in tr_gene_caller_ids_of_interest:
                gene_call = gene_calls_in_TR_contig[gene_callers_id]
                gene_call['gene_callers_id'] = gene_callers_id

                if function_sources_found:
                    where_clause = '''gene_callers_id IN (%s)''' % (', '.join([f"{str(g)}" for g in tr_gene_caller_ids_of_interest]))
                    hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
                    gene_call['functions'] = [h for h in hits if h['gene_callers_id'] == gene_callers_id]

                dna_sequence = self.contig_sequences[TR_contig_name]['sequence'][gene_call['start']:gene_call['stop']]
                rev_compd = None
                if gene_call['direction'] == 'f':
                    gene_call['DNA_sequence'] = dna_sequence
                    rev_compd = False
                else:
                    gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)
                    rev_compd = True

                where_clause = f'''gene_callers_id == {gene_callers_id}'''
                aa_sequence = contigs_db.db.get_some_rows_from_table_as_dict(t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                gene_call['AA_sequence'] = aa_sequence[gene_callers_id]['sequence']

                gene_call['length'] = gene_call['stop'] - gene_call['start']

                header = '|'.join([f"contig:{gene_call['contig']}",
                                f"start:{gene_call['start']}",
                                f"stop:{gene_call['stop']}",
                                f"direction:{gene_call['direction']}",
                                f"rev_compd:{rev_compd}",
                                f"length:{gene_call['length']}"])
                gene_call['header'] = ' '.join([str(gene_callers_id), header])

                tr_context_genes.append(gene_call)

            self.genomic_context_surrounding_dgrs[dgr_id] = copy.deepcopy(tr_context_genes)


            if TR_contig_name not in gene_calls_per_RT_contig:
                where_clause = f'''contig="{TR_contig_name}" and source="{self.gene_caller_to_consider_in_context}"'''
                gene_calls_per_RT_contig[TR_contig_name] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

            gene_calls_in_RT_contig = gene_calls_per_RT_contig[TR_contig_name]

            if not len(gene_calls_in_RT_contig):
                rts_with_no_gene_calls_around.add(dgr_id)
                print(f'No gene calls found around RT {dgr_id}')
                continue

            # Process RT gene calls
            min_distance_to_RT_start, min_distance_to_RT_end = float('inf'), float('inf')
            closest_gene_call_to_RT_start, closest_gene_call_to_RT_end = None, None
            for gene_callers_id, gene_call in gene_calls_in_RT_contig.items():
                if abs(gene_call['start'] - RT_start) < min_distance_to_RT_start:
                    closest_gene_call_to_RT_start = gene_callers_id
                    min_distance_to_RT_start = abs(gene_call['start'] - RT_start)

                if abs(gene_call['start'] - RT_end) < min_distance_to_RT_end:
                    closest_gene_call_to_RT_end = gene_callers_id
                    min_distance_to_RT_end = abs(gene_call['start'] - RT_end)

            RT_range = range(closest_gene_call_to_RT_start - self.num_genes_to_consider_in_context,
                            closest_gene_call_to_RT_end + self.num_genes_to_consider_in_context)
            rt_gene_caller_ids_of_interest = [c for c in RT_range if c in gene_calls_in_RT_contig]

            for gene_callers_id in rt_gene_caller_ids_of_interest:
                gene_call = gene_calls_in_RT_contig[gene_callers_id]
                gene_call['gene_callers_id'] = gene_callers_id

                if function_sources_found:
                    where_clause = '''gene_callers_id IN (%s)''' % (', '.join([f"{str(g)}" for g in rt_gene_caller_ids_of_interest]))
                    hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
                    gene_call['functions'] = [h for h in hits if h['gene_callers_id'] == gene_callers_id]

                dna_sequence = self.contig_sequences[TR_contig_name]['sequence'][gene_call['start']:gene_call['stop']]
                rev_compd = None
                if gene_call['direction'] == 'f':
                    gene_call['DNA_sequence'] = dna_sequence
                    rev_compd = False
                else:
                    gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)
                    rev_compd = True

                where_clause = f'''gene_callers_id == {gene_callers_id}'''
                aa_sequence = contigs_db.db.get_some_rows_from_table_as_dict(t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                gene_call['AA_sequence'] = aa_sequence[gene_callers_id]['sequence']

                gene_call['length'] = gene_call['stop'] - gene_call['start']

                header = '|'.join([f"contig:{gene_call['contig']}",
                                f"start:{gene_call['start']}",
                                f"stop:{gene_call['stop']}",
                                f"direction:{gene_call['direction']}",
                                f"rev_compd:{rev_compd}",
                                f"length:{gene_call['length']}"])
                gene_call['header'] = ' '.join([str(gene_callers_id), header])

                rt_context_genes.append(gene_call)

            self.genomic_context_surrounding_rts[dgr_id] = copy.deepcopy(rt_context_genes)


            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key

                self.progress.update(f"{dgr_id} , {vr_id}", increment=True)

                VR_contig = vr_data.get('VR_contig')
                VR_start = vr_data.get('VR_start_position')
                VR_end = vr_data.get('VR_end_position')

                # Initialize VR context
                vr_context_genes = []

                if VR_contig not in gene_calls_per_VR_contig:
                    where_clause = f'''contig="{VR_contig}" and source="{self.gene_caller_to_consider_in_context}"'''
                    gene_calls_per_VR_contig[VR_contig] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

                gene_calls_in_VR_contig = gene_calls_per_VR_contig[VR_contig]

                if not len(gene_calls_in_VR_contig):
                    vrs_with_no_gene_calls_around.add(vr_id)
                    print(f'No gene calls found around DGR {dgr_id} VR {vr_id}')
                    continue

                min_distance_to_VR_start, min_distance_to_VR_end = float('inf'), float('inf')
                closest_gene_call_to_VR_start, closest_gene_call_to_VR_end = None, None
                for gene_callers_id, gene_call in gene_calls_in_VR_contig.items():
                    if abs(gene_call['start'] - VR_start) < min_distance_to_VR_start:
                        closest_gene_call_to_VR_start = gene_callers_id
                        min_distance_to_VR_start = abs(gene_call['start'] - VR_start)

                    if abs(gene_call['start'] - VR_end) < min_distance_to_VR_end:
                        closest_gene_call_to_VR_end = gene_callers_id
                        min_distance_to_VR_end = abs(gene_call['start'] - VR_end)

                VR_range = range(closest_gene_call_to_VR_start - self.num_genes_to_consider_in_context,
                                closest_gene_call_to_VR_end + self.num_genes_to_consider_in_context)
                vr_gene_caller_ids_of_interest = [c for c in VR_range if c in gene_calls_in_VR_contig]

                for gene_callers_id in vr_gene_caller_ids_of_interest:
                    gene_call = gene_calls_in_VR_contig[gene_callers_id]
                    gene_call['gene_callers_id'] = gene_callers_id

                    if function_sources_found:
                        where_clause = '''gene_callers_id IN (%s)''' % (', '.join([f"{str(g)}" for g in vr_gene_caller_ids_of_interest]))
                        hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
                        gene_call['functions'] = [h for h in hits if h['gene_callers_id'] == gene_callers_id]

                    dna_sequence = self.contig_sequences[VR_contig]['sequence'][gene_call['start']:gene_call['stop']]
                    rev_compd = None
                    if gene_call['direction'] == 'f':
                        gene_call['DNA_sequence'] = dna_sequence
                        rev_compd = False
                    else:
                        gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)
                        rev_compd = True

                    where_clause = f'''gene_callers_id == {gene_callers_id}'''
                    aa_sequence = contigs_db.db.get_some_rows_from_table_as_dict(t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                    gene_call['AA_sequence'] = aa_sequence[gene_callers_id]['sequence']

                    gene_call['length'] = gene_call['stop'] - gene_call['start']

                    header = '|'.join([f"contig:{gene_call['contig']}",
                                    f"start:{gene_call['start']}",
                                    f"stop:{gene_call['stop']}",
                                    f"direction:{gene_call['direction']}",
                                    f"rev_compd:{rev_compd}",
                                    f"length:{gene_call['length']}"])
                    gene_call['header'] = ' '.join([str(gene_callers_id), header])

                    vr_context_genes.append(gene_call)

                self.genomic_context_surrounding_dgrs[vr_id] = copy.deepcopy(vr_context_genes)

        contigs_db.disconnect()
        self.progress.end()
        print('Completed recovering genomic context surrounding the DGRs')

        self.run.info(f"[Genomic Context] Searched for {PL('DGR', len(dgrs_dict))}",
                    f"Recovered for {PL('TR', len(self.genomic_context_surrounding_dgrs[dgr_id]))}",
                    f"And recovered for {PL('VR', len(self.genomic_context_surrounding_dgrs[vr_id]))}",
                    f"And recovered for {PL('RT',len(self.genomic_context_surrounding_rts[dgr_id]))}",
                    nl_before=1,
                    lc="yellow")

        if len(trs_with_no_gene_calls_around):
            print('No gene calls around the following TRs:', trs_with_no_gene_calls_around, "Here is the list in case you would like to track them down: "
                            f"{', '.join(trs_with_no_gene_calls_around)}.")
        if len(vrs_with_no_gene_calls_around):
            print('No gene calls around the following VRs:', vrs_with_no_gene_calls_around, "Here is the list in case you would like to track them down: "f"{', '.join(vrs_with_no_gene_calls_around)}.")

        if len(rts_with_no_gene_calls_around):
            print('No gene calls around the following RTs:', rts_with_no_gene_calls_around, "Here is the list in case you would like to track them down: "f"{', '.join(rts_with_no_gene_calls_around)}.")

        if not len(self.genomic_context_surrounding_dgrs):
            self.run.warning(f"Even though the tool went through all {PL('DGR', len(dgrs_dict))} "
                            f"it was unable to recover any genomic context for any of them. So your final reports will "
                            f"not include any insights into the surrounding genomic context of all your DGRs ( a little bit sad "
                            f"but otherwise you will be fine).")



    def report_genomic_context_surrounding_dgrs(self, dgrs_dict):
        """
        Reports two long-format output files for genes and functions around inversion
        STOLEN (modified) FROM INVERSIONS CODE (line 1925)
        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======


        """

        if self.skip_recovering_genomic_context:
            self.run.info_single('Skipping reporting genomic context due to self.skip_recovering_genomic_context being True')
            return

        if not len(self.genomic_context_surrounding_dgrs):
            print("No genomic context data available to report")
            return

        # we are in business
        genes_output_headers = ["gene_callers_id", "start", "stop", "direction", "partial", "call_type", "source", "version", "contig"]
        functions_output_headers = ["gene_callers_id", "source", 'accession', 'function']

        # Process each DGR and its VRs
        for dgr_key, dgr_data in dgrs_dict.items():
            # Assuming dgr_key itself is the dgr_id or a dictionary containing it
            dgr_id = dgr_key  # If dgr_key is the dgr_id itself

            # Create output directory for DGR
            if dgrs_dict == self.DGRs_found_dict:
                dgr_directory = os.path.join(self.output_directory, "PER_DGR", dgr_id)
            elif dgrs_dict == self.dgrs_in_collections:
                dgr_directory = os.path.join(self.collections_dir, "PER_DGR", dgr_id)

            filesnpaths.gen_output_directory(dgr_directory, delete_if_exists=False)

            # TR output paths
            tr_genes_output_path = os.path.join(dgr_directory, 'TR_SURROUNDING-GENES.txt')
            tr_functions_output_path = os.path.join(dgr_directory, 'TR_SURROUNDING-FUNCTIONS.txt')

            with open(tr_genes_output_path, 'w') as tr_genes_output, open(tr_functions_output_path, 'w') as tr_functions_output:
                tr_genes_output.write("dgr_id\tentry_type\t%s\n" % '\t'.join(genes_output_headers))
                tr_functions_output.write("dgr_id\t%s\n" % '\t'.join(functions_output_headers))

                # Create fake gene call entries for TR and VRs:
                d = dict([(h, '') for h in genes_output_headers])

                # Fill in non-empty data for the TR in the DGR and insert it:
                d['contig'] = dgr_data.get('TR_contig')  # Use dgr_data to access TR info
                d['start'] = dgr_data.get('TR_start_position')
                d['stop'] = dgr_data.get('TR_end_position')
                tr_genes_output.write(f"{dgr_id}_TR\tTEMPLATE_REGION\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                #Check if there are surrounding genes for the TR and write them
                if dgr_id in self.genomic_context_surrounding_dgrs:
                    for gene_call in self.genomic_context_surrounding_dgrs[dgr_id]:
                        tr_genes_output.write(f"{dgr_id}_TR\tGENE\t%s\n" % '\t'.join([f"{gene_call[h]}" for h in genes_output_headers]))

                        if 'functions' in gene_call:
                            for hit in gene_call['functions']:
                                tr_functions_output.write(f"{dgr_id}_TR\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                        else:
                            tr_functions_output.write(f"{dgr_id}_TR\t{gene_call['gene_callers_id']}\t\t\t\n")

                self.run.info(f"Reporting file on gene context for {dgr_id} TR", tr_genes_output_path)
                self.run.info(f"Reporting file on functional context for {dgr_id} TR", tr_functions_output_path, nl_after=1)

            # RT output paths
            rt_genes_output_path = os.path.join(dgr_directory, 'RT_SURROUNDING-GENES.txt')
            rt_functions_output_path = os.path.join(dgr_directory, 'RT_SURROUNDING-FUNCTIONS.txt')

            with open(rt_genes_output_path, 'w') as rt_genes_output, open(rt_functions_output_path, 'w') as rt_functions_output:
                rt_genes_output.write("dgr_id\tentry_type\t%s\n" % '\t'.join(genes_output_headers))
                rt_functions_output.write("dgr_id\t%s\n" % '\t'.join(functions_output_headers))

                # Create fake gene call entries for RTs:
                d = dict([(h, '') for h in genes_output_headers])

                # Fill in non-empty data for the TR in the DGR and insert it:
                d['contig'] = dgr_data.get('TR_contig')  # Use dgr_data to access TR info
                d['start'] = dgr_data.get('HMM_start')
                d['stop'] = dgr_data.get('HMM_stop')
                rt_genes_output.write(f"{dgr_id}_RT\tREVERSE_TRANSCRIPTASE\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                #Check if there are surrounding genes for the TR and write them
                if dgr_id in self.genomic_context_surrounding_rts:
                    for gene_call in self.genomic_context_surrounding_rts[dgr_id]:
                        rt_genes_output.write(f"{dgr_id}_RT\tGENE\t%s\n" % '\t'.join([f"{gene_call[h]}" for h in genes_output_headers]))

                        if 'functions' in gene_call:
                            for hit in gene_call['functions']:
                                rt_functions_output.write(f"{dgr_id}_RT\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                        else:
                            rt_functions_output.write(f"{dgr_id}_RT\t{gene_call['gene_callers_id']}\t\t\t\n")

                self.run.info(f"Reporting file on gene context for {dgr_id} RT", rt_genes_output_path)
                self.run.info(f"Reporting file on functional context for {dgr_id} RT", rt_functions_output_path, nl_after=1)


            # Fill in non-empty data for each VR in the DGR and insert it:
            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key
                vr_directory = os.path.join(dgr_directory, f"VR_{vr_id}")
                filesnpaths.gen_output_directory(vr_directory, delete_if_exists=False)

                vr_genes_output_path = os.path.join(vr_directory, 'VR_SURROUNDING-GENES.txt')
                vr_functions_output_path = os.path.join(vr_directory, 'VR_SURROUNDING-FUNCTIONS.txt')

                with open(vr_genes_output_path, 'w') as vr_genes_output, open(vr_functions_output_path, 'w') as vr_functions_output:
                    vr_genes_output.write("dgr_id\tentry_type\t%s\n" % '\t'.join(genes_output_headers))
                    vr_functions_output.write("dgr_id\t%s\n" % '\t'.join(functions_output_headers))

                    # Create fake gene call entry for VR:
                    d['contig'] = vr_data.get('VR_contig')
                    d['start'] = vr_data.get('VR_start_position')
                    d['stop'] = vr_data.get('VR_end_position')
                    vr_genes_output.write(f"{dgr_id} {vr_id}\tVARIABLE_REGION\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                    # Check if there are surrounding genes for the VR and write them
                    if vr_id in self.genomic_context_surrounding_dgrs:
                        for gene_call in self.genomic_context_surrounding_dgrs[vr_id]:
                            vr_genes_output.write(f"{dgr_id} VR_{vr_id}\tGENE\t%s\n" % '\t'.join([f"{gene_call[h]}" for h in genes_output_headers]))

                            if 'functions' in gene_call:
                                for hit in gene_call['functions']:
                                    vr_functions_output.write(f"{dgr_id} VR_{vr_id}\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                            else:
                                vr_functions_output.write(f"{dgr_id} VR_{vr_id}\t{gene_call['gene_callers_id']}\t\t\t\n")

                    self.run.info(f'Reporting file on gene context for {dgr_id} {vr_id}', vr_genes_output_path)
                    self.run.info(f'Reporting file on functional context for {dgr_id} {vr_id}', vr_functions_output_path, nl_after=1)



    # Function to get the consensus base
    def get_consensus_base(row):
        for nucleotide in nucleotides:
            if row[nucleotide] > 0.5:  # Assuming a threshold for consensus, adjust as necessary
                return nucleotide
        return None



    @staticmethod
    def compute_dgr_variability_profiling_per_vr(input_queue, output_queue, samples_dict, primers_dict, output_directory_path, run=run_quiet, progress=progress_quiet):
        """
        Go back to the raw metagenomic reads to compute the variability profiles of the variable regions for each single Sample.

        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions
        Reads_1 : fasta file
            A fasta file containing the forward reads used to create the merged PROFILE.db
        Reads_2 : fasta file
            A fasta file containing the reverse reads used to create the merged PROFILE.db
        Samples.txt : txt file
            Contains sample information and the path to both of the read files

        Returns
        =======


        """
        while True:
            sample_name = input_queue.get(True)
            print(f"(I am in compute per sample func) Processing sample: {sample_name}")
            if sample_name is None:
                print('Sample is none, loop will be broken.')
                break

            samples_dict_for_sample = {sample_name: samples_dict[sample_name]}

            # setup the args object
            args = argparse.Namespace(samples_dict= samples_dict_for_sample,
                                    primers_dict= primers_dict,
                                    #min_remainder_length= self.primer_remainder_lengths,
                                    #min_frequency=min_frequency,
                                    output_dir=output_directory_path,
                                    only_report_primer_matches = True
                                    )

            s = PrimerSearch(args, run=run, progress=progress)
            sample_dict, primer_hits = s.process(return_dicts = True)

            output_queue.put(sample_name)


    def generate_primers_for_vrs(self, dgrs_dict):
        """
        A function to generate primers for each and every VR. These are composed of not only an initial primer sequence before the VR
        but also of anchor points in the VR. These anchor points are the bases in the TR that are not A bases and only at the places the TR and VR match.

                    TR:GCTAACTGACATAATT
        Anchor_primer :GCT..C.G.C.T..TT
                    VR:GCTCACGGACTTCATT

        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======
        primers_dict : dict
            A dictionary containing the various primers for each sample including the compositional sections of each primer (i.e. anchor primer and the initial primer)
        """

        #create primers dictionary
        primers_dict = {}

        #get contigs.db and contig sequences
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)

        # need to get the length of each consensus dgr's tr to have a set length for each VR profile in every sample so that they are the same
        #self.primer_remainder_lengths = {}

        for dgr_id, dgr_data in dgrs_dict.items():
            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key

                #save original frames
                VR_FRAME = vr_data['VR_frame']
                TR_FRAME = vr_data['TR_frame']
                vr_data['original_VR_frame'] = VR_FRAME
                vr_data['original_TR_frame'] = TR_FRAME
                print(f"dgrs dit: {dgrs_dict}")

                VR_frame = vr_data['VR_frame']
                TR_frame = vr_data['TR_frame']

                #CHECK if VR sequence is not at the start of a contig so you can get the initial primer sequence
                if vr_data['VR_start_position'] >= self.variable_region_primer_base_length:
                    vr_primer_region_start = vr_data['VR_start_position'] - self.variable_region_primer_base_length
                    vr_primer_region_end = (vr_data['VR_start_position'] -1)
                else:
                    #this will take the start of the contig to create a shorter initial primer length
                    vr_primer_region_start = vr_data['VR_start_position']
                    vr_primer_region_end = (vr_data['VR_start_position'] -1)
                    self.run.warning(f"The primer sequence for this VR {dgr_id}_{vr_id}. Is going to fall off the start of the Contig it is one: {vr_data['VR_contig']}.")

                contig_sequence = self.contig_sequences[vr_data['VR_contig']]['sequence']
                vr_primer_region = contig_sequence[vr_primer_region_start:vr_primer_region_end]
                #add every primer sequence to the dgrs_dict
                vr_data['vr_primer_region'] = vr_primer_region

                print(f"original dgr dict: {dgrs_dict}")
                print('\n')
                VR_sequence = vr_data['VR_sequence']
                TR_sequence = vr_data['TR_sequence']
                vr_start = vr_data.get('VR_start_position')
                vr_end = vr_data.get('VR_end_position')
                vr_contig = vr_data.get('VR_contig')

                #FIRST CHECK THAT YOU HAVEN'T FOUND THE REVERSE COMPLEMENT OF THE STRAND ANYWAYS
                #flip them so that they are the original way that blast found them!
                #this will be gross but you need to keep the original keys so therefore create objects to flip
                if vr_data['TR_reverse_complement'] == True and vr_data['VR_frame'] == -1 and vr_data['TR_frame'] == -1:
                    VR_frame = vr_data['VR_frame'] * -1
                    VR_sequence = utils.rev_comp(vr_data['VR_sequence'])
                    TR_frame = vr_data['TR_frame'] * -1
                    TR_sequence = utils.rev_comp(vr_data['TR_sequence'])


                print(f"AFTER TR REV COMP CHECK:\n{vr_id} VR sequence: {VR_sequence}\n VR frame {VR_frame}")

                print(f"AFTER TR REV COMP CHECK:\n{dgr_id} TR sequence: {TR_sequence}\n  TR frame {TR_frame}")

                #make every frame positive so that the initial primer is always on the left and then both the vr and tr are being compared on the same strand.
                #Always make the VR strand +1 so that you can compare the primer to the fasta file by definition
                if vr_data['VR_frame'] == -1:
                    print(f"I am reverse complementing {dgr_id} {vr_id} AND THE TR because both are -1")
                    vr_data['rev_comp_VR_seq'] = utils.rev_comp(vr_data['VR_sequence'])
                    VR_frame = VR_frame * -1
                    vr_data['VR_reverse_comp_for_primer'] = True
                    VR_sequence = str(vr_data['rev_comp_VR_seq'])
                    vr_data['rev_comp_TR_seq'] = utils.rev_comp(vr_data['TR_sequence'])
                    TR_frame = TR_frame * -1
                    vr_data['TR_reverse_comp_for_primer'] = True
                    TR_sequence = str(vr_data['rev_comp_TR_seq'])
                    print(f"AFTER VR + TR = -1:\n{vr_id} VR sequence: {VR_sequence}\n VR frame {VR_frame}")
                    print(f"AFTER VR + TR = -1:\n{dgr_id} TR sequence: {TR_sequence}\n  TR frame {TR_frame}")


                print(f"FINAL:\n{vr_id} VR sequence: {VR_sequence}\n VR frame {VR_frame}")
                print(f"FINAL:\n{dgr_id} TR sequence: {TR_sequence}\n  TR frame {TR_frame}")

                #print(f"Updated dgr dict after frame handling: {dgrs_dict}\n")

                #check the TR and VR sequence are the same length
                if len(TR_sequence) == len(VR_sequence):
                    profile_db = dbops.ProfileDatabase(self.profile_db_path)
                    snvs_table = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
                    profile_db.disconnect()

                    # Extract the contig name from split_name
                    snvs_table['contig_name'] = snvs_table['split_name'].str.split('_split_', expand=True)[0]


                    print("snvs table")
                    print(snvs_table)

                    vr_positions = set(range(vr_start, vr_end))

                    # Filter snvs_table to include only rows with the matching VR_contig
                    filtered_snvs_table = snvs_table[(snvs_table['contig_name'] == vr_contig) & (snvs_table['pos_in_contig'].isin(vr_positions))]

                    # Convert the positions to a set for faster lookups
                    primer_snv_positions = set(filtered_snvs_table['pos_in_contig'])


                    # Optional: Print to verify
                    print(f"VR_contig: {vr_contig}")
                    print(f"Filtered positions for VR_contig: {primer_snv_positions}")
                    print(f"vr positions lists: {vr_positions}")


                    vr_primer = []
                    # Create the vr_primer sequence
                    # First loop: Create the initial vr_primer sequence
                    for tr_base, vr_base in zip(TR_sequence, VR_sequence):
                        if tr_base == 'A':
                            vr_primer.append('.')
                        elif tr_base != vr_base:
                            vr_primer.append('.')
                        else:
                            vr_primer.append(tr_base)

                    #Second loop: Apply SNV-based modifications
                    for i in range(len(vr_primer)):
                        # Calculate the current position in the VR sequence
                        current_position = vr_start + i

                        # Check if the position should be replaced with '.'
                        if current_position in primer_snv_positions:
                            vr_primer[i] = '.'

                    # Convert list back to string
                    vr_anchor_primer = ''.join(vr_primer)

                    # Add the vr_anchor_primer sequence to the dgrs_dict
                    vr_data['vr_anchor_primer'] = vr_anchor_primer

                    print(dgrs_dict)

                elif len(TR_sequence) != len(VR_sequence):
                    print(f"Thats weird! The {vr_id} does not have the same length as the {dgr_id}'s TR :( so you can't create an anchor primer sequence")

        ###########################################
        # UPDATED DGRs dict with Primer Sequences #
        ###########################################

        for dgr_id, dgr_data in dgrs_dict.items():
            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key
                primers_dict[dgr_id + '_' + vr_id + '_Primer'] = {'initial_primer_sequence': vr_data['vr_primer_region'],
                                                                'vr_anchor_primer': vr_data['vr_anchor_primer']}
        print('\n')
        print(f'original primers dictionary:' , primers_dict)

        return primers_dict



    def print_primers_dict_to_csv(self, primers_dict):
        output_directory_path = self.output_directory
        output_path_for_genes_found = os.path.join(output_directory_path, "DGR_Primers_used_for_VR_diversity.csv")

        # Define the header for the CSV file
        csv_header = ['Primer_ID', 'Initial Primer', 'Anchor Primer', 'Whole Primer']

        # Open the CSV file in write mode
        with open(output_path_for_genes_found, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(csv_header)  # Write the header row
            # Iterate through the dictionary and write each gene's information to the CSV file
            for primer_name, primer_info in primers_dict.items():
                    initial_primer = primer_info['initial_primer_sequence']
                    anchor_primer = primer_info['vr_anchor_primer']
                    whole_primer = primer_info['primer_sequence']


                    writer.writerow([
                    primer_name,
                    initial_primer,
                    anchor_primer,
                    whole_primer
                    ])
        return

    def compute_dgr_variability_profiling(self, dgrs_dict):
        """
        Go back to the raw metagenomic reads to compute the variability profiles of the variable regions
        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions
        Reads_1 : fasta file
            A fasta file containing the forward reads used to create the merged PROFILE.db
        Reads_2 : fasta file
            A fasta file containing the reverse reads used to create the merged PROFILE.db
        Samples.txt : txt file
            Contains sample information and the path to both of the read files

        Returns
        =======

        """

        if self.skip_compute_DGR_variability_profiling or not self.raw_r1_r2_reads_are_present:
            return
        print(f'self.samples_txt_dict:', self.samples_txt_dict)

        if not len(self.DGRs_found_dict):
            self.run.info_single("Compute DGR variability profile function speaking: There are no DGRs to "
                                "compute in-sample variability :/", mc="red")

        sample_names = list(self.samples_txt_dict.keys())
        num_samples = len(sample_names)

        print('sample names below')
        print(sample_names)

        # let the user know what is going on
        msg = (f"Now anvi'o will compute in-sample activity of {PL('DGR VR', len(self.DGRs_found_dict))} "
            f"across {PL('sample', num_samples)}. Brace yourself and please note that this can "
            f"take a very long time since for each sample, anvi'o will go through each short read to search for ever variable region "
            f"sequence/s per DGR. You can always skip this step and search for individual primers "
            f"listed in the output file using the program `anvi-search-primers` with the parameter "
            f"`--min-remainder-length` set to the length of the consensus template region of the DGRs' VR you are interested in "
            "and the flag `--only-report-remainders` to explore the variable region activity of that one variable region "
            f"manually. Maybe this makes no sense? See the documentation for `anvi-report-dgrs` (and hope for "
            f"the best)")
        self.run.warning(None, header="PERFORMANCE NOTE", lc="yellow")

        if num_samples > self.num_threads:
            self.run.info_single(f"You have {PL('sample', num_samples)} but {PL('thread', self.num_threads)}. Therefore, not all samples will be processed "
                                f"in parallel. Just an FYI. {msg}.", level=0, nl_after=1)
        elif self.num_threads > num_samples:
            self.run.info_single(f"You have {PL('sample', num_samples)} but {PL('thread', self.num_threads)}. Since only samples are run in "
                                f"parallel, the additional {PL('thread', self.num_threads - num_samples)} you have there are not really "
                                f"useful for anything. Just an FYI. {msg}.", level=0, nl_after=1)
            self.num_threads = num_samples
        else:
            self.run.info_single(f"{msg}.", level=0, nl_after=1)

        # here we will need to reconstruct a samples_dict and primers_dict to pass to the class
        # `PrimerSearch`. for this we first need to generate a list of primers. for each
        # DGRs VR, which at this point are described in a dict like this,
        # {'DGR_001': {'TR_sequence': 'GCGGCTCCTGGAACAACTATCCTAGGAGGTGTCGCTCTGCGAACCGCAACAACTATAACTCGGACGAGGCGGACAACAACAATATTGGTTTTCGTCTTGTGAGT',
        #                'base': 'A',
        #               'TR_reverse_complement': True,
        #               'VRs': {'VR_003':
        #                          {'VR_sequence': 'GCGGCTCCTGGTACGACTTTCCTTGGTGGTGTCGCTCTGCGTTCCGCGGCTACTATTTCTCGGTCGAGGCGGTCAACGACTTTGTTGGTTTTCGTCTTGTGAGT',
        #                           'TR_sequence': 'GCGGCTCCTGGAACAACTATCCTAGGAGGTGTCGCTCTGCGAACCGCAACAACTATAACTCGGACGAGGCGGACAACAACAATATTGGTTTTCGTCTTGTGAGT',
        #                           'midline': '||||||||||| || ||| |||| || ||||||||||||||  ||||  | |||||  ||||| |||||||| |||| ||  | ||||||||||||||||||||',
        #                           'percentage_of_mismatches': 1.0,
        #                           'VR_start_position': 2631629,
        #                           'VR_end_position': 2631732,
        #                           'VR_contig': 'T_erythraeum_IMS101_000000000001',
        #                           'VR_sequence_found': 'query',
        #                           'TR_start_position': 2632660,
        #                           'TR_end_position': 2632763,
        #                           'VR_bin':'DGR2_meren'},
        #                       'VR_004':
        #                           {'VR_sequence': 'GCGGCTCCTGGCTCAACTATCCTTGGTGGTGTCGCTCTGCGTACCGCTACGACTTTAGCTCGGACGGGGCGGTCATCATCAATTTTGGTTTTCGTCTTGTGAGT',
        #                            'TR_sequence': 'GCGGCTCCTGGAACAACTATCCTAGGAGGTGTCGCTCTGCGAACCGCAACAACTATAACTCGGACGAGGCGGACAACAACAATATTGGTTTTCGTCTTGTGAGT',
        #                            'midline': '|||||||||||  |||||||||| || |||||||||||||| ||||| || ||| || |||||||| ||||| || || |||| ||||||||||||||||||||',
        #                             'percentage_of_mismatches': 1.0,
        #                             'VR_start_position': 2634431,
        #                             'VR_end_position': 2634534,
        #                             'VR_contig': 'T_erythraeum_IMS101_000000000001',
        #                             'VR_sequence_found': 'query',
        #                             'TR_start_position': 2632660,
        #                             'TR_end_position': 2632763,
        #                             'VR_bin': 'DGR2_meren'}},
        #               'TR_start_position': 2632660,
        #               'TR_end_position': 2632763,
        #               'TR_contig': 'T_erythraeum_IMS101_000000000001',
        #               'TR_sequence_found': 'subject',
        #               'TR_bin': 'DGR2_meren',
        #               'HMM_gene_callers_id': '1688',
        #               'distance_to_HMM': 434.5,
        #               'HMM_start': 2632879,
        #               'HMM_stop': 2633413,
        #               'HMM_direction': 'r',
        #               'HMM_source': 'Reverse_Transcriptase',
        #               'HMM_gene_name': 'Clean_DGR_clade_3',
        #               'HMM_gene_source': 'prodigal'}}
        #
        # need primers for each VR
        # either we have the same number of primers as VRs
        # OR we have one primer per VR. I think this might be more suitable for ease of code and what we want the output to look like

        # FIRST we need to get the primers sequences!
        primers_dict = self.generate_primers_for_vrs(dgrs_dict)

        #function to find variability in primers and form consensus primer for each primer in each sample!
        if not self.skip_primer_variability:
            sample_primers_dict = {}
            sample_names = list(self.samples_txt_dict.keys())
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            snvs_table = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
            profile_db.disconnect()

            for sample_name in sample_names:
                # Filter SNVs for the specific sample
                sample_snvs = self.snv_panda[self.snv_panda['sample_id'] == sample_name]

                # Iterate through dgrs_dict to get the VR_start and VR_contig
                for dgr_id, dgr_data in dgrs_dict.items():
                    for vr_key, vr_data in dgr_data['VRs'].items():
                        vr_start = vr_data.get('VR_start_position')
                        vr_end = vr_data.get('VR_end_position')
                        vr_contig = vr_data.get('VR_contig')

                        #sanity check for mismatch between samples given and samples in snv table
                        sample_names_given = set(sample_names)
                        print(f"sample names given from samples.txt: {sample_names_given}")
                        sample_names_in_snv_table = set(snvs_table['sample_id'])
                        print(f"sample names in snv table {sample_names_in_snv_table}")
                        samples_missing_in_snv_table = sample_names_given.difference(sample_names_in_snv_table)
                        print(f"samples_missing_in_snv_table {samples_missing_in_snv_table}")
                        if anvio.DEBUG:
                            self.run.info("Samples given", ", ".join(list(sample_names_given)))
                            self.run.info("Samples in profile.db's nucleotide variability table", ", ".join(list(sample_names_in_snv_table)))
                            self.run.info("Missing samples from profile.db's nucleotide variability table", ", ".join(list(samples_missing_in_snv_table)))
                        if sample_names_given == samples_missing_in_snv_table:
                            raise ConfigError(f"Anvi'o is not angry just disappointed :/ You gave 'anvi-report-dgrs' these samples ({list(sample_names_given)}), yet you have none of these samples in your profile.db. "
                                                "This is fatal, anvio will now quit. Either you recreate your profile.db with the samples you would like to search for the DGR VRs variability, "
                                                "or you give 'anvi-report-dgrs' samples that were used to create your profile.db.")

                        for sample_name in sample_names:
                            if sample_name not in snvs_table['sample_id']:
                                self.run.warning(f"This is unacceptable, not all samples are in the SNV table, missing {sample_name}. BUT we will keep going anyways, just for you.")
                                continue
                            # Filter SNVs for the specific sample
                            sample_snvs = snvs_table[snvs_table['sample_id'] == sample_name]

                            # Further filter SNVs within the primer region
                            primer_snvs = sample_snvs[
                                (sample_snvs['contig_name'] == vr_contig) &
                                (sample_snvs['pos_in_contig'] >= vr_start) &
                                (sample_snvs['pos_in_contig'] < vr_end)
                            ]

                        # Get the original primer sequence
                        original_primer_key = f'{dgr_id}_{vr_key}_Primer'
                        original_initial_primer_sequence = primers_dict[original_primer_key]['initial_primer_sequence']
                        new_primer_sequence = list(original_initial_primer_sequence)

                        # Vectorized operation to find consensus SNVs and update the primer sequence
                        consensus_snvs = primer_snvs[primer_snvs['departure_from_reference'] > 0.5].apply(self.get_consensus_base, axis=1)
                        positions_in_primer = vr_start - primer_snvs[primer_snvs['departure_from_reference'] > 0.5]['pos_in_contig'] - 1

                        for position, consensus_base in zip(positions_in_primer, consensus_snvs):
                            if consensus_base and 0 <= position < len(new_primer_sequence):
                                new_primer_sequence[position] = consensus_base


                        ## Add the new primer sequence to the sample-specific primers dictionary
                        dgr_vr_key = f'{dgr_id}_{vr_key}'
                        if dgr_vr_key not in sample_primers_dict:
                            sample_primers_dict[dgr_vr_key] = {}
                        sample_primers_dict[dgr_vr_key][sample_name] = {
                            'initial_primer_sequence': ''.join(new_primer_sequence),
                            #'primer_remainder_TR_length': primers_dict[original_primer_key]['primer_remainder_TR_length']
                        }

                print(f'Sample {sample_name} updated primers dictionary:', sample_primers_dict)

        #update the primers dictionary with the total primer sequence including any edits the initial primer accounting for SNVs
        for dgr_id, dgr_data in dgrs_dict.items():
            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key
                primer_key = f'{dgr_id}_{vr_id}_Primer'
                # Set the final primer sequence in primers_dict based on initial sequence and variability analysis
                primers_dict[primer_key]['primer_sequence'] = primers_dict[primer_key]['initial_primer_sequence'] + primers_dict[primer_key]['vr_anchor_primer']

                #if the primer sequence is too long than cut it to an appropriate length for the length of short reads
                if len(primers_dict[primer_key]['primer_sequence']) > 71:
                    self.chosen_primer_length = 65
                    print(f"The primer for {dgr_id} {vr_id} is above the desired length for getting enough diversity. The primer is being cut to {self.chosen_primer_length}. Nothing is wrong anvi'o just wants to make you aware.")
                    primers_dict[primer_key]['primer_sequence'] = primers_dict[primer_key]['primer_sequence'][:self.chosen_primer_length]

        print('\n')
        print(f'Updated primers dictionary:', primers_dict)

        self.run.info_single("Computing the Variable Regions Primers and creating a 'DGR_Primers_used_for_VR_diversity.csv' file.")
        self.print_primers_dict_to_csv(primers_dict)

        ##################
        # MULTITHREADING #
        ##################

        # setup the input/output queues
        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        self.dgr_activity = []
        #create directory for Primer matches
        primer_folder= os.path.join(self.output_directory, "PRIMER_MATCHES")


        # put all the sample names in our input queue
        for sample_name in sample_names:
            input_queue.put(sample_name)

        # engage the proletariat, our hard-working wage-earner class
        workers = []
        for i in range(self.num_threads):
            print(f"starting worker {i}")
            worker = multiprocessing.Process(target=DGR_Finder.compute_dgr_variability_profiling_per_vr,
                                            args=(input_queue,
                                                output_queue,
                                                self.samples_txt_dict,
                                                primers_dict,
                                                primer_folder),

                                            kwargs=({'progress': self.progress if self.num_threads == 1 else progress_quiet}))
            workers.append(worker)
            worker.start()

        self.progress.new('DGR variability profile', progress_total_items=num_samples)
        if self.num_threads > 1:
            self.progress.update(f"Processing {PL('sample', num_samples)} and {PL('primer', len(primers_dict))} in {PL('thread', self.num_threads)}.")

        num_samples_processed = 0
        while num_samples_processed < num_samples:
            try:
                sample_finished_processing = output_queue.get()
                if anvio.DEBUG:
                    self.progress.reset()
                    self.run.info_single(f"Sample {sample_finished_processing} has finished processing.")
                num_samples_processed += 1
                self.progress.increment(increment_to=num_samples_processed)
                if self.num_threads > 1:
                    if num_samples_processed < num_samples:
                        self.progress.update(f"Samples processed: {num_samples_processed} of {num_samples}. Still working ...")
                    else:
                        self.progress.update("All done!")
            except KeyboardInterrupt:
                self.run.info_single("Received kill signal, terminating all processes... Don't believe anything you see "
                                    "below this and destroy all the output files with fire.", nl_before=1, nl_after=1)
                break

        if self.num_threads > 1:
            self.progress.end()

        # always double-tap?
        for worker in workers:
            worker.terminate()

        print(f'FINAL PRIMERS DICT AFTER PROCESSING', primers_dict)

        ######################
        # END MULTITHREADING #
        ######################



    def process_dgr_data_for_HTML_summary(self):
        """Take everything that is known, turn them into data that can be used from Django templates.

        A lot of ugly/frightening things happening here to prepare coordinates for SVG objects to be displayed
        or store boolean variables to use the Django template engine effectively. IF YOU DON'T LIKE
        IT DON'T LOOK AT IT. IT MIGHT MAKE YOU CRY
        """

        # in which we will store all the static HTML output related stuff
        self.summary = {}

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.summary['meta'] = {'summary_type': 'dgrs',
                                #'num_dgrs': len(self.dgrs_in_collections) if self.metagenomics_contigs_mode else len(self.DGRs_found_dict),
                                'num_dgrs': len(self.DGRs_found_dict),
                                #'num_samples': len(self.profile_db_paths) if self.metagenomics_contigs_mode else len(self.collections_given),
                                'output_directory': self.output_directory,
                                'genomic_context_recovered': not self.skip_recovering_genomic_context,
                                #'inversion_activity_computed': not self.skip_compute_inversion_activity,
                                # if no function source, it says 'the contigs.db' because it fits with the message
                                # displayed in the final index.html. See the inversion template, line 215
                                # if it works, it works
                                'gene_function_sources': contigs_db.meta['gene_function_sources'] or ['the contigs.db']}
        contigs_db.disconnect()

        print("Configured Summary Type:", self.summary['meta']['summary_type'])


        self.summary['files'] = {'Putative_DGRs': 'Putative-DGRs.txt'}
        self.summary['dgrs'] = {}

        for dgr_key, dgr_data in self.DGRs_found_dict.items():
            # Assuming dgr_key itself is the dgr_id or a dictionary containing it
            dgr_id = dgr_key

            self.summary['dgrs'][dgr_id] = {'dgr_data': copy.deepcopy(dgr_data)}

            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key

                # Ensure the dictionary structure for the VR exists
                if 'VRs' not in self.summary['dgrs'][dgr_id]['dgr_data']:
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'] = {}

                self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id] = {}

                if self.skip_recovering_genomic_context:
                    pass
                else:
                    # we will get a deepcopy of the gene context associated with the inversion
                    tr_genes = copy.deepcopy(self.genomic_context_surrounding_dgrs[dgr_id])
                    vr_genes = copy.deepcopy(self.genomic_context_surrounding_dgrs[dgr_id])

                    # then we will learn about these so we can transform the coordinates of anything we wish
                    # to display in the output
                    genomic_context_start_tr = tr_genes[0]['start'] - 100
                    genomic_context_end_tr = tr_genes[-1]['stop'] + 100
                    genomic_context_start_vr = vr_genes[0]['start'] - 100
                    genomic_context_end_vr = vr_genes[-1]['stop'] + 100

                    # this is our magic number, which is matching to the actual width of the genomic context
                    # display in the static HTML output. we will have to transform start-stop coordinates
                    # of each gene to this value.
                    new_context_length = 1000

                    # how big the gene arrows should be (in an ideal world -- see below the real world, Neo)
                    default_gene_arrow_width = 20

                    # before we start working on the genes, we will figure out the location of the inverted site
                    # in the genomic context. here we quickly identify the transformed start and the end position
                    # and store it in the dgr data dict
                    tr_start = (int(dgr_data['TR_start_position']) - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length
                    tr_end  = (int(dgr_data['TR_end_position']) - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length
                    self.summary['dgrs'][dgr_id]['dgr_data']['TX'] = tr_start
                    self.summary['dgrs'][dgr_id]['dgr_data']['TW'] = tr_end - tr_start
                    self.summary['dgrs'][dgr_id]['dgr_data']['TT'] = tr_start + (tr_end - tr_start) / 2

                    vr_start = (int(vr_data['VR_start_position']) - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length
                    vr_end  = (int(vr_data['VR_end_position']) - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['VX'] = vr_start
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['VW'] = vr_end - vr_start
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['VT'] = vr_start + (vr_end - vr_start) / 2

                    # here we will add transformed gene coordinates to the genes dict
                    for gene in tr_genes:
                        gene['start_tr_g'] = (gene['start'] - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length
                        gene['stop_tr_g'] = (gene['stop'] - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length

                        if (gene['stop_tr_g'] - gene['start_tr_g']) < default_gene_arrow_width:
                            # if we are here, it means the transformed length of the gene is already
                            # shorter than the space we assign for the arrow to display gene calls.
                            # this means we will only will be able to show an arrow, but even in that
                            # case the `gene_arrow_width` may be too long to display (i.e., if the
                            # transformed gene length is 10 and arrow is 15, we are already showing
                            # too much). The solution is to make the gene nothing more but the arrow
                            # but make the arrow width equal to the gene width
                            gene_arrow_width = gene['stop_tr_g'] - gene['start_tr_g']
                            gene['stop_tr_g'] = gene['start_tr_g']
                            gene['TRW'] = 0
                        else:
                            gene_arrow_width = default_gene_arrow_width
                            gene['TRW'] = (gene['stop_tr_g'] - gene['start_tr_g']) - gene_arrow_width

                        if 'functions' in gene.keys():
                            gene['has_functions'] = True
                            gene['COLOR'] = '#008000'
                        else:
                            gene['has_functions'] = False
                            gene['COLOR'] = '#c3c3c3'

                        gene['TRX'] = gene['start_tr_g']
                        gene['TCX'] = (gene['start_tr_g'] + (gene['stop_tr_g'] - gene['start_tr_g']) / 2)
                        gene['TGY'] = gene['TRX'] + gene['TRW'] + gene_arrow_width
                        gene['TGTRANS'] = gene['TRX'] + gene['TRX'] + gene['TRW'] + gene_arrow_width
                        gene['TRX_TRW'] = gene['TRX'] + gene['TRW'] - 0.5 # <-- minus 0.5 makes the arrow nicely cover the rest of the gene

                        for gene in vr_genes:
                            gene['start_vr_g'] = (gene['start'] - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length
                            gene['stop_vr_g'] = (gene['stop'] - genomic_context_start_vr) / (genomic_context_end_tr - genomic_context_start_vr) * new_context_length

                            if (gene['stop_vr_g'] - gene['start_vr_g']) < default_gene_arrow_width:
                                # if we are here, it means the transformed length of the gene is already
                                # shorter than the space we assign for the arrow to display gene calls.
                                # this means we will only will be able to show an arrow, but even in that
                                # case the `gene_arrow_width` may be too long to display (i.e., if the
                                # transformed gene length is 10 and arrow is 15, we are already showing
                                # too much). The solution is to make the gene nothing more but the arrow
                                # but make the arrow width equal to the gene width
                                gene_arrow_width = gene['stop_vr_g'] - gene['start_vr_g']
                                gene['stop_vr_g'] = gene['start_vr_g']
                                gene['VRW'] = 0
                            else:
                                gene_arrow_width = default_gene_arrow_width
                                gene['VRW'] = (gene['stop_vr_g'] - gene['start_vr_g']) - gene_arrow_width

                            if 'functions' in gene.keys():
                                gene['has_functions'] = True
                                gene['COLOR'] = '#008000'
                            else:
                                gene['has_functions'] = False
                                gene['COLOR'] = '#c3c3c3'

                            gene['VRX'] = gene['start_vr_g']
                            gene['VCX'] = (gene['start_vr_g'] + (gene['stop_vr_g'] - gene['start_vr_g']) / 2)
                            gene['VGY'] = gene['VRX'] + gene['VRW'] + gene_arrow_width
                            gene['VGTRANS'] = gene['VRX'] + gene['VRX'] + gene['VRW'] + gene_arrow_width
                            gene['VRX_VRW'] = gene['VRX'] + gene['VRW'] - 0.5 # <-- minus 0.5 makes the arrow nicely cover the rest of the gene

                            # and finally we will store this hot mess in our dictionary
                            self.summary['dgrs'][dgr_id]['tr_genes'] = tr_genes
                            self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['vr_genes'] = vr_genes

                            # also we need the path to the output files
                            self.summary['files'][dgr_id] = {'tr_genes': os.path.join('PER_DGR', dgr_id, 'SURROUNDING-GENES.txt'),
                                                                'functions': os.path.join('PER_DGR', dgr_id, 'SURROUNDING-FUNCTIONS.txt')}
                            self.summary['files'][dgr_id] = {'vr_genes': os.path.join('PER_DGR', vr_id, 'SURROUNDING-GENES.txt'),
                                                                        'functions': os.path.join('PER_DGR', vr_id, 'SURROUNDING-FUNCTIONS.txt')}

        # Ensure the destination directory does not exist before generating the summary HTML
        destination_dir = 'summary_html_output'
        if os.path.exists(destination_dir):
            shutil.rmtree(destination_dir)

        summary_html_output = SummaryHTMLOutput(self.summary)
        summary_html_output.generate('summary_html_output')

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
                ("Skip '-'", self.skip_dashes if self.skip_dashes else "FALSE"),
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
                ("Metagenomics Contigs Mode", self.metagenomics_contigs_mode if self.metagenomics_contigs_mode else "FALSE"),
                ("Output Directory", self.output_directory if self.output_directory else "default")
            ]

            csv_writer.writerows(parameters)
        return


    def process(self,args):
        """Here we process all of the functions in our class and call upon different functions depending on the inputs used"""
        self.sanity_check()
        self.get_blast_results()
        self.filter_blastn_for_none_identical()
        self.filter_for_TR_VR()
        if self.metagenomics_contigs_mode:
            self.run.info_single("Running metagenomics mode")
            self.collections_mode()
        if args.parameter_output:
            self.run.info_single("Writing to Parameters used file.")
            self.run.info_single('\n')
            self.parameter_output_sheet()

        if self.fasta_file_path:
            return

        else:
            self.run.info_single("Computing the Genes the Variable Regions occur in and creating a 'DGR_genes_found.csv'.")
            print('\n')
            self.get_gene_info()
            self.run.info_single("Computing the closest HMMs to the Template Regions and printing them in your output csv.")
            self.get_hmm_info(self.DGRs_found_dict)
            self.create_found_tr_vr_csv(self.DGRs_found_dict) #TODO: check if this works if the dictionary is empty or do you need a different method?
            self.recover_genomic_context_surrounding_dgrs(self.DGRs_found_dict)
            self.report_genomic_context_surrounding_dgrs(self.DGRs_found_dict)
            self.compute_dgr_variability_profiling(self.DGRs_found_dict) # add if statement to this for the metagenomics mode DGRs
            #self.process_dgr_data_for_HTML_summary()
        return