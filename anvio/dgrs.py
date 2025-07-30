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
import pytrf # type: ignore

from Bio import SeqIO # type: ignore
from Bio.SeqRecord import SeqRecord # type: ignore
from Bio.Seq import Seq # type: ignore

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths
import anvio.tables as t
import multiprocess as multiprocessing # type: ignore


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
        self.word_size = A('word_size')
        self.skip_Ns = A('skip_Ns')
        self.skip_dashes = A('skip_dashes')
        self.num_threads = A('num-threads')
        self.number_of_mismatches = A('number_of_mismatches')
        self.percentage_mismatch = A('percentage_mismatch')
        self.min_mismatching_base_types_vr = A('min_mismatching_base_types_vr') or 2
        self.min_base_types_vr = A('min_base_types_vr') or 2
        self.min_base_types_tr = A('min_base_types_tr') or 2
        self.only_a_bases =A('only_a_bases')
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
        self.collections_mode = A('collections_mode')
        self.collections_given = A('collection_name')
        self.skip_recovering_genomic_context = A('skip_recovering_genomic_context')
        self.num_genes_to_consider_in_context = A('num_genes_to_consider_in_context') or 3
        self.samples_txt = A('samples_txt')
        self.whole_primer_length = A('whole_primer_length') or 65
        self.initial_primer_length = A('initial_variable_primer_length') or 12 #TODO test different values for this. If Illumina reads are 250 bases then depends on length of VR
        self.skip_compute_DGR_variability_profiling = A('skip_compute_DGR_variability_profiling')
        self.skip_primer_variability = A('skip_primer_variability')
        self.numb_imperfect_tandem_repeats = A('numb_imperfect_tandem_repeats') or 10
        self.repeat_motif_coverage = A('repeat_motif_coverage') or 0.8
        self.snv_matching_proportion = A('snv_matching_proportion') or None
        self.snv_codon_position = A('snv_codon_position') or 0.33 #default is 33% of SNVs in the third codon position

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
        if self.only_a_bases:
            self,run.info('only A base mismatches', self.only_a_bases)
        self.run.info('Minimum Mismatching Base Types in VR', self.min_mismatching_base_types_vr)
        self.run.info('Minimum Base Types in VR', self.min_base_types_vr)
        self.run.info('Minimum Base Types in VR', self.min_base_types_tr)
        self.run.info('Number of imperfect tandem repeats', self.numb_imperfect_tandem_repeats)
        self.run.info('Collections Mode', self.collections_mode)
        if self.collections_mode:
            self.run.info('Collection(s) Provided', (self.collections_given))
        self.run.info('Output Directory', self.output_directory)
        self.run.info('Gene Caller Provided', self.gene_caller_to_consider_in_context)
        self.run.info('Contigs.db', self.contigs_db_path)
        self.run.info('Profile.db', self.profile_db_path)
        self.run.info('Minimum distance between SNVs', self.min_dist_bw_snvs)
        self.run.info('Minimum length of SNV window', self.min_range_size)
        self.run.info('Variable buffer length', self.variable_buffer_length)
        self.run.info('Departure from reference percentage', self.departure_from_reference_percentage)
        if self.snv_matching_proportion:
            self.run.info('SNV matching proportion', self.snv_matching_proportion)
        self.run.info('HMM(s) Provided', ", ".join(self.hmm))
        if not self.skip_recovering_genomic_context:
            self.run.info('Number of genes to consider in context', self.num_genes_to_consider_in_context)
        #computing variability profiling for every VR in every DGR by searching through raw reads?
        if not self.skip_compute_DGR_variability_profiling:
            self.run.info('Samples.txt', self.samples_txt)
            self.run.info('Initial Primer Length', self.initial_primer_length)
            self.run.info('Variable Region Primer Length', self.whole_primer_length)
            #self.run.info("R1/R2 for raw reads present?", "True" if self.raw_r1_r2_reads_are_present else "False")



    def sanity_check(self):
        """Basic checks for a smooth operation"""

        #First check the contigs.db and profile.db exists
        utils.is_contigs_db(self.contigs_db_path)
        utils.is_profile_db(self.profile_db_path)

        try:
            output_dir = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=False or self.just_do_it)

            if output_dir:
                filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=False)
        except:
            raise ConfigError(f"Hold up your directory ({self.output_directory}) already exists. To avoid overwriting data we are stopping "
                            "you here. You have three options if you want to continue: rename your directory, delete your existing directory, "
                            "or rerun with the flag `--just-do-it` (which will overwrite your directory. You have been warned).")

        if (self.contigs_db_path and self.profile_db_path) and not self.hmm:
            raise ConfigError("Just so you know in order to report a DGR you need to have Reverse Transcriptases "
                            "annotated or else you are not finding DGRs. If you would like to do this please run "
                            "your 'contigs.db' with `anvi-run-hmms -I Reverse_Transcriptase` (type: 6 clades of DGR "
                            "Retroelements from doi.org/10.1038/s41467-021-23402-7 including other known reverse "
                            "transcriptases). Then re-run `anvi-report-dgrs` with the `--hmm-usage` flag. This is "
                            "one option, but you can use other HMM profiles if you like.")

        if int(self.word_size) < 0:
            raise ConfigError('The word size value you are trying to input should be positive integer.')

        if self.variable_buffer_length < 0:
            raise ConfigError('The variable buffer length value you are trying to input should be positive integer.')

        if self.departure_from_reference_percentage < 0:
            raise ConfigError('The departure from reference percentage value you are trying to input should be a positive decimal number.')

        if self.min_mismatching_base_types_vr >= 5:
            raise ConfigError('The number of mismatching base types of the VR sequence cannot exceed 4 this is because there are only 4 bases in our DNA alphabet')

        if self.min_base_types_tr >= 5:
            raise ConfigError('The number of base types of the TR sequence cannot exceed 4 this is because there are only 4 bases in our DNA alphabet')

        if self.min_base_types_vr >= 5:
            raise ConfigError('The number of base types of the VR sequence cannot exceed 4 this is because there are only 4 bases in our DNA alphabet')

        if self.collections_mode:
            if not self.collections_given:
                raise ConfigError("You must provide a collection name for collections mode to work. If you want to know about "
                                "all the collections in your profile database you can use the program "
                                "`anvi-show-collections-and-bins` or run the same command with the flag "
                                "`--list-collections`.")
            # Ensure collections_given is a single string (collection name)
            if not isinstance(self.collections_given, str):
                raise ValueError("'collection-name' must be a single collection name")

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
            raise ConfigError("Files with .html suffix exist in the directory. Please delete them before rerunning and we will keep calm and carry on. (later this will delete them for you)")

        if not self.skip_compute_DGR_variability_profiling and not self.samples_txt:
            raise ConfigError("No samples.txt declared and no skip-compute-DGR-variability-profiling flag used. Either use the skip flag or "
                            "instruct anvi'o where to find the samples.txt that you want to profile are.")

        if not self.skip_compute_DGR_variability_profiling and self.samples_txt:
            self.samples_txt_dict = utils.get_samples_txt_file_as_dict(self.samples_txt)
            self.raw_r1_r2_reads_are_present = False
            self.raw_r1_r2_reads_are_present = all('r1' in v and 'r2' in v and isinstance(v['r1'], str) and isinstance(v['r2'], str) and v['r1'] and v['r2'] for v in self.samples_txt_dict.values())

            if not self.raw_r1_r2_reads_are_present and not self.skip_compute_DGR_variability_profiling:
                    raise ConfigError("You asked anvi'o to calculate DGR profiling variability across samples, but your samples-txt "
                                    "does not include raw R1/R2 reads :(")



    def get_blast_results(self):
        """
        This function runs the BLASTn search, which uses sequences with high SNVs as the query and the contigs.db sequences as the target sequences.

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

        if self.collections_mode:
            self.run.info_single("Collections mode activated. Get ready to see as many BLASTn as bins in your collection. Big things be happenin'.")
            #I know we don't need contigs sequences here as everything is done for splits, but I think I need this here fore the rest of the code to function as normal
            #TODO: Need to see if this works or if I need to change the rest of the code to work with splits
            #TODO: add warning if a contig is over multiple splits. Because you blast the splits gaianst contigs.
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
            contigs_db.disconnect()

            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            #check that the collections provided are in the collections available in the profile.
            where_clause = f'''collection_name == "{self.collections_given}"'''
            self.split_collections_dict = profile_db.db.get_some_rows_from_table_as_dict(t.collections_splits_table_name, where_clause=where_clause,
                error_if_no_data=True)
            profile_db.disconnect()

            # Dictionary to store splits by bin_name
            bin_splits_dict = {}

            # Iterate through the entries in self.split_collections_dict
            for key, value in self.split_collections_dict.items():
                bin_name = value['bin_name']
                split = value['split']

                # Initialize the list for the bin if it doesn't exist yet
                if bin_name not in bin_splits_dict:
                    bin_splits_dict[bin_name] = []

                # Append the split to the corresponding bin's list
                bin_splits_dict[bin_name].append(split)
            self.bin_names_list = list(bin_splits_dict.keys())

            # Now process each bin and create subsets based on the splits list
            for bin_name, bin_splits_list in bin_splits_dict.items():

                # Extract the contig names from split names
                bin_contigs = [split.split('_split_')[0] for split in bin_splits_list]

                # Subset self.contig_sequences using the extracted contig names
                bin_contig_sequences = {contig: self.contig_sequences[contig] for contig in bin_contigs if contig in self.contig_sequences}

                profile_db = dbops.ProfileDatabase(self.profile_db_path)
                #Sort pandas data-frame of SNVs by contig name and then by position of SNV within contig
                self.snv_panda_bin = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
                self.snv_panda_bin['contig_name'] = self.snv_panda_bin.split_name.str.split('_split_').str[0]

                self.split_names_unique = bin_splits_list

                profile_db.disconnect()

                # Filter snv_panda by the list of splits associated with this bin
                self.snv_panda = self.snv_panda_bin[self.snv_panda_bin['split_name'].isin(bin_splits_list)]

                sample_id_list = list(set(self.snv_panda.sample_id.unique()))

                self.all_possible_windows = {} # we will keep this as a dictionary that matches contig name to list of window tuples within that contig
                # structure of self.all_possible_windows: {'contig_0' : [(start0, stop0), (start1, stop1), (start2, stop2), ....... ],
                #                                           'contig_1': [(start0, stop0), (start1, stop1), (start2, stop2), ....... ],
                #                                           ....}
                # the windows will not necessarily be sorted within each inner list (yet) because we add windows from one sample at a time

                # filter snv_panda for min departure from ref and codon position
                if self.discovery_mode:
                    self.run.info("Running discovery mode. Search for SNVs in all possible locations. You go Dora the explorer!")
                    self.snv_panda = self.snv_panda.query("departure_from_reference >= @self.departure_from_reference_percentage")
                else:
                    self.snv_panda = self.snv_panda.query(
                        "departure_from_reference >= @self.departure_from_reference_percentage and base_pos_in_codon in (1, 2)"
                    )

                # First, group the DataFrame by 'split_name' and 'sample_id' upfront
                grouped = self.snv_panda.groupby(['split_name', 'sample_id'])

                # Now iterate over the grouped data
                for (split, sample), group in grouped:
                    # Only process groups that match the desired split and sample
                    if split in self.split_names_unique and sample in sample_id_list:
                        if group.shape[0] == 0:
                            continue

                        # Extract the contig name and positions for the group
                        contig_name = group.contig_name.unique()[0]
                        pos_list = group.pos_in_contig.to_list()

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

                            contig_len = len(bin_contig_sequences[contig_name]['sequence'])

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
                utils.export_sequences_from_contigs_db(self.contigs_db_path, self.target_file_path, seq_names_to_export=sorted(list(bin_contig_sequences.keys())))

                #get short sequences from all_merged_snv_window and create new fasta from self.target_file_path
                contig_records = []
                for record in SeqIO.parse(self.target_file_path, "fasta"):
                    contig_name = record.id
                    if contig_name in all_merged_snv_windows:
                        self.positions= all_merged_snv_windows[contig_name]
                        for i, (start, end) in enumerate(self.positions):
                            section_sequence = record.seq[start:end]
                            section_id = f"{contig_name}_section_{i}_start_bp{start}_end_bp{end}"
                            contig_records.append(SeqRecord(section_sequence, id=section_id, description=""))

                if len(contig_records) == 0:
                    self.run.warning(f"No sequences with SNVs were found with the parameters minimum distance between SNVs:{self.min_dist_bw_snvs} "
                                    f"and range size of SNVs:{self.min_range_size}, this means there are no variable region candidates for a blast search"
                                    , header="NO SEQUENCES WITH SUBSTANTIAL SNVS FOUND")
                    raise ConfigError(f"Therefore, we will exit here because anvi'o has nothing to search for DGRs in, "
                                    "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                                    "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")


                # Save the subset sequences to a temporary FASTA file
                output_fasta_path = os.path.join(self.temp_dir, f"bin_{bin_name}_subsequences.fasta")
                with open(output_fasta_path, "w") as output_handle:
                    SeqIO.write(contig_records, output_handle, "fasta")

                self.run.info('Temporary (SNV window) query input for blast', output_fasta_path)

                self.blast_output = os.path.join(tmp_directory_path,f"blast_output_for_bin_{bin_name}_step_{self.step}_wordsize_{self.word_size}.xml")
                blast = BLAST(output_fasta_path, target_fasta = self.target_file_path, search_program = 'blastn',
                            output_file=self.blast_output, additional_params = '-dust no', num_threads=self.num_threads)
                blast.evalue = 10 #set Evalue to be same as blastn default
                blast.makedb(dbtype = 'nucl')
                blast.blast(outputfmt = '5', word_size = self.word_size)

            profile_db.disconnect()
        #else not in collections mode
        else:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
            contigs_db.disconnect()

            #open merged profile-db and get the variable nucleotide table as a dictionary then access the split names as a list to use in get_snvs
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            #Sort pandas data-frame of SNVs by contig name and then by position of SNV within contig
            self.snv_panda = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
            self.snv_panda['contig_name'] = self.snv_panda['split_name'].apply(lambda x: x.split('_split_')[0])
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
                self.run.info("Running discovery mode. Search for SNVs in all possible locations. You go Dora the explorer!")
                self.snv_panda = self.snv_panda.query("departure_from_reference >= @self.departure_from_reference_percentage")
            else:
                self.snv_panda = self.snv_panda.query(
                    "departure_from_reference >= @self.departure_from_reference_percentage and base_pos_in_codon in (1, 2)"
                )

                # First, group the DataFrame by 'split_name' and 'sample_id' upfront
                grouped = self.snv_panda.groupby(['split_name', 'sample_id'])

                # Now iterate over the grouped data
                for (split, sample), group in grouped:
                    # Only process groups that match the desired split and sample
                    if split in self.split_names_unique and sample in sample_id_list:
                        if group.shape[0] == 0:
                            continue

                        # Extract the contig name and positions for the group
                        contig_name = group.contig_name.unique()[0]
                        pos_list = group.pos_in_contig.to_list()

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
                    for i, (start, end) in enumerate(self.positions):
                        section_sequence = record.seq[start:end]
                        section_id = f"{contig_name}_section_{i}_start_bp{start}_end_bp{end}"
                        contig_records.append(SeqRecord(section_sequence, id=section_id, description=""))

            if len(contig_records) == 0:
                self.run.warning(f"No sequences with SNVs were found with the parameters minimum distance between SNVs:{self.min_dist_bw_snvs} "
                                f"and range size of SNVs:{self.min_range_size}, this means there are no variable region candidates for a blast search"
                                , header="NO SEQUENCES WITH SUBSTANTIAL SNVS FOUND")
                raise ConfigError(f"Therefore, we will exit here because anvi'o has nothing to search for DGRs in, "
                                    "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                                    "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")

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
        #NOTE: An old function that is here just in case she is ever needed again.
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



    def process_blast_results(self, max_percent_identity=100):
        """
        Process BLAST output depending on whether collections_mode is enabled.
        If collections_mode is True, process multiple bins separately and merge results.
        If collections_mode is False, process a single BLAST output normally.
        This function takes the BLASTn xml output and refines the results to those with less than 100% identity.

        Takes the xml file and filters for hits with less than 100% identity, then gives every hit a name
        with its original position in the sequence, counts the bases that are mismatching and on which strand they occur.
        Finally initializes all of these within a dictionary.

        Parameters
        ==========
        blast_output : xml file
            BLASTn results
        max_percent_identity: number
            Maximum percentage identity threshold (default 100%).

        Returns
        mismatch_hits : dict
            A dictionary of all of the BLASTn hits that are less than 100%
        """

        # If collections_mode is enabled, process multiple bins separately
        if self.collections_mode:
            self.merged_mismatch_hits = {}  # Store combined results

            tmp_directory_path = self.temp_dir

            for bin_name in self.bin_names_list:

                # Reset mismatch hits for each bin
                self.mismatch_hits = {}

                blast_file = os.path.join(
                    tmp_directory_path,
                    f"blast_output_for_bin_{bin_name}_step_{self.step}_wordsize_{self.word_size}.xml"
                )

                if not os.path.exists(blast_file):
                    self.run.warning(f"Warning: BLAST output file for {bin_name} not found. Skipping...")
                    continue

                if os.stat(f"{self.blast_output}").st_size == 0:
                    self.run.warning(f"No DGR like sequences are being found via BLAST.", header="NO DGRS FOUND")
                    raise ConfigError(f"Therefore, we will exit here because anvi'o has found no DGRs in your data, "
                                        "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                                        "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")

                # Parse XML file for the current bin
                tree = ET.parse(blast_file)
                root = tree.getroot()

                # Process mismatches from this bin
                self.process_single(root, bin_name, max_percent_identity)

                # Merge results into `merged_mismatch_hits`
                for hit_identity, hit_data in self.mismatch_hits.items():
                    self.merged_mismatch_hits.setdefault(hit_identity, []).append(hit_data)

            self.run.info_single(f"Total unique mismatches: {len(self.merged_mismatch_hits)}")
            return self.merged_mismatch_hits

        #run in normal none collections mode
        else:
            # Single BLAST output processing
            self.mismatch_hits = {}

            if os.stat(f"{self.blast_output}").st_size == 0:
                self.run.warning(f"No DGR like sequences are being found via BLAST.", header="NO DGRS FOUND")
                raise ConfigError(f"Therefore, we will exit here because anvi'o has found no DGRs in your data, "
                                        "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                                        "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")

            # Parse the standard BLAST output
            tree = ET.parse(self.blast_output)
            root = tree.getroot()

            # Process the BLAST output normally
            self.process_single(root, bin_name=None, max_percent_identity=max_percent_identity)

            return self.mismatch_hits



    def process_single(self, root, bin_name, max_percent_identity):
        """
        Process a single bin's BLAST XML data and store mismatch hits.
        If bin_name is None, it means we're processing a single BLAST output.
        """

        for iteration in root.findall(".//Iteration"):
            for hit in iteration.findall(".//Hit"):
                for hsp in hit.findall(".//Hsp"):
                    identical_positions = int(hsp.find('Hsp_identity').text)
                    alignment_length = int(hsp.find('Hsp_align-len').text)
                    percentage_identity = (identical_positions / alignment_length) * 100

                    if percentage_identity < max_percent_identity:
                        section_id = iteration.find('Iteration_query-def').text
                        hsp_num = hsp.find('Hsp_num').text

                        hit_identity = '_'.join([section_id, f'_BLAST_hsp_is_{hsp_num}'])
                        pattern = r"start_bp(\d+)_end_bp(\d+)"

                        match = re.search(pattern, section_id)
                        query_start_position = int(match.group(1))
                        query_end_position = int(match.group(2))

                        qseq = str(hsp.find('Hsp_qseq').text)
                        hseq = str(hsp.find('Hsp_hseq').text)
                        midline = str(hsp.find('Hsp_midline').text)

                        has_repeat = False
                        # here we use pytrf for removing tandem and not exact tandem repeats, this is a python compiled version of tandem repeats finder:
                        # the idea is that there are a lot of false positives found in metagenomic samples that are found with `anvi-report-dgrs` due to
                        # the nature of the two sequences being very similar and we only want those that are similar due to the TR/VR constraints and not
                        # because they are repeated sequences. The first loop uses the approximate repeat finder to search for imperfect tandem repeats
                        # that are above a threshold number of times repeated, by default this is 10.
                        #
                        # The second loop loops for sequential tandem repeats and has the parameter to catch any mono or di motif (so single base, or pair of bases)
                        # that repeats 6 times after each other and default settings for tri=5, tetra=4, penta=4, hexa=4. If any are caught the sequence is removed.
                        #
                        # The third loop for approximate repeat finder works by looking for any repeated motif that is over 4 and under 10 bp long that has a minimum
                        # identity fo up to 1 bp difference. then works out a coverage of how much of the sequence the repeat covers and if it is above the user defined
                        # value or 80% as default then the sequence is removed.
                        # The fourth loop is there because the pytrf code has the quirk that the if you change the min motif to smaller it misses the larger motif repeats,
                        # this loop also has the seed len as 9 not 10 which is the default to remove some more repeats

                        for seq in [hseq, qseq]:
                            #remove - in sequences to remove repeats properly
                            seq = seq.replace("-", "")
                            for atr1 in pytrf.ATRFinder('name', seq):
                                if int(atr1.repeat) > self.numb_imperfect_tandem_repeats:
                                    has_repeat = True

                            #look for tandem homopolymers and 2 base short tandem repeats that are over 4 (so occur 5 times) times in the sequence
                            for ssr in pytrf.STRFinder('name', seq, mono=5, di=5):
                                #if ((len(ssr.motif) == 1) or (len(ssr.motif) == 2)) and ssr.repeat > 6:
                                if ssr:
                                    has_repeat = True

                            #look for approximate tandem repeats that in the VR, using a coverage value of the motif length times by the number of repeats
                            # divided by the sequence length
                            for atr in pytrf.ATRFinder('name', seq, min_motif=4, max_motif=10, min_seedrep=2, min_identity=70):
                                coverage = (len(atr.motif)*atr.repeat) / len(seq)
                                if str(coverage) > self.repeat_motif_coverage:
                                    has_repeat = True

                            #look for approximate tandem repeats that in the VR, using a coverage value of the motif length times by the number of repeats
                            # divided by the sequence length
                            # Need to do this for shorter motifs too because pytrf misses them otherwise
                            for atr2 in pytrf.ATRFinder('name', seq, min_motif=3, max_motif=10, min_seedrep=2, min_seedlen=9, min_identity=70):
                                coverage = (len(atr2.motif)*atr2.repeat) / len(seq)
                                if coverage > self.repeat_motif_coverage:
                                    has_repeat = True

                        if has_repeat == True:
                            #breakout of loop (don't add as dgr)
                            continue

                        subject_genome_start_position = min(
                            [int(hsp.find('Hsp_hit-from').text) - 1, int(hsp.find('Hsp_hit-to').text)])
                        subject_genome_end_position = max(
                            [int(hsp.find('Hsp_hit-from').text) - 1, int(hsp.find('Hsp_hit-to').text)])
                        query_genome_start_position = query_start_position + min(
                            [int(hsp.find('Hsp_query-from').text) - 1, int(hsp.find('Hsp_query-to').text)])
                        query_genome_end_position = query_start_position + max(
                            [int(hsp.find('Hsp_query-from').text) - 1, int(hsp.find('Hsp_query-to').text)])
                        query_frame = str(hsp.find('Hsp_query-frame').text)
                        subject_frame = str(hsp.find('Hsp_hit-frame').text)

                        query_mismatch_positions = []
                        all_possible_characters = set(qseq + hseq)

                        query_mismatch_counts = {char: 0 for char in all_possible_characters}
                        subject_mismatch_counts = {char: 0 for char in all_possible_characters}

                        chars_to_skip = []
                        if self.skip_Ns:
                            chars_to_skip.append('N')
                        if self.skip_dashes:
                            chars_to_skip.append('-')

                        for idx in range(len(qseq)):
                            if qseq[idx] in chars_to_skip or hseq[idx] in chars_to_skip:
                                continue
                            if qseq[idx] != hseq[idx]:
                                query_mismatch_counts[qseq[idx]] += 1
                                query_mismatch_positions.append(idx)
                                subject_mismatch_counts[hseq[idx]] += 1

                        query_contig = section_id.split('_section', 1)[0]
                        subject_contig = hit.find('Hit_def').text

                        self.mismatch_hits[hit_identity] = {
                            'bin': bin_name if bin_name else "single",
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



    def add_new_DGR(self, DGR_number, bin, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig, base,
                    is_reverse_complement, TR_frame, VR_sequence, VR_frame, subject_genome_start_position, subject_genome_end_position, subject_contig, midline,
                    percentage_of_mismatches, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative, snv_VR_positions,
                    numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs):
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
        self.DGRs_found_dict[DGR_key] = {}
        # TR stuff
        self.DGRs_found_dict[DGR_key]['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['base'] = base
        self.DGRs_found_dict[DGR_key]['TR_bin'] = bin
        self.DGRs_found_dict[DGR_key]['TR_start_position'] = subject_genome_start_position
        self.DGRs_found_dict[DGR_key]['TR_end_position'] = subject_genome_end_position
        self.DGRs_found_dict[DGR_key]['TR_contig'] = subject_contig
        self.DGRs_found_dict[DGR_key]['TR_sequence_found'] = 'subject'
        #self.DGRs_found_dict[DGR_key]['next_VR_id'] = 2  # next available VR index so we dont overwrite VR ids when we delete them if there are overlapping ones in update_existing_DGR

        # VR stuff
        self.DGRs_found_dict[DGR_key]['VRs'] = {'VR_001':{}}
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_sequence'] = VR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['midline'] = midline
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['percentage_of_mismatches'] = percentage_of_mismatches
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_frame'] = VR_frame
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_bin'] = bin
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_start_position'] = query_genome_start_position
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_end_position'] =   query_genome_end_position
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_contig'] = query_contig
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_sequence_found'] = 'query'
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['DGR_looks_snv_false'] = DGR_looks_snv_false
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['snv_at_3_codon_over_a_third'] = snv_at_3_codon_over_a_third
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['VR_TR_mismatch_positions'] = mismatch_pos_contig_relative
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['snv_VR_positions'] = snv_VR_positions
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['numb_of_snv_in_matches_not_mutagen_base']= numb_of_snv_in_matches_not_mutagen_base
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['numb_of_mismatches']= numb_of_mismatches
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['numb_of_SNVs']= numb_of_SNVs

        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_start_position'] = subject_genome_start_position
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_end_position'] = subject_genome_end_position
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_reverse_complement'] = is_reverse_complement
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_frame'] = TR_frame



    def update_existing_DGR(self, existing_DGR_key, bin, TR_frame, VR_sequence, VR_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement,
                            query_genome_start_position, query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                            subject_contig, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative, snv_VR_positions,
                            numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs):
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
        if existing_DGR_key not in self.DGRs_found_dict:
            raise KeyError(f"Existing DGR key {existing_DGR_key} not found in DGRs_found_dict")

        #keep a list of the VRs that are overlapping and therefore need deleting
        #VRs_to_delete=[]
        #check if the VR overlaps with another one present and only keep the longest one.
        #for vr_key, existing_vr in list(self.DGRs_found_dict[existing_DGR_key]["VRs"].items()):
            #print(f"{existing_vr}")
            #existing_vr_start = existing_vr["VR_start_position"]
            #existing_vr_end = existing_vr["VR_end_position"]


            #if existing_vr['VR_contig'] == query_contig and self.range_overlapping(query_genome_start_position,
                                                                                    #query_genome_end_position,
                                                                                    #existing_vr["VR_start_position"],
                                                                                    #existing_vr["VR_end_position"]):

                #new_length = query_genome_end_position - query_genome_start_position
                #existing_length = existing_vr_end - existing_vr_start

                #if new_length > existing_length:
                    #only add the longest VR and remove existing vr and add new one
                    #VRs_to_delete.append(vr_key)
                    #self.run.warning(f"Heads up. An existing VR is being discarded because it has a range that overlaps with another VR in the same DGR. The new VR is positioned here: {query_genome_start_position}:{query_genome_end_position} "
                                    #f"the existing VR is here: {existing_vr_start}:{existing_vr_end} both on this contig: {existing_vr['VR_contig']}. This could mean that the DGR is less likely to be a true positive, or that it was found twice! Maybe go check "
                                    #"your read recruitment for this contig, which is always a good idea", header="Existing VR is being discarded because it is shorter")
                #else:
                    #self.run.warning(f"Heads up. The new VR is being discarded because it has a range that overlaps with another VR in the same DGR. The new VR is positioned here: {query_genome_start_position}:{query_genome_end_position} "
                                    #f"the existing VR is here: {existing_vr_start}:{existing_vr_end} both on this contig: {existing_vr['VR_contig']}. This could mean that the DGR is less likely to be a true positive, or that it was found twice! Maybe go check "
                                    #"your read recruitment for this contig, which is always a good idea", header="New VR is being discarded because it is shorter")
                    #return skips adding the shorter vr.
                    #return
        # Now delete the old, overlapping VRs safely
        #for vr_key in VRs_to_delete:
            #del self.DGRs_found_dict[existing_DGR_key]['VRs'][vr_key]

        # VR info
        #vr_id = self.DGRs_found_dict[existing_DGR_key].get('next_VR_id', 2)
        #new_VR_key = f'VR_{vr_id:03d}'
        #self.DGRs_found_dict[existing_DGR_key]['next_VR_id'] = vr_id + 1
        num_VR = len(self.DGRs_found_dict[existing_DGR_key]['VRs']) + 1
        new_VR_key = f'VR_{num_VR:03d}'
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key] = {}
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_sequence'] = VR_sequence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_bin'] = bin
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['midline'] = midline
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['percentage_of_mismatches'] = percentage_of_mismatches
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_frame'] = VR_frame
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_start_position'] = query_genome_start_position
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_end_position'] = query_genome_end_position
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_contig'] = query_contig
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_sequence_found'] = 'query'

        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_reverse_complement'] = is_reverse_complement
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_frame'] = TR_frame
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_start_position'] = subject_genome_start_position
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['TR_end_position'] = subject_genome_end_position
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['DGR_looks_snv_false'] = DGR_looks_snv_false
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['snv_at_3_codon_over_a_third'] = snv_at_3_codon_over_a_third
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['VR_TR_mismatch_positions'] = mismatch_pos_contig_relative
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['snv_VR_positions'] = snv_VR_positions
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['numb_of_snv_in_matches_not_mutagen_base']= numb_of_snv_in_matches_not_mutagen_base
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['numb_of_mismatches']= numb_of_mismatches
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['numb_of_SNVs']= numb_of_SNVs

        self.DGRs_found_dict[existing_DGR_key]['TR_start_position'] = min(subject_genome_start_position, self.DGRs_found_dict[existing_DGR_key]['TR_start_position'])
        self.DGRs_found_dict[existing_DGR_key]['TR_end_position'] = max(subject_genome_end_position, self.DGRs_found_dict[existing_DGR_key]['TR_end_position'])



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

        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        #Sort pandas data-frame of SNVs by contig name and then by position of SNV within contig
        self.snv_panda_bin = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
        self.snv_panda_bin['contig_name'] = self.snv_panda_bin.split_name.str.split('_split_').str[0]
        profile_db.disconnect()

        if self.only_a_bases:
                #This is here so that every potential VR doesn't get a new warning and clog up the terminal
                self.run.warning("Just a note to say that we are only looking for DGRs that have A bases as their site of mutagenesis.", header="Searching for only A mutagenesis based DGRs")

        if self.collections_mode:
            hits_item = self.merged_mismatch_hits
        else:
            hits_item = self.mismatch_hits

        for sequence_component, hit_list in hits_item.items():
            # Ensure hit_list is always a list
            if not isinstance(hit_list, list):
                hit_list = [hit_list]  # Wrap non-list items in a list

            for hit_data in hit_list:
                bin = hit_data['bin']
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
                for letter, count in subject_mismatch_counts.items():
                    percentage_of_mismatches = (count / mismatch_length_bp)
                    if (percentage_of_mismatches > self.percentage_mismatch) and (mismatch_length_bp > self.number_of_mismatches):

                        # if the letter is T, then we assume that it is an A base and we reverse complement EVERYTHING
                        if letter == 'T':
                            TR_sequence = str(subject_sequence.reverse_complement())
                            VR_sequence = str(query_sequence.reverse_complement())
                            midline = ''.join(reversed(original_midline))
                            base = 'A'
                            is_reverse_complement = True
                            query_frame = query_frame * -1
                            subject_frame = subject_frame * -1
                        else:
                            TR_sequence = str(subject_sequence)
                            VR_sequence = str(query_sequence)
                            midline = original_midline
                            base = letter
                            is_reverse_complement = False


                        #### SNV logic based removal of DGRs.
                        #first remove dgrs that have more than 34% of the SNVs in the VR range to be coming from the 3rd codon position
                        # - because this is not very DGR like and means that there is more non-specific read mapping happening nothing biological
                        # Secondly, we want to check that the SNVs are coming from the sites of matches in the mutagenesis base (usually A), and
                        # all of the mismatching bases (remember we have already cut out those that have more than 20% of the mismatches coming
                        # from the non-mutagenesis site). This means that we filter out those that have SNVs at the site of matches in the non-mutagenesis base (usually C,G,T).
                        # Currently the code (30.04.2025) adds a column to the DGR_looks_snv_false if there is ONE single SNV at a non-mutagenesis match site.
                        # Going to look and create a threshold - for populations etc

                        #subset snv df by query contig (vr contig) and VR range
                        matching_snv_rows = self.snv_panda_bin[(self.snv_panda_bin['contig_name'] == query_contig) &
                        (self.snv_panda_bin['pos_in_contig'].between(query_genome_start_position, query_genome_end_position))]

                        #make snv vr positions a list so we can print it in the output
                        snv_VR_positions = sorted(set(matching_snv_rows['pos_in_contig'].to_list()))

                        #currently position of mismatches is relative to the VR needs to be relative to the contig so that they match the vr snv ones
                        mismatch_pos_contig_relative = [x + query_genome_start_position for x in position]

                        if matching_snv_rows.empty:
                            #skip because you are useless to us
                            #i.e. there are no SNVs in the VR)
                            continue
                        else:
                            # first check majority of SNVs come from 1&2 codon pos
                            # get the counts of the snv codon positions
                            codon_pos_sorted = matching_snv_rows['base_pos_in_codon'].value_counts()
                            count_1 = codon_pos_sorted.get(1, 0)
                            count_2 = codon_pos_sorted.get(2, 0)
                            count_3 = codon_pos_sorted.get(3, 0)

                            total = count_1 + count_2 + count_3
                            percent_3 = count_3 / total * 100

                            # Apply threshold of 66% because we want less than a third of snvs to be at the third codon position
                            is_3_over_a_third = percent_3 > (self.snv_codon_position * 100)
                            if is_3_over_a_third:
                                self.run.info_single(f"3rd codon position is over a third of the total VR SNVs, so we remove you. {percent_3}")
                                snv_at_3_codon_over_a_third = True
                                self.run.warning("Skipping candidate DGR due to SNV filters. Specifically, in this case the candidate DGR has a high "
                                                "likelihood of being a false positive due to the fact that there are a high proportion of SNVs that are coming "
                                                f"from the third codon base position which is unexpected behaviour. This percentage of SNVs comes from the third codon position: "
                                                f"{percent_3} are in the third codon position of the VR. The percentage cut off for third codon position SNVs is 33% due to the "
                                                "likelihood of a third of the SNVs in an ORF coming from the third base position. If you would like to change this use the '--snv-codon-position' parameter "
                                                "to give it another value.", header="WARNING: DGR REMOVED", lc='yellow')
                            else:
                                snv_at_3_codon_over_a_third = False

                            #look if matches of VR have SNVs (exclude mismatch positions and matches of base of mutagenesis in VR (normally A).)
                            #letter to skip is mutagenesis base (usually A), if mismatch majority A IF rev comp then T, else make mismatch majority use base
                            letter_to_skip = base.upper()

                            if is_reverse_complement:
                                if letter_to_skip == "A":
                                    letter_to_skip = "T"
                                elif letter_to_skip == "T":
                                    letter_to_skip = "A"
                                elif letter_to_skip == "C":
                                    letter_to_skip = "G"
                                elif letter_to_skip == "G":
                                    letter_to_skip = "C"
                                elif letter_to_skip == "N":
                                    letter_to_skip = "N"

                            # subset VR snv df, by matches in VR and TR *AND* reference in VR being mutagenesis base (usually A)
                            snv_in_matches_not_mutagen_base = matching_snv_rows[~matching_snv_rows['pos_in_contig'].isin(mismatch_pos_contig_relative) & (matching_snv_rows['reference'] != letter_to_skip)]
                            #this needs to be a set of the posisiont in contig so that there are not multiple reported
                            numb_of_snv_in_matches_not_mutagen_base = len(set(snv_in_matches_not_mutagen_base['pos_in_contig']))
                            numb_of_mismatches = len(position)
                            numb_of_SNVs = len(snv_VR_positions)

                            # Report proportion of non-mutagenesis SNVs vs all SNVs
                            if numb_of_SNVs > 0:
                                prop_non_mutagen_snv = numb_of_snv_in_matches_not_mutagen_base/numb_of_SNVs

                            # Flag to store final DGR-like decision based on SNVs
                            DGR_looks_snv_false = False

                            # Apply filters based on SNV count and thresholds
                            if self.snv_matching_proportion:
                                max_allowed_bad_snv_fraction = self.snv_matching_proportion
                            else:
                                if numb_of_SNVs >= 30:
                                    max_allowed_bad_snv_fraction = 0.30  # 30%
                                elif numb_of_SNVs < 30:
                                    max_allowed_bad_snv_fraction = 0.25  # 25%

                                # Evaluate
                                if prop_non_mutagen_snv >= max_allowed_bad_snv_fraction:
                                    DGR_looks_snv_false = True
                                    self.run.warning("Skipping candidate DGR due to SNV filters. Specifically, in this case the candidate DGR has a high "
                                                "likelihood of being a false positive due to the fact that there are a high proportion of SNVs that are coming "
                                                f"from the non mismatching or mutagenesis bases which is where the SNVs are expoected to be. There are this many SNVs total: {numb_of_SNVs} of which "
                                                f"{numb_of_snv_in_matches_not_mutagen_base} are in matching positions of the TR and VR (the proportion is therefore: {prop_non_mutagen_snv:.2%}). "
                                                "The cut off for these SNVs is proportional to the total number of SNVs in the VR if there are >30 than 30% SNVs in matching positions are allowed, "
                                                "if less than 30 SNVs than 25% are allowed, this is by default. If you think that this is incorrect please change the '--snv-matching-proportion' parameter "
                                                "to give it a blanket value, this is what we found to be most effective based on our short read metagenome testing.", header="WARNING: DGR REMOVED", lc='yellow')

                        #to test for VR diversity of base types in the protein sequence
                        for letter, count in query_mismatch_counts.items():
                            non_zero_bases = sum(1 for count in query_mismatch_counts.values() if count > 0)
                        if not non_zero_bases >= self.min_mismatching_base_types_vr:
                            continue

                        #to test for VR diversity of base types in the sequence
                        # Count the distinct base types in the sequence
                        vr_unique_bases = set(query_sequence) - {"-", "N"}  # Ignore gaps and ambiguous bases if needed

                        # Ensure the sequence has at least the required number of distinct base types
                        if len(vr_unique_bases) <= self.min_base_types_vr:
                            continue

                        #to test for TR diversity of base types in the sequence
                        # Count the distinct base types in the sequence
                        tr_unique_bases = set(subject_sequence) - {"-", "N"}  # Ignore gaps and ambiguous bases if needed

                        # Ensure the sequence has at least the required number of distinct base types
                        if len(tr_unique_bases) <= self.min_base_types_tr:
                            continue

                        if self.only_a_bases:
                            # Filter out bases with count > 0 before checking
                            nonzero_mismatch_bases = [base for base, count in subject_mismatch_counts.items() if count > 0]

                            # If reverse complemented, treat 'T' as 'A'
                            if is_reverse_complement:
                                all_mismatches_are_A = all(base in ('A', 'T') for base in nonzero_mismatch_bases)
                            else:
                                all_mismatches_are_A = all(base == 'A' for base in nonzero_mismatch_bases)

                            # Skip if any mismatching base is not valid
                            if not all_mismatches_are_A:
                                continue

                        # Here we dont add VR candidates based on SNV parameters.
                        # Skip DGR if flagged due to SNV-based filters
                        if DGR_looks_snv_false or snv_at_3_codon_over_a_third:
                            continue

                        #need to check if the new TR you're looping through exists in the DGR_found_dict, see if position overlap
                        if not self.DGRs_found_dict:
                            # add first DGR
                            num_DGR += 1
                            self.run.warning(f"Adding new DGR {num_DGR} in the bin: {bin}", header="NEW DGR", lc='yellow')
                            self.add_new_DGR(num_DGR, bin, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                        base, is_reverse_complement, query_frame, VR_sequence, subject_frame, subject_genome_start_position, subject_genome_end_position,
                                        subject_contig, midline, percentage_of_mismatches, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative,
                                        snv_VR_positions, numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs)
                        else:
                            was_added = False
                            for dgr in self.DGRs_found_dict:
                                if self.DGRs_found_dict[dgr]['TR_contig'] == subject_contig and self.range_overlapping(subject_genome_start_position,
                                                                                                                subject_genome_end_position,
                                                                                                                self.DGRs_found_dict[dgr]['TR_start_position'],
                                                                                                                self.DGRs_found_dict[dgr]['TR_end_position']):
                                    was_added = True
                                    #TODO can rename consensus_TR
                                    self.update_existing_DGR(dgr, bin, query_frame, VR_sequence, subject_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement, query_genome_start_position,
                                                    query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                                                    subject_contig, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative, snv_VR_positions,
                                                    numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs)
                                    break
                            if not was_added:
                                # add new TR and its first VR
                                self.run.warning(f"Adding new DGR {num_DGR} in the bin: {bin}", header="NEW DGR", lc='yellow')
                                num_DGR += 1
                                self.add_new_DGR(num_DGR, bin, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                        base, is_reverse_complement, query_frame, VR_sequence, subject_frame, subject_genome_start_position, subject_genome_end_position,
                                        subject_contig, midline, percentage_of_mismatches, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative,
                                        snv_VR_positions, numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs)

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
        if len(self.DGRs_found_dict) == 0:
            return

        self.run.info_single(f"Computing the Genes the Variable Regions occur in and creating a '{self.output_directory}_DGR_genes_found.csv'.")

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

                        self.vr_gene_info[dgr][vr] = gene_call
                        break

        #TODO: MAKE into WRITE DGR_genes_found write function
        # define output path
        output_directory_path = self.output_directory
        output_path_for_genes_found = os.path.join(output_directory_path, f"{self.output_directory}_DGR_genes_found.csv")

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
        dgrs_dict = self.DGRs_found_dict

        if len(dgrs_dict) == 0:
            return

        self.run.info_single("Computing the closest HMMs to the Template Regions and printing them in your output csv.")
        self.run.info_single('\n')

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
                self.run.warning(f"No HMM Reverse Transcriptase was found for {DGR_id}. This does not affect the "
                                "outcome of the tool, other than meaning that there is no RT found on the same contig "
                                "as the Template Region meaning it might be a good idea to manually inspect your data "
                                "to see if any RT's are present.", header="NO REVERSE TRANSCRIPTASE HMMs")

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
        output_directory_path = self.output_directory

        dgrs_dict = self.DGRs_found_dict
        output_path_dgrs = os.path.join(output_directory_path, f"{self.output_directory}_DGRs_found.csv")
        headers = [
            "DGR", "VR", "VR_contig", "VR_frame", "VR_sequence", "Midline",
            "VR_start_position", "VR_end_position", "VR_bin", "Mismatch %",
            "TR_contig", "TR_frame", "TR_sequence", "Base", "Reverse Complement",
            "TR_start_position", "TR_end_position", "TR_bin", "TR_in_gene", "HMM_source",
            "distance_to_HMM", "HMM_gene_name", "HMM_direction", "HMM_start",
            "HMM_stop", "HMM_gene_callers_id", "DGR_looks_snv_false", "snv_at_3_codon_over_a_third",
            "numb_of_snv_in_matches_not_mutagen_base", "numb_of_mismatches", "numb_of_SNVs",
            "VR_TR_mismatch_positions", "snv_VR_positions"
        ]

        # Check if either dictionary is empty or lacks meaningful keys
        if not any(dgrs_dict.values()):
            self.run.warning("No DGRS were found so no output file will be written :( However, you can go "
                            "back and tinker with the parameters of this tool if you believe this should not "
                            "be the case. Anvi'o wishes you a nice day :)", header="NO DIVERSITY-GENERATING RETROELEMENTS")
            return

        # Open the CSV file and write headers and rows
        with open(output_path_dgrs, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)

            for dgr, tr in dgrs_dict.items():
                for vr, vr_data in tr['VRs'].items():
                    # Populate csv_row for DGRs_found.csv format
                    csv_row = [
                        dgr, vr, vr_data['VR_contig'], vr_data.get('VR_frame', 'N/A'),
                        vr_data['VR_sequence'], vr_data['midline'], vr_data['VR_start_position'],
                        vr_data['VR_end_position'], vr_data.get('VR_bin'),
                        vr_data['percentage_of_mismatches'], tr['TR_contig'],
                        vr_data.get('TR_frame', 'N/A'), vr_data['TR_sequence'],
                        tr['base'], vr_data['TR_reverse_complement'],
                        vr_data['TR_start_position'], vr_data['TR_end_position'],
                        tr.get('TR_bin', 'N/A'),tr.get('TR_in_gene', 'N/A'), tr['HMM_source'], tr["distance_to_HMM"],
                        tr["HMM_gene_name"], tr["HMM_direction"], tr["HMM_start"],
                        tr["HMM_stop"], tr["HMM_gene_callers_id"], vr_data["DGR_looks_snv_false"],
                        vr_data["snv_at_3_codon_over_a_third"], vr_data["numb_of_snv_in_matches_not_mutagen_base"],
                        vr_data["numb_of_mismatches"], vr_data["numb_of_SNVs"], vr_data["VR_TR_mismatch_positions"], vr_data["snv_VR_positions"]
                    ]
                    csv_writer.writerow(csv_row)
        return



    def recover_genomic_context_surrounding_dgrs(self):
        """Learn about what surrounds the variable region sites of each found DGR"""

        # in which we will store the genomic context that surrounds dgrs for downstream fun
        self.genomic_context_surrounding_dgrs = {}

        dgrs_dict = self.DGRs_found_dict

        if len(dgrs_dict) == 0:
            return

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
                            "for the genes. PITY.")
            return

        self.progress.new('Recovering genomic context surrounding the DGRs', progress_total_items=len(self.DGRs_found_dict))
        self.progress.update('...')

        # now we will go through each dgr to populate `self.genomic_context_surrounding_dgrs`
        # with gene calls and functions
        gene_calls_per_TR_contig = {}
        gene_calls_per_VR_contig = {}
        trs_with_no_gene_calls_around = set([])
        vrs_with_no_gene_calls_around = set([])


        for dgr_key, dgr_data in dgrs_dict.items():
            dgr_id = dgr_key
            self.progress.update(f"{dgr_id}", increment=True)

            TR_contig_name = dgr_data.get('TR_contig')
            TR_start = dgr_data.get('TR_start_position')
            TR_end = dgr_data.get('TR_end_position')

            # Initialize TR context
            tr_context_genes = {}

            is_tr_in_gene = False

            if TR_contig_name not in gene_calls_per_TR_contig:
                where_clause = f'''contig="{TR_contig_name}" and source="{self.gene_caller_to_consider_in_context}"'''
                gene_calls_per_TR_contig[TR_contig_name] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

            gene_calls_in_TR_contig = gene_calls_per_TR_contig[TR_contig_name]

            if not len(gene_calls_in_TR_contig):
                trs_with_no_gene_calls_around.add(dgr_id)
                self.run.info_single(f'No gene calls found around TR {dgr_id}')
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

                # Store the gene call in the dictionary using gene_callers_id as the key
                tr_context_genes[gene_callers_id] = gene_call

                # Check if the TR is inside this gene
                if TR_start >= gene_call['start'] and TR_end <= gene_call['stop']:
                    is_tr_in_gene = True

            # After processing all the gene calls for the TR, add the result to dgrs_dict
            dgrs_dict[dgr_id]['TR_in_gene'] = is_tr_in_gene

            self.genomic_context_surrounding_dgrs[dgr_id] = copy.deepcopy(tr_context_genes)

            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key

                self.progress.update(f"{dgr_id} , {vr_id}", increment=True)

                VR_contig = vr_data.get('VR_contig')
                VR_start = vr_data.get('VR_start_position')
                VR_end = vr_data.get('VR_end_position')

                # Initialize VR context
                vr_context_genes = {}

                if VR_contig not in gene_calls_per_VR_contig:
                    where_clause = f'''contig="{VR_contig}" and source="{self.gene_caller_to_consider_in_context}"'''
                    gene_calls_per_VR_contig[VR_contig] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

                gene_calls_in_VR_contig = gene_calls_per_VR_contig[VR_contig]

                if not len(gene_calls_in_VR_contig):
                    vrs_with_no_gene_calls_around.add(vr_id)
                    self.run.info_single(f'No gene calls found around DGR {dgr_id} VR {vr_id}')
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

                    # Store the gene call in the dictionary using gene_callers_id as the key
                    vr_context_genes[gene_callers_id] = gene_call

                    self.genomic_context_surrounding_dgrs[dgr_id][vr_id] = copy.deepcopy(vr_context_genes)

        contigs_db.disconnect()
        self.progress.end()
        self.run.info_single('Completed recovering genomic context surrounding the DGRs')

        self.run.info(f"[Genomic Context] Searched for {PL('DGR', len(dgrs_dict))}",
                    f"Recovered for {PL('TR', len(self.genomic_context_surrounding_dgrs[dgr_id]))}",
                    f"And recovered for {PL('VR', len(self.genomic_context_surrounding_dgrs[dgr_id][vr_id]))}",
                    nl_before=1,
                    lc="yellow")

        if len(trs_with_no_gene_calls_around):
            self.run.warning('No gene calls around the following TRs:', trs_with_no_gene_calls_around, "Here is the list in case you would like to track them down: "
                            f"{', '.join(trs_with_no_gene_calls_around)}.")
        if len(vrs_with_no_gene_calls_around):
            self.run.warning('No gene calls around the following VRs:', vrs_with_no_gene_calls_around, "Here is the list in case you would like to track them down: "f"{', '.join(vrs_with_no_gene_calls_around)}.")

        if not len(self.genomic_context_surrounding_dgrs):
            self.run.warning(f"Even though the tool went through all {PL('DGR', len(dgrs_dict))} "
                            f"it was unable to recover any genomic context for any of them. So your final reports will "
                            f"not include any insights into the surrounding genomic context of all your DGRs ( a little bit sad "
                            f"but otherwise you will be fine).")



    def report_genomic_context_surrounding_dgrs(self):
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
            self.run.warning("No genomic context data available to report on any of the DGRs :(")
            return

        dgrs_dict = self.DGRs_found_dict

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

                # Check if there are surrounding genes for the TR and write them
                if dgr_id in self.genomic_context_surrounding_dgrs:
                    # Iterate over the gene_callers_id and the associated gene_call dictionary
                    for gene_callers_id, gene_call in self.genomic_context_surrounding_dgrs[dgr_id].items():
                        # Ensure the gene_call is a dictionary
                        if isinstance(gene_call, dict):
                            # Write gene information to the output, ensuring we access the correct fields
                            tr_genes_output.write(f"{dgr_id}_TR\tGENE\t%s\n" % '\t'.join([f"{gene_call.get(h, '')}" for h in genes_output_headers]))

                            # If 'functions' is present in the gene_call, write functions to the output
                            if 'functions' in gene_call and gene_call['functions']:
                                for hit in gene_call['functions']:
                                    tr_functions_output.write(f"{dgr_id}_TR\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                            else:
                                # Write placeholder if no functions are found
                                tr_functions_output.write(f"{dgr_id}_TR\t{gene_call.get('gene_callers_id', '')}\t\t\t\n")
                        else:
                            # Log if gene_call is not a dictionary
                            self.run.info_single(f"Unexpected type for gene_call: {gene_call} (expected dict but got {type(gene_call)})")

                # Log information about the reporting files
                self.run.info(f"Reporting file on gene context for {dgr_id} TR", tr_genes_output_path)
                self.run.info(f"Reporting file on functional context for {dgr_id} TR", tr_functions_output_path, nl_after=1)

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
                                    vr_functions_output.write(f"{dgr_id} {vr_id}\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                            else:
                                vr_functions_output.write(f"{dgr_id} {vr_id}\t{gene_call['gene_callers_id']}\t\t\t\n")

                    self.run.info(f'Reporting file on gene context for {dgr_id} {vr_id}', vr_genes_output_path)
                    self.run.info(f'Reporting file on functional context for {dgr_id} {vr_id}', vr_functions_output_path, nl_after=1)



    # Function to get the consensus base
    @staticmethod
    def get_consensus_base(row):
        for nucleotide in nucleotides:
            if row[nucleotide] > 0.5:  # Assuming a threshold for consensus, adjust as necessary
                return nucleotide
        return None



    @staticmethod
    def compute_dgr_variability_profiling_per_vr(self, input_queue, output_queue, samples_dict, primers_dict, output_directory_path, run=run_quiet, progress=progress_quiet, sample_primers_dict=None, use_sample_primers=False):
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
        Sample_primers_dict : dict, optional
            Dictionary of sample-specific primers.
        Use_sample_primers : bool, optional
            Whether to use sample-specific primers.

        Returns
        =======


        """
        while True:
            sample_name = input_queue.get(True)
            self.run.info_single(f"(I am in compute per sample func) Processing sample: {sample_name}")
            if sample_name is None:
                self.run.warning('Sample {sample_name} is none, loop will be broken.')
                break

            self.run.single_info(f"Processing sample: {sample_name}")
            #Extract sample-specific primers if the flag is set
            if use_sample_primers and sample_primers_dict:
                primers_for_sample = sample_primers_dict.get(sample_name, primers_dict)
            else:
                primers_for_sample = primers_dict
            samples_dict_for_sample = {sample_name: samples_dict[sample_name]}

            # setup the args object
            args = argparse.Namespace(samples_dict= samples_dict_for_sample,
                                    primers_dict= primers_dict,
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

                VR_frame = vr_data['VR_frame']
                TR_frame = vr_data['TR_frame']

                #CHECK if VR sequence is not at the start of a contig so you can get the initial primer sequence
                if vr_data['VR_start_position'] >= self.initial_primer_length:
                    vr_primer_region_start = vr_data['VR_start_position'] - self.initial_primer_length
                    vr_primer_region_end = (vr_data['VR_start_position'] -1)
                else:
                    #this will take the start of the contig to create a shorter initial primer length
                    vr_primer_region_start = vr_data['VR_start_position']
                    vr_primer_region_end = (vr_data['VR_start_position'] -1)
                    self.run.warning(f"The primer sequence for this VR {dgr_id}_{vr_id}, is going to fall off the start of the Contig. This is the pesky VR contig: {vr_data['VR_contig']}.")

                contig_sequence = self.contig_sequences[vr_data['VR_contig']]['sequence']
                vr_primer_region = contig_sequence[vr_primer_region_start:vr_primer_region_end + 1]
                #add every primer sequence to the dgrs_dict
                vr_data['vr_primer_region'] = vr_primer_region

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

                if anvio.DEBUG:
                    self.run.info_single(f"After TR Reverse Complementary Check:{dgr_id} "
                        f"TR sequence: {TR_sequence} "
                        f"TR frame {TR_frame}")

                    self.run.info_single(f"After VR Reverse Complementary Check: {vr_id} "
                        f"VR sequence: {VR_sequence} "
                        f"VR frame {VR_frame}")

                    print('\n')

                #make every frame positive so that the initial primer is always on the left and then both the vr and tr are being compared on the same strand.
                #Always make the VR strand +1 so that you can compare the primer to the fasta file by definition
                if vr_data['VR_frame'] == -1:
                    self.run.info_single(f"I am reverse complementing {dgr_id} {vr_id} AND THE TR because both are -1")
                    vr_data['rev_comp_VR_seq'] = utils.rev_comp(vr_data['VR_sequence'])
                    VR_frame = VR_frame * -1
                    vr_data['VR_reverse_comp_for_primer'] = True
                    VR_sequence = str(vr_data['rev_comp_VR_seq'])
                    vr_data['rev_comp_TR_seq'] = utils.rev_comp(vr_data['TR_sequence'])
                    TR_frame = TR_frame * -1
                    vr_data['TR_reverse_comp_for_primer'] = True
                    TR_sequence = str(vr_data['rev_comp_TR_seq'])
                    if anvio.DEBUG:
                        self.run.info_single(f"AFTER VR + TR = -1:\n{vr_id} VR sequence: {VR_sequence}\n VR frame {VR_frame}")
                        self.run.info_single(f"AFTER VR + TR = -1:\n{dgr_id} TR sequence: {TR_sequence}\n  TR frame {TR_frame}")
                        print('\n')

                if anvio.DEBUG:
                    self.run.info_single(f"FINAL:\n{vr_id} VR sequence: {VR_sequence}\n VR frame {VR_frame}")
                    self.run.info_single(f"FINAL:\n{dgr_id} TR sequence: {TR_sequence}\n  TR frame {TR_frame}")
                    print('\n')

                #check the TR and VR sequence are the same length
                if len(TR_sequence) == len(VR_sequence):
                    profile_db = dbops.ProfileDatabase(self.profile_db_path)
                    snvs_table = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
                    profile_db.disconnect()

                    # Extract the contig name from split_name
                    snvs_table['contig_name'] = snvs_table['split_name'].str.split('_split_', expand=True)[0]

                    vr_positions = set(range(vr_start, vr_end))

                    # Filter snvs_table to include only rows with the matching VR_contig
                    filtered_snvs_table = snvs_table[(snvs_table['contig_name'] == vr_contig) & (snvs_table['pos_in_contig'].isin(vr_positions))]

                    # Convert the positions to a set for faster lookups
                    primer_snv_positions = set(filtered_snvs_table['pos_in_contig'])

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

                elif len(TR_sequence) != len(VR_sequence):
                    self.run.warning(f"Thats weird! The {vr_id} does not have the same length as the {dgr_id}'s TR :( so you can't create an anchor primer sequence")

        ###########################################
        # UPDATED DGRs dict with Primer Sequences #
        ###########################################

        for dgr_id, dgr_data in dgrs_dict.items():
            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key
                primers_dict[dgr_id + '_' + vr_id + '_Primer'] = {'initial_primer_sequence': vr_data['vr_primer_region'],
                                                                'vr_anchor_primer': vr_data['vr_anchor_primer']}
        return primers_dict



    def print_primers_dict_to_csv(self, primers_dict):
        self.run.info_single("Storing the primers used to calculate VR diversity in 'DGR_Primers_used_for_VR_diversity.csv'.")
        output_directory_path = self.output_directory
        output_path= os.path.join(output_directory_path, "DGR_Primers_used_for_VR_diversity.csv")

        if self.skip_primer_variability:
            # Define the header for the CSV file
            csv_header = ['Primer_ID', 'Initial_Primer', 'Anchor_Primer', 'Whole_Primer']

            # Open the CSV file in write mode
            with open(output_path, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(csv_header)  # Write the header row
                # Iterate through the dictionary and write each gene's information to the CSV file
                for primer_name, primer_info in primers_dict.items():
                        used_original_primer = primer_info['used_original_primer']
                        initial_primer = primer_info['initial_primer_sequence']
                        anchor_primer = primer_info['vr_anchor_primer']
                        whole_primer = primer_info['primer_sequence']

                        writer.writerow([
                        primer_name,
                        used_original_primer,
                        initial_primer,
                        anchor_primer,
                        whole_primer
                        ])
        else:
            primers_dict = self.sample_primers_dict
            # Define the header for the CSV file
            csv_header = ['Primer_ID', 'Sample_ID', 'No_SNV_Primer', 'Initial_Primer', 'Anchor_Primer', 'Whole_Primer']

            # Open the CSV file in write mode
            with open(output_path, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(csv_header)  # Write the header row

                # Iterate through the dictionary and write each primer's information to the CSV file
                for primer_name, samples in primers_dict.items():
                    for sample_id, primer_info in samples.items():
                        used_original_primer = primer_info['used_original_primer']
                        initial_primer = primer_info['initial_primer_sequence']
                        anchor_primer = primer_info['vr_anchor_primer']
                        whole_primer = primer_info['primer_sequence']

                        writer.writerow([
                        primer_name,
                        sample_id,
                        used_original_primer,
                        initial_primer,
                        anchor_primer,
                        whole_primer
                        ])
        return



    def compute_dgr_variability_profiling(self):
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
        #TODO: Double check if these are already checked in sanity check
        if self.skip_compute_DGR_variability_profiling or not self.raw_r1_r2_reads_are_present:
            return

        dgrs_dict = self.DGRs_found_dict

        if not len(dgrs_dict):
            self.run.info_single("Compute DGR variability profile function speaking: There are no DGRs to "
                                "compute in-sample variability :/", mc="red")

        sample_names = list(self.samples_txt_dict.keys())
        num_samples = len(sample_names)

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

        if not self.skip_primer_variability:
            self.sample_primers_dict = {}
            sample_names = list(self.samples_txt_dict.keys())
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            snvs_table = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name).sort_values(by=['split_name', 'pos_in_contig'])
            profile_db.disconnect()

            # Sanity check for mismatch between samples given and samples in SNV table
            sample_names_given = set(sample_names)
            sample_names_in_snv_table = set(snvs_table['sample_id'])
            samples_missing_in_snv_table = sample_names_given.difference(sample_names_in_snv_table)

            if anvio.DEBUG:
                self.run.info("Samples given", ", ".join(list(sample_names_given)))
                self.run.info("Samples in profile.db's nucleotide variability table", ", ".join(list(sample_names_in_snv_table)))
                self.run.info("Missing samples from profile.db's nucleotide variability table", ", ".join(list(samples_missing_in_snv_table)))

            if sample_names_given == samples_missing_in_snv_table:
                raise ConfigError(f"Anvi'o is not angry, just disappointed :/ You gave 'anvi-report-dgrs' these samples ({list(sample_names_given)}), but you have none of them in your profile.db. "
                                "This is fatal; Anvi'o will now quit. Either recreate your profile.db with the samples you would like to search for the DGR VRs variability, "
                                "or give 'anvi-report-dgrs' the correct samples.")

            # Iterate over the samples
            for sample_name in sample_names:
                if sample_name not in sample_names_in_snv_table:
                    self.run.warning(f"Sample {sample_name} is missing from the SNV table, skipping this sample.")
                    continue  # Skip this sample if it's not in the SNV table

                # Filter SNVs for the current sample
                sample_snvs = snvs_table[snvs_table['sample_id'] == sample_name]

                # Iterate through DGRs and VRs
                for dgr_id, dgr_data in dgrs_dict.items():
                    for vr_key, vr_data in dgr_data['VRs'].items():
                        vr_id = vr_key
                        vr_start = vr_data.get('VR_start_position')
                        vr_end = vr_data.get('VR_end_position')
                        vr_contig = vr_data.get('VR_contig')

                        # Filter SNVs within the primer region for the current VR
                        primer_snvs = sample_snvs[
                            (sample_snvs['split_name'].apply(lambda x: x.split('_split')[0]) == vr_contig) &
                            (sample_snvs['pos_in_contig'] >= vr_start - self.initial_primer_length) &
                            (sample_snvs['pos_in_contig'] < vr_end + 1)
                        ]

                        dgr_vr_key = f'{dgr_id}_{vr_key}'
                        original_primer_key = f'{dgr_id}_{vr_key}_Primer'
                        if dgr_vr_key not in self.sample_primers_dict:
                            self.sample_primers_dict[dgr_vr_key] = {}

                        if anvio.DEBUG:
                            self.run.info_single(f"Processing sample {sample_name} for DGR {dgr_id} VR {vr_id}")
                            print(primers_dict[original_primer_key])
                            print(primers_dict)
                        if not primer_snvs.empty:
                            # Get original sequences separately
                            original_initial_primer = primers_dict[original_primer_key]['initial_primer_sequence']
                            vr_anchor_primer = primers_dict[original_primer_key]['vr_anchor_primer']

                            # Work on a copy of the initial primer sequence only
                            new_initial_primer = list(original_initial_primer)

                            # Vectorized operation to find consensus SNVs and update the primer sequence
                            # TODO: make departure from ref a parameter?
                            consensus_snvs = primer_snvs[primer_snvs['departure_from_reference'] > 0.5].apply(DGR_Finder.get_consensus_base, axis=1)
                            positions_in_primer = (vr_start - self.initial_primer_length) - primer_snvs[primer_snvs['departure_from_reference'] > 0.5]['pos_in_contig'] - 1

                            for position, consensus_base in zip(positions_in_primer, consensus_snvs):
                                if consensus_base and 0 <= position < len(new_initial_primer):
                                    new_initial_primer[position] = consensus_base

                            # Store the initial + anchor separately, combine only later
                            self.sample_primers_dict[dgr_vr_key][sample_name] = {
                                'initial_primer_sequence': ''.join(new_initial_primer),
                                'vr_anchor_primer': vr_anchor_primer,
                                'used_original_primer': False,
                            }

                            self.run.info_single(f"Updated sample {sample_name} primer for {dgr_vr_key}: {''.join(new_initial_primer)}")
                        else:
                            # Use the original primer sequence since no SNVs were found
                            self.sample_primers_dict[dgr_vr_key][sample_name] = {
                                'initial_primer_sequence': original_initial_primer,
                                'vr_anchor_primer': vr_anchor_primer,
                                'used_original_primer': True,
                            }
                            self.run.warning(f"No valid SNVs for primer region in sample {sample_name} for {dgr_vr_key}, skipping sample consensus.")
                            continue

                if anvio.DEBUG:
                    self.run.info_single(f"Sample {sample_name} processed. Sample-specific primers dict: {self.sample_primers_dict}")

                # Update the sample-specific primers dictionary with the new primer sequence
                for dgr_vr_key, samples in self.sample_primers_dict.items():
                    for sample_name, primer_data in samples.items():
                        primer_key = f"{dgr_vr_key}_Primer"
                        vr_anchor_primer = primers_dict[primer_key]['vr_anchor_primer']
                        primer_data['vr_anchor_primer'] = vr_anchor_primer

                        # Combine the initial primer sequence and vr_anchor_primer for each sample
                        initial_primer = primer_data['initial_primer_sequence']
                        vr_anchor_primer = primer_data['vr_anchor_primer']
                        primer_sequence = initial_primer + vr_anchor_primer

                        # Add the combined primer sequence to the sample data
                        primer_data['primer_sequence'] = primer_sequence

                        if anvio.DEBUG:
                            # Print the updated primer sequence for debugging
                            self.run.info_single(f"Sample: {sample_name}, DGR: {dgr_vr_key}, Primer Sequence: {primer_sequence}")

            # Final processing to ensure primer sequence consistency and length restrictions
            for dgr_id, dgr_data in dgrs_dict.items():
                for vr_key, vr_data in dgr_data['VRs'].items():
                    vr_id = vr_key
                    primer_key = f'{dgr_id}_{vr_id}_Primer'
                    # Set the final primer sequence in primers_dict based on initial sequence and variability analysis
                    primers_dict[primer_key]['primer_sequence'] = primers_dict[primer_key]['initial_primer_sequence'] + primers_dict[primer_key]['vr_anchor_primer']

                    # Ensure the primer sequence does not exceed the desired length
                    if len(primers_dict[primer_key]['primer_sequence']) > self.whole_primer_length:
                        print(f"The primer for {sample_name} {dgr_id} {vr_id} is above the desired length. Trimming to {self.whole_primer_length}.")
                        primers_dict[primer_key]['primer_sequence'] = primers_dict[primer_key]['primer_sequence'][:self.whole_primer_length]

            if anvio.DEBUG:
                self.run.info_single(f"Updated sample primers dictionary: {self.sample_primers_dict}")
            # Output the final primers dictionary
            self.run.info_single("Computing the Variable Regions Primers and creating a 'DGR_Primers_used_for_VR_diversity.csv' file.")
            self.print_primers_dict_to_csv(self.sample_primers_dict if not self.skip_primer_variability else primers_dict)

            if not self.skip_primer_variability:
                self.run.info_single("Primer variability analysis is enabled. Using sample-specific primers.")

                # Ensure `sample_primers_dict` is updated and passed during computation
                use_sample_primers = True
            else:
                self.run.info_single("Skipping primer variability analysis. Using default primers.")
                use_sample_primers = False
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
            self.run.info_single(f"starting worker {i}")
            worker = multiprocessing.Process(target=DGR_Finder.compute_dgr_variability_profiling_per_vr,
                                            args=(input_queue,
                                                output_queue,
                                                self.samples_txt_dict,
                                                primers_dict,
                                                primer_folder),

                                            kwargs=({'progress': self.progress if self.num_threads == 1 else progress_quiet,
                                                    'use_sample_primers': use_sample_primers,  # Pass flag based on `self.skip_primer_variability`
                                                    'sample_primers_dict': self.sample_primers_dict if use_sample_primers else None,  # Pass only if needed
                                                }))
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

        ######################
        # END MULTITHREADING #
        ######################



    def process_dgr_data_for_HTML_summary(self):
        """Take everything that is known, turn them into data that can be used from Django templates.

        A lot of ugly/frightening things happening here to prepare coordinates for SVG objects to be displayed
        or store boolean variables to use the Django template engine effectively. IF YOU DON'T LIKE
        IT DON'T LOOK AT IT. IT MIGHT MAKE YOU CRY
        """
        if not len(self.DGRs_found_dict):
            self.run.warning("No DGRs were found so no HTML file will be written :(", header="No HTML OUTPUT")
            return

        # in which we will store all the static HTML output related stuff
        self.summary = {}

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.summary['meta'] = {'summary_type': 'dgrs',
                                'num_dgrs': len(self.DGRs_found_dict), # len(self.dgrs_in_collections) if self.collections_mode else
                                #'num_samples': len(self.profile_db_paths) if self.collections_mode else len(self.collections_given),
                                'output_directory': self.output_directory,
                                'genomic_context_recovered': not self.skip_recovering_genomic_context,
                                # if no function source, it says 'the contigs.db' because it fits with the message
                                # displayed in the final index.html. See the inversion template, line 215
                                # if it works, it works
                                'gene_function_sources': contigs_db.meta['gene_function_sources'] or ['the contigs.db']}
        contigs_db.disconnect()

        self.summary['files'] = {'Putative_DGRs': 'Putative-DGRs.txt'}
        self.summary['dgrs'] = {}

        dgrs_dict = self.DGRs_found_dict

        for dgr_key, dgr_data in dgrs_dict.items():
            # Assuming dgr_key itself is the dgr_id or a dictionary containing it
            dgr_id = dgr_key

            self.summary['dgrs'][dgr_id] = {'dgr_data': copy.deepcopy(dgr_data), 'tr_genes': {}, 'vr_genes': {}}

            # Handle genomic context recovery
            if not self.skip_recovering_genomic_context:
                # Deep copy the genomic context for TR and VR genes
                tr_genes = {gene_id: gene_info
                    for gene_id, gene_info in self.genomic_context_surrounding_dgrs.get(dgr_id, {}).items()
                    if isinstance(gene_id, int)}

                # then we will learn about these so we can transform the coordinates of anything we wish
                # to display in the output
                # Get the start and stop positions of the first and last genes in tr_genes
                genomic_context_start_tr = min(gene['start'] for gene in tr_genes.values()) - 100
                genomic_context_end_tr = max(gene['stop'] for gene in tr_genes.values()) + 100

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

                for gene_id, gene in tr_genes.items():
                    # Transform start and stop coordinates for the gene
                    gene['start_tr_g'] = (gene['start'] - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length
                    gene['stop_tr_g'] = (gene['stop'] - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length

                    if (gene['stop_tr_g'] - gene['start_tr_g']) < default_gene_arrow_width:
                        # If the transformed length of the gene is smaller than the default arrow width
                        gene_arrow_width = gene['stop_tr_g'] - gene['start_tr_g']
                        gene['stop_tr_g'] = gene['start_tr_g']
                        gene['TRW'] = 0
                    else:
                        gene_arrow_width = default_gene_arrow_width
                        gene['TRW'] = (gene['stop_tr_g'] - gene['start_tr_g']) - gene_arrow_width

                    # Determine color based on gene functions first
                    if gene['functions']:
                        gene['has_functions'] = True
                        gene['COLOR'] = '#008000'  # Green for genes with functions
                    else:
                        gene['has_functions'] = False
                        gene['COLOR'] = '#c3c3c3'  # Grey for genes without functions

                    # Check for hmm_id and gene_id match if provided in the dgr_data
                    try:
                        hmm_id = int(dgr_data.get('HMM_gene_callers_id'))
                        current_gene_id = int(gene.get('gene_callers_id'))
                    except (TypeError, ValueError):
                        hmm_id = None
                        current_gene_id = None

                    if hmm_id is not None and current_gene_id is not None and hmm_id == current_gene_id:
                        gene['COLOR'] = '#c366e8'  # Purple if there's a match based on hmm_id and gene_id

                    # Additional transformations for coordinates
                    gene['TRX'] = gene['start_tr_g']
                    gene['TCX'] = (gene['start_tr_g'] + (gene['stop_tr_g'] - gene['start_tr_g']) / 2)
                    gene['TGY'] = gene['TRX'] + gene['TRW'] + gene_arrow_width
                    gene['TGTRANS'] = gene['TRX'] + gene['TRX'] + gene['TRW'] + gene_arrow_width
                    gene['TRX_TRW'] = gene['TRX'] + gene['TRW'] - 0.5

                    # Append the gene to the summary
                    self.summary['dgrs'][dgr_id]['tr_genes'][gene_id] = gene

                for vr_key, vr_data in dgr_data.get('VRs', {}).items():
                    vr_id = vr_key
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'] = self.summary['dgrs'][dgr_id]['dgr_data'].get('VRs', {})

                    # Extract VR genes (string keys like 'VR_001')
                    vr_genes = self.genomic_context_surrounding_dgrs.get(dgr_id, {})[vr_key]

                    # Calculate genomic context for VR genes
                    genomic_context_start_vr = min(
                        vr_genes[gene]['start']
                        for gene in vr_genes
                    ) - 100

                    genomic_context_end_vr = max(
                        vr_genes[gene]['stop']
                        for gene in vr_genes
                    ) + 100

                    #transform coordinates for VR data
                    vr_start = (int(vr_data['VR_start_position']) - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length
                    vr_end  = (int(vr_data['VR_end_position']) - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length

                    # store transformed VR info
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['VX'] = vr_start
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['VW'] = vr_end - vr_start
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['VT'] = vr_start + (vr_end - vr_start) / 2
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['midline'] = vr_data['midline']

                    for gene in vr_genes.values():
                        gene['start_vr_g'] = (gene['start'] - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length
                        gene['stop_vr_g'] = (gene['stop'] - genomic_context_start_vr) / (genomic_context_end_vr - genomic_context_start_vr) * new_context_length

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

                        # Determine color based on gene functions first
                        if gene['functions']:
                            gene['has_functions'] = True
                            gene['COLOR'] = '#008000'  # Green for genes with functions
                        else:
                            gene['has_functions'] = False
                            gene['COLOR'] = '#c3c3c3'  # Grey for genes without functions

                        # Check for hmm_id and gene_id match if provided in the dgr_data
                        try:
                            hmm_id = int(dgr_data.get('HMM_gene_callers_id'))
                            current_gene_id = int(gene.get('gene_callers_id'))
                        except (TypeError, ValueError):
                            hmm_id = None
                            current_gene_id = None

                        if hmm_id is not None and current_gene_id is not None and hmm_id == current_gene_id:
                            gene['COLOR'] = '#c366e8'  # Purple if there's a match based on hmm_id and gene_id


                        gene['VRX'] = gene['start_vr_g']
                        gene['VCX'] = (gene['start_vr_g'] + (gene['stop_vr_g'] - gene['start_vr_g']) / 2)
                        gene['VGY'] = gene['VRX'] + gene['VRW'] + gene_arrow_width
                        gene['VGTRANS'] = gene['VRX'] + gene['VRX'] + gene['VRW'] + gene_arrow_width
                        gene['VRX_VRW'] = gene['VRX'] + gene['VRW'] - 0.5 # <-- minus 0.5 makes the arrow nicely cover the rest of the gene

                    # Append transformed VR genes to the summary
                    self.summary['dgrs'][dgr_id]['vr_genes'][vr_id] = gene

                    # and finally we will store this hot mess in our dictionary
                    self.summary['dgrs'][dgr_id]['tr_genes'] = tr_genes
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['vr_genes'] = vr_genes

                    # also we need the path to the output files
                    # Store output file paths for DGRs and VRs
                    self.summary['files'][dgr_id] = {
                        'tr_genes': os.path.join('PER_DGR', dgr_id, 'SURROUNDING-GENES.txt'),
                        'functions': os.path.join('PER_DGR', dgr_id, 'SURROUNDING-FUNCTIONS.txt'),
                        'vr_genes': os.path.join('PER_DGR', vr_id, 'SURROUNDING-GENES.txt'),
                        'vr_functions': os.path.join('PER_DGR', vr_id, 'SURROUNDING-FUNCTIONS.txt')
                    }

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
            #TODO:MAKE SURE THIS LIST IS COMPLETE AND IN THE RIGHT ORDER
            ##############
            parameters = [
                ("Contig.db", self.contigs_db_path if self.contigs_db_path else None),
                ("Profile.db", self.profile_db_path if self.profile_db_path else None),
                ("Word Size of BLAST", self.word_size if self.word_size else "8"),
                ("Number of Threads for BLAST", self.num_threads if self.num_threads else None),
                ("Skip 'Ns'", self.skip_Ns if self.skip_Ns else "FALSE"),
                ("Skip '-'", self.skip_dashes if self.skip_dashes else "FALSE"),
                ("Number of Mismatches", self.number_of_mismatches if self.number_of_mismatches else "7"),
                ("Percentage of Mismatches", self.percentage_mismatch if self.percentage_mismatch else "0.8"),
                ("Minimum Mismatching Base Types in VR", self.min_base_types_vr if self.min_base_types_vr else "2"),
                ("Minimum Base Types in VR", self.min_base_types_vr if self.min_base_types_vr else "2"),
                ("Minimum Base Types in TR", self.min_base_types_tr if self.min_base_types_tr else "2"),
                ("Temporary Directory", self.temp_dir if self.temp_dir else None),
                ("Distance between SNVs", self.min_dist_bw_snvs if self.min_dist_bw_snvs else "5"),
                ("Variable Buffer Length", self.variable_buffer_length if self.variable_buffer_length else "20"),
                ("Departure from Reference (Percentage)", self.departure_from_reference_percentage if self.departure_from_reference_percentage else "0.1"),
                ("Minimum Range size of High Density SNVs", self.min_range_size if self.min_range_size else "5"),
                ("Gene caller", self.gene_caller_to_consider_in_context if self.gene_caller_to_consider_in_context else "prodigal"),
                ("HMMs Provided to Search through", self.hmm if self.hmm else "Reverse_Transcriptase"),
                ("Discovery mode", self.discovery_mode if self.discovery_mode else "FALSE"),
                ("Collections Mode", self.collections_mode if self.collections_mode else "FALSE"),
                ("Output Directory", self.output_directory if self.output_directory else "default")
            ]

            csv_writer.writerows(parameters)
        return



    def process(self,args):
        """Here we process all of the functions in our class and call upon different functions depending on the inputs used"""
        self.sanity_check()
        self.get_blast_results()
        self.process_blast_results()
        self.filter_for_TR_VR()
        if args.parameter_output:
            self.run.info_single("Writing to Parameters used file.")
            print('\n')
            self.parameter_output_sheet()
        self.get_gene_info()
        self.get_hmm_info()
        self.recover_genomic_context_surrounding_dgrs()
        self.report_genomic_context_surrounding_dgrs()
        self.create_found_tr_vr_csv()
        self.compute_dgr_variability_profiling() # add if statement to this for the metagenomics mode DGRs
        self.process_dgr_data_for_HTML_summary()
        return