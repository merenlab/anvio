# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity-Generating Retroelements"""

import re
import csv
import os
import shutil
import argparse
import copy
import bisect
import numpy as np
import pytantan

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths
import anvio.tables as t
import anvio.fastalib as fastalib

import multiprocess as multiprocessing # type: ignore
import xml.etree.ElementTree as ET

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.summaryhtml import SummaryHTMLOutput
from anvio.sequencefeatures import PrimerSearch
from anvio.constants import nucleotides
from anvio.artifacts.samples_txt import SamplesTxt

from Bio.Seq import Seq # type: ignore
from collections import defaultdict

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

# pre-compiled regex for extracting start/end positions from section IDs (used in BLAST parsing loop)
SECTION_ID_PATTERN = re.compile(r"start_bp(\d+)_end_bp(\d+)")


def execute_blast(query_records, target_sequences, temp_dir, word_size, num_threads=1,
                  query_fasta_filename="blast_query.fasta", target_fasta_filename="blast_target.fasta",
                  blast_output_filename="blast_output.xml"):
    """
    Execute BLASTn search with given query and target sequences.

    This is a standalone function (not a method) that handles the core BLAST execution
    logic. It can be called from both instance methods and static methods, avoiding
    code duplication.

    Parameters
    ==========
    query_records : dict
        Query sequences as {section_id: sequence}.
    target_sequences : dict
        Target sequences as {contig_name: {'sequence': seq}} or {contig_name: sequence}.
    temp_dir : str
        Directory for temporary BLAST files.
    word_size : int
        BLAST word size parameter.
    num_threads : int
        Number of threads for BLAST (default 1).
    query_fasta_filename : str
        Filename for query FASTA file.
    target_fasta_filename : str
        Filename for target FASTA file.
    blast_output_filename : str
        Filename for BLAST output XML file.

    Returns
    =======
    str
        Path to the BLASTn output XML file.

    Raises
    ======
    ConfigError
        If no query or target sequences are provided.
    """
    if not query_records:
        raise ConfigError("No query sequences provided for BLAST search.")

    if not target_sequences:
        raise ConfigError("No target sequences found for BLAST search.")

    query_fasta_path = os.path.join(temp_dir, query_fasta_filename)
    target_fasta_path = os.path.join(temp_dir, target_fasta_filename)
    blast_output_path = os.path.join(temp_dir, blast_output_filename)

    # Write query FASTA
    query_fasta = fastalib.FastaOutput(query_fasta_path)
    for seq_id, seq in query_records.items():
        query_fasta.write_id(seq_id)
        query_fasta.write_seq(seq)
    query_fasta.close()

    # Write target FASTA
    # Handle both dict formats: {'sequence': seq} and just seq
    target_fasta = fastalib.FastaOutput(target_fasta_path)
    for name, sequence in target_sequences.items():
        if isinstance(sequence, dict):
            seq = sequence['sequence']
        else:
            seq = sequence
        target_fasta.write_id(name)
        target_fasta.write_seq(seq)
    target_fasta.close()

    # Run BLAST
    blast = BLAST(query_fasta_path, target_fasta=target_fasta_path, search_program='blastn',
                  output_file=blast_output_path, additional_params='-dust no', num_threads=num_threads)
    blast.evalue = 10
    blast.makedb(dbtype='nucl')
    blast.blast(outputfmt='5', word_size=word_size)

    return blast_output_path


class DGR_Finder:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x):
            return args.__dict__[x] if x in args.__dict__ else None

        self.contigs_db_path = A('contigs_db')
        self.profile_db_path= A('profile_db')
        self.word_size = A('word_size')
        self.skip_Ns = A('skip_Ns')
        self.skip_dashes = A('skip_dashes')
        self.number_of_mismatches = A('number_of_mismatches')
        self.initial_mismatch_bias_threshold = A('initial_mismatch_bias_threshold') or 0.6
        self.max_non_dominant = A('max_non_dominant') or 1
        self.minimum_vr_length = A('minimum_vr_length') or 50
        self.min_mismatching_base_types_vr = A('min_mismatching_base_types_vr') or 2
        self.min_base_types_vr = A('min_base_types_vr') or 2
        self.min_base_types_tr = A('min_base_types_tr') or 2
        self.allow_any_base = A('allow_any_base')
        self.temp_dir = A('temp_dir') or filesnpaths.get_temp_directory_path()
        self.variable_buffer_length = A('variable_buffer_length')
        self.departure_from_reference_percentage = A('departure_from_reference_percentage')
        self.gene_caller_to_consider_in_context = A('gene_caller') or 'prodigal'
        self.snv_window_size = A('snv_window_size') or 50
        self.snv_window_step = A('snv_window_step') or 10
        self.minimum_snv_density = A('minimum_snv_density') or 0.1
        self.repeat_threshold = A('repeat_threshold') or 0.5
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
        self.skip_compute_DGR_variability_profiling = A('skip_compute_DGR_variability_profiling')
        self.pre_computed_dgrs_path = A('pre_computed_dgrs')
        self.initial_primer_length = A('initial_variable_primer_length') or  12 #TODO test different values for this. If Illumina reads are 250 bases then depends on length of VR
        self.numb_imperfect_tandem_repeats = A('numb_imperfect_tandem_repeats') or 10
        self.repeat_motif_coverage = A('repeat_motif_coverage') or 0.8
        self.snv_matching_proportion = A('snv_matching_proportion') or None
        self.snv_codon_position = A('snv_codon_position') or 0.33 # default is 33% of SNVs in the third codon position
        self.max_alignment_gaps = A('max_alignment_gaps') or 0

        # Detection mode parameters
        self.detection_mode = A('detection_mode')  # Will be resolved in sanity_check based on available inputs
        self.rt_window_size = A('rt_window_size') or 2000  # bp on each side of RT gene

        # Per-sample SNV analysis thresholds (hardcoded for now, no CLI args)
        self.min_snvs_per_sample = 5
        self.max_pct_snv_codon_3_per_sample = 33  # percent
        self.min_pct_snvs_explained_per_sample = 80  # percent

        # High confidence thresholds for SNV analysis
        self.high_conf_max_pct_snv_codon_3 = 10  # percent
        self.high_conf_min_pct_snvs_explained = 90  # percent (i.e., < 10% unexplained)

        # performance
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        if self.num_threads:
            self.run.info('Threads Used', self.num_threads)
        # be talkative or not
        self.verbose = A('verbose')
        self.run.info('Detection mode', self.detection_mode if self.detection_mode else "(will be determined based on inputs)")
        self.run.info('RT window size (bp each side)', self.rt_window_size)
        self.run.info('BLASTn word size', self.word_size)
        self.run.info('Skip "N" characters', self.skip_Ns)
        self.run.info('Skip "-" characters', self.skip_dashes)
        self.run.info('Discovery mode', self.discovery_mode)
        self.run.info('minimum_snv_density', self.minimum_snv_density)
        self.run.info('Number of Mismatches', self.number_of_mismatches)
        self.run.info('Initial mismatch bias threshold', self.initial_mismatch_bias_threshold)
        self.run.info('Max Non-Dominant Base Threshold', self.max_non_dominant)
        self.run.info('Minimum VR length after trimming', self.minimum_vr_length)
        self.run.info('Max alignment gaps', self.max_alignment_gaps)
        if self.allow_any_base:
            self.run.info('Allow any dominant base', self.allow_any_base)
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
        self.run.info('SNV window size', self.snv_window_size)
        self.run.info('SNV window step', self.snv_window_step)
        self.run.info('Variable buffer length', self.variable_buffer_length)
        self.run.info('Departure from reference percentage', self.departure_from_reference_percentage)
        if self.snv_matching_proportion:
            self.run.info('SNV matching proportion', self.snv_matching_proportion)
        self.run.info('Per-sample SNV: min SNVs', self.min_snvs_per_sample)
        self.run.info('Per-sample SNV: max % codon 3', self.max_pct_snv_codon_3_per_sample)
        self.run.info('Per-sample SNV: min % explained', self.min_pct_snvs_explained_per_sample)
        self.run.info('HMM(s) Provided', ", ".join(self.hmm) if self.hmm else "(will use default: Reverse_Transcriptase)")
        if not self.skip_recovering_genomic_context:
            self.run.info('Number of genes to consider in context', self.num_genes_to_consider_in_context)
        # computing variability profiling for every VR in every DGR by searching through raw reads?
        if not self.skip_compute_DGR_variability_profiling:
            self.run.info('Samples.txt', self.samples_txt)
            self.run.info('Initial Primer Length', self.initial_primer_length)
            self.run.info('Variable Region Primer Length', self.whole_primer_length)

        # these are the keys we are interested in finding in input files offered to reconstruct
        # DGR profiles via the --pre-computed-dgrs flag. NOTE that these keys are not ALL
        # keys that are used to build dgr profiles in the code, this way, the user can
        # attempt to characterize the activity of dgrs found in a single sample if they wish:
        self.essential_keys_to_describe_dgrs = [('DGR', str), ('VR', str), ('VR_contig', str), ('VR_frame_reported', int),
                                                ('VR_sequence', str), ('Midline', str), ('VR_start_position', int),
                                                ('VR_end_position', int), ('VR_bin', str), ('Mismatch %', float),
                                                ('TR_contig', str), ('TR_frame_Reported', int), ('TR_sequence', str),
                                                ('Base', str), ('Reverse Complemented_from_BLAST', bool), ('TR_start_position', int),
                                                ('TR_end_position', int), ('TR_bin', str), ('TR_in_gene', bool), ('HMM_source', str),
                                                ('distance_to_HMM', float), ('HMM_gene_name', str), ('HMM_direction', str), ('HMM_start', int),
                                                ('HMM_stop', int), ('HMM_gene_callers_id', int),
                                                ('numb_of_mismatches', int), ('numb_of_SNVs', int), ('VR_TR_mismatch_positions', list),
                                                ('snv_VR_positions', list), ('best_amongst_multiple_TRs_for_one_VR', bool),
                                                ('mismatch_codon_1', int), ('mismatch_codon_2', int), ('mismatch_codon_3', int),
                                                ('pct_mismatch_codon_3', float), ('vr_gene_id', str), ('n_snvs_total', int),
                                                ('n_snvs_explained', int), ('n_snvs_unexplained', int), ('n_explained_diverse', int), ('pct_snvs_explained', float),
                                                ('snv_codon_1', int), ('snv_codon_2', int), ('snv_codon_3', int),
                                                ('pct_snv_codon_3', float), ('confidence', str), ('confidence_reasons', list),
                                                ('detection_method', str)]



    def sanity_check(self):
        """
        Basic checks for a smooth operation
        """

        # Contigs.db is always required
        utils.is_contigs_db(self.contigs_db_path)

        # Check if profile.db exists (may be optional depending on detection mode)
        profile_db_exists = self.profile_db_path and os.path.exists(self.profile_db_path)
        if profile_db_exists:
            utils.is_profile_db(self.profile_db_path)

        # ========== DETECTION MODE RESOLUTION ==========
        # Resolve detection_mode based on user input and available data
        if self.detection_mode is None:
            # User didn't specify, determine based on available inputs
            if profile_db_exists:
                self.detection_mode = 'both'
                self.run.info_single("No --detection-mode specified. Since a profile.db is provided, "
                                    "anvi'o will use 'both' activity and homology-based detection.", nl_before=1)
            else:
                self.detection_mode = 'homology'
                self.run.warning("No --detection-mode specified and no profile.db provided. Anvi'o will use "
                               "'homology' mode only, which searches for DGRs near Reverse Transcriptase genes. "
                               "To use activity-based detection (which finds active DGRs using SNV patterns), "
                               "please provide a profile.db.",
                               header="DEFAULTING TO HOMOLOGY-ONLY MODE")
        elif self.detection_mode == 'activity':
            # Activity mode requires profile.db
            if not profile_db_exists:
                raise ConfigError("You requested --detection-mode 'activity', but no profile.db was provided. "
                                "Activity-based detection requires SNV data from a profile database. Please provide "
                                "a profile.db with --profile-db, or use --detection-mode 'homology' instead.")
        elif self.detection_mode == 'both':
            # Both mode prefers profile.db, but can fall back to homology-only
            if not profile_db_exists:
                self.detection_mode = 'homology'
                self.run.warning("You requested --detection-mode 'both', but no profile.db was provided. "
                               "Anvi'o is switching to 'homology' mode only. To use both detection modes, "
                               "please provide a profile.db with --profile-db.",
                               header="SWITCHING TO HOMOLOGY-ONLY MODE")
        # else: detection_mode == 'homology', no profile.db needed

        # Update the displayed detection mode now that it's resolved
        self.run.info('Detection mode (resolved)', self.detection_mode)

        # Validate rt_window_size
        if self.rt_window_size <= 0:
            raise ConfigError("The RT window size (--rt-window-size) must be a positive integer.")

        # ========== OUTPUT DIRECTORY ==========
        try:
            output_dir = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=False or self.just_do_it)

            if output_dir:
                filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=False)
        except Exception:
            raise ConfigError(f"Hold up your directory ({self.output_directory}) already exists. To avoid overwriting data we are stopping "
                            "you here. You have three options if you want to continue: rename your directory, delete your existing directory, "
                            "or rerun with the flag `--just-do-it` (which will overwrite your directory. You have been warned).")

        # ========== HMM CONFIGURATION ==========
        # Set default HMM to Reverse_Transcriptase (required for all modes)
        if not self.hmm:
            self.hmm = ['Reverse_Transcriptase']
            self.run.warning("No HMM source was specified with `--hmm-usage`. Anvi'o will use the default "
                           "'Reverse_Transcriptase' HMM to identify DGRs. If you haven't run this HMM yet, "
                           "please run `anvi-run-hmms -I Reverse_Transcriptase` on your contigs database first.",
                           header="USING DEFAULT HMM")

        if int(self.word_size) < 0:
            raise ConfigError('The word size value you are trying to input should be positive integer.')

        if self.variable_buffer_length < 0:
            raise ConfigError('The variable buffer length value you are trying to input should be positive integer.')

        if self.snv_window_size <= 0:
            raise ConfigError('The SNV window size must be a positive integer.')

        if self.snv_window_step <= 0:
            raise ConfigError('The SNV window step must be a positive integer.')

        if not (0 < self.initial_mismatch_bias_threshold <= 1):
            raise ConfigError('The initial mismatch bias threshold must be between 0 and 1 (exclusive of 0).')

        if not (self.max_non_dominant > 0):
            raise ConfigError('The max non-dominant base threshold needs to be a positive number.')

        if self.minimum_vr_length <= 0:
            raise ConfigError('The minimum VR length must be a positive integer.')

        if self.departure_from_reference_percentage < 0:
            raise ConfigError('The departure from reference percentage value you are trying to input should be a positive decimal number.')

        if self.min_mismatching_base_types_vr >= 5:
            raise ConfigError('The number of mismatching base types of the VR sequence cannot exceed 4 this is because there are only 4 bases in our DNA alphabet')

        if self.min_base_types_tr >= 5:
            raise ConfigError('The number of base types of the TR sequence cannot exceed 4 this is because there are only 4 bases in our DNA alphabet')

        if self.min_base_types_vr >= 5:
            raise ConfigError('The number of base types of the VR sequence cannot exceed 4 this is because there are only 4 bases in our DNA alphabet')

        if self.collections_mode:
            # Collections mode requires activity-based detection (needs profile.db)
            if self.detection_mode == 'homology':
                raise ConfigError("Collections mode (--collections-mode) requires activity-based detection, "
                                "but you are running in 'homology' mode only. Please provide a profile.db "
                                "and use --detection-mode 'activity' or 'both' to use collections mode.")
            if not self.collections_given:
                raise ConfigError("You must provide a collection name for collections mode to work. If you want to know about "
                                "all the collections in your profile database you can use the program "
                                "`anvi-show-collections-and-bins` or run the same command with the flag "
                                "`--list-collections`.")
            # ensure collections_given is a single string (collection name)
            if not isinstance(self.collections_given, str):
                raise ValueError("'collection-name' must be a single collection name")

        # ========== HMM VALIDATION ==========
        # HMM is required for all detection modes (for RT location in activity mode post-hoc,
        # and as the search anchor in homology mode)
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        hmm_hits_info_dict = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
        hmm_hits_info_dict = contigs_db.db.smart_get(t.hmm_hits_info_table_name, column = 'source')
        self.hmms_provided = set(hmm_hits_info_dict.keys())
        contigs_db.disconnect()

        bad_hmm = []
        for i in self.hmm:
            if i not in self.hmms_provided:
                bad_hmm.append(i)

        if bad_hmm == ['Reverse_Transcriptase']:
            raise ConfigError("The classic reverse transcriptase HMM has not been run with these contigs. "
                                "Please run 'anvi-run-hmms' with the '-I' flag and the 'Reverse_Transcriptase' HMM, "
                                f"so that anvi'o knows where there are Reverse Transcriptases in your {self.contigs_db_path}.")
        elif len(bad_hmm) > 0:
            if 'Reverse_Transcriptase' in bad_hmm:
                raise ConfigError("The classic reverse transcriptase HMM has not been run with these contigs. "
                                "Please run 'anvi-run-hmms' with the '-I' flag and the 'Reverse_Transcriptase' HMM, "
                                f"so that anvi'o knows where there are Reverse Transcriptases in your {self.contigs_db_path}. "
                                f"Also you don't have these HMMs you requested: {bad_hmm} in {self.hmms_provided}")
            else:
                raise ConfigError(f"You requested these HMMs to be searched through: {bad_hmm} in these HMMs {self.hmms_provided} "
                                f"that are in your {self.contigs_db_path}. The HMMs you give 'anvi-report-dgrs' need to be in your "
                                "contigs.db.")

        # Load gene information once (validates gene caller source and creates lookup structures)
        # This is needed for codon position analysis in all modes
        self.load_gene_info()

        html_files_exist = any(file.endswith('.html') for file in os.listdir(self.output_directory) if os.path.isfile(os.path.join(self.output_directory, file)))
        if html_files_exist:
            raise ConfigError("Files with .html suffix exist in the directory. Please delete them before rerunning and we will keep calm and carry on. (later this will delete them for you)")

        if self.initial_primer_length <  0:
            raise ConfigError("The initial primer length is set to a negative value or zero. This is not allowed. Please set the initial primer length to a positive value.")

        if self.whole_primer_length <= 0:
                raise ConfigError("The whole primer length is set to a negative value or zero. This is not allowed. Please set the whole primer length to a positive value.")

        if not self.samples_txt:
            # No samples.txt provided - skip variability profiling automatically
            self.skip_compute_DGR_variability_profiling = True
            self.run.warning("No samples.txt was provided, and THAT IS FINE. Anvi'o will skip computing DGR variability "
                           "profiling. If you want to compute variability profiling later, you can re-run the program "
                           "with a samples.txt file and use the `--pre-computed-dgrs` flag pointing to the DGRs_found.tsv "
                           "output from this run, so anvi'o only computes the variability profiling step.",
                           header="SKIPPING VARIABILITY PROFILING")
        elif not self.skip_compute_DGR_variability_profiling:
            # samples.txt provided and user wants variability profiling
            self.samples_artifact = SamplesTxt(self.samples_txt)

        if self.pre_computed_dgrs_path:
            filesnpaths.is_file_tab_delimited(self.pre_computed_dgrs_path)

            columns = utils.get_columns_of_TAB_delim_file(self.pre_computed_dgrs_path, include_first_column=True)
            # Find missing essential keys
            missing_keys = [k[0] for k in self.essential_keys_to_describe_dgrs if k[0] not in columns]

            if missing_keys:
                raise ConfigError(
                    f"The pre-computed DGR output file you provided does not look like one generated by "
                    f"`anvi-report-dgrs` :/\n\n"
                    f"missing required columns: {', '.join(missing_keys)}\n"
                    f"columns found in file ({len(columns)} total): {', '.join(columns)}\n\n"
                )



    def get_blast_results(self):
        """
        This function runs either the standard or the collections mode pipeline for the overall BLASTn xml file of results.
        It first identifies candidate variable regions with high SNV density and runs the BLASTn search, which
        uses sequences with high SNVs as the query and the contigs.db sequences as the target sequences.

        Returns
        =======
        blast_output : str
            Path to the BLASTn output file in XML format. The filename depends
            on whether collections mode is activated.
        Raises
        =======
        ConfigError
            If no SNV windows are found that meet the filtering criteria.

        """

        # initialise the SNV table
        self.init_snv_table()

        if self.collections_mode:
            self.run.info_single("Collections mode activated. Get ready to see as many BLASTn as bins in your collection. Big things be happenin'.", nl_before=1)
            return self.process_collections_mode()
        else:
            return self.process_standard_mode()



    def load_data_and_setup(self, bin_splits_list=None):
        """
        Load contig sequences and SNV data from databases, apply filtering,
        and optionally restrict to a specific bin.

        Parameters
        ==========
        bin_splits_list : list of str, optional
            Split names to restrict SNV and contig data to (used in collections mode).

        Returns
        =======
        sample_id_list : list of str
            Unique sample IDs with SNVs after filtering.
        contig_sequences : dict
            Dictionary mapping contig names to their sequences.
        """
        # load contig sequences
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
        contigs_db.disconnect()

        # filter for collections if needed
        if bin_splits_list:
            self.split_names_unique = bin_splits_list
            # get bin-specific contig sequences
            bin_contigs = [split.split('_split_')[0] for split in bin_splits_list]
            contig_sequences = {contig: contig_sequences[contig]
                                    for contig in bin_contigs if contig in contig_sequences}
        else:
            self.split_names_unique = utils.get_all_item_names_from_the_database(self.profile_db_path)

        sample_id_list = self.snv_panda.sample_id.unique().tolist()

        return sample_id_list, contig_sequences



    def init_snv_table(self):
        """
        Initialise the snv table but only the columns that we need.

        Returns
        =======
        self.snv_panda : pandas df
            Dataframe containing all of the snv information

        """
        # load SNV data
        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        columns_of_interest = [
            'sample_id', 'split_name', 'pos_in_contig', 'base_pos_in_codon',
            'departure_from_reference', 'reference'
        ] + nucleotides
        self.snv_panda = profile_db.db.get_table_as_dataframe(
            t.variable_nts_table_name,
            columns_of_interest=columns_of_interest
        ).sort_values(by=['split_name', 'pos_in_contig'])
        self.snv_panda['contig_name'] = self.snv_panda['split_name'].str.split('_split_').str[0]
        profile_db.disconnect()

        # apply SNV filters
        if self.discovery_mode:
            self.run.info_single("Running discovery mode. Search for SNVs in all possible locations. You go Dora the explorer!")
            self.snv_panda = self.snv_panda.query("departure_from_reference >= @self.departure_from_reference_percentage")
        else:
            self.snv_panda = self.snv_panda.query("departure_from_reference >= @self.departure_from_reference_percentage and base_pos_in_codon in (1, 2)")


    def load_gene_info(self):
        """
        Load gene information from contigs database once and create optimized lookup structures.

        Creates two data structures:
        - self.genes_in_contigs: dict keyed by gene_callers_id for direct lookups
        - self.gene_positions: dict of sorted lists per contig for efficient overlap queries
                               {contig_name: [(start, stop, direction, gene_id), ...]}

        Also validates that the requested gene caller source exists in the database.

        Raises
        ======
        ConfigError
            If the requested gene caller source is not found in the database.
        """
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        all_genes = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        contigs_db.disconnect()

        # Validate gene caller source exists
        unique_sources = set(gene_info['source'] for gene_info in all_genes.values())
        if self.gene_caller_to_consider_in_context not in unique_sources:
            raise ConfigError(f"Anvi'o can't find {self.gene_caller_to_consider_in_context} in your {self.contigs_db_path}. "
                            f"Here are the sources of your genes: {list(unique_sources)}.")

        # Filter by gene caller and build both data structures
        self.genes_in_contigs = {}
        self.gene_positions = defaultdict(list)

        for gene_id, gene_info in all_genes.items():
            if gene_info['source'] != self.gene_caller_to_consider_in_context:
                continue

            # Store full gene info for direct lookups by gene_id
            self.genes_in_contigs[gene_id] = gene_info

            # Build position-sorted structure for overlap queries
            contig = gene_info['contig']
            start = gene_info['start']
            stop = gene_info['stop']
            direction = gene_info['direction']
            self.gene_positions[contig].append((start, stop, direction, gene_id))

        # Sort by start position for binary search
        for contig in self.gene_positions:
            self.gene_positions[contig].sort(key=lambda x: x[0])


    def get_snvs_for_region(self, contig_name, start_pos, end_pos):
        """
        Fetch all SNVs (including 3rd codon position) for a specific genomic region.

        This is called only for candidate DGRs that pass initial filtering,
        so it's called ~10-50 times, not millions.

        Parameters
        ==========
        contig_name : str
            Name of the contig
        start_pos : int
            Start position in contig (inclusive)
        end_pos : int
            End position in contig (inclusive)

        Returns
        =======
        dict or None
            Dictionary with numpy arrays:
                'positions': np.array of pos_in_contig
                'codon_pos': np.array of base_pos_in_codon
                'reference': np.array of reference base
                'sample_ids': np.array of sample_id
                'A', 'C', 'G', 'T': np.array of nucleotide coverage counts
            Returns None if no SNVs found.
        """
        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        # Query for specific region using split_name LIKE pattern (no codon position filter)
        where_clause = f'''split_name LIKE "{contig_name}_split_%" AND pos_in_contig >= {start_pos} AND pos_in_contig <= {end_pos} AND departure_from_reference >= {self.departure_from_reference_percentage}'''

        snv_data = profile_db.db.get_some_rows_from_table_as_dict(
            t.variable_nts_table_name,
            where_clause=where_clause,
            error_if_no_data=False
        )
        profile_db.disconnect()

        if not snv_data:
            return None

        # Convert to sorted numpy arrays
        rows = sorted(snv_data.values(), key=lambda x: x['pos_in_contig'])

        return {
            'positions': np.array([r['pos_in_contig'] for r in rows]),
            'codon_pos': np.array([r['base_pos_in_codon'] for r in rows]),
            'reference': np.array([r['reference'] for r in rows]),
            'sample_ids': np.array([r['sample_id'] for r in rows]),
            'A': np.array([r['A'] for r in rows]),
            'C': np.array([r['C'] for r in rows]),
            'G': np.array([r['G'] for r in rows]),
            'T': np.array([r['T'] for r in rows]),
        }


    def get_codon_position(self, contig_pos, gene_start, gene_stop, direction):
        """
        Calculate the codon position (1, 2, or 3) for a given contig position.

        Parameters
        ==========
        contig_pos : int
            Position in the contig (0-indexed)
        gene_start : int
            Gene start position (0-indexed)
        gene_stop : int
            Gene stop position (0-indexed)
        direction : str
            'f' for forward, 'r' for reverse

        Returns
        =======
        int
            1, 2, or 3 indicating codon position, or 0 if position is not in gene
        """
        if contig_pos < gene_start or contig_pos >= gene_stop:
            return 0

        if direction == 'f':
            offset = contig_pos - gene_start
        else:
            # Reverse strand: count from the other end
            offset = gene_stop - 1 - contig_pos

        return (offset % 3) + 1  # Returns 1, 2, or 3


    def find_overlapping_gene(self, contig, position):
        """
        Find gene that overlaps with a given position using binary search.

        Parameters
        ==========
        contig : str
            Contig name
        position : int
            Position in contig (0-indexed)

        Returns
        =======
        tuple or None
            (gene_start, gene_stop, direction, gene_id) if found, None otherwise
        """
        if contig not in self.gene_positions:
            return None

        genes = self.gene_positions[contig]

        # Binary search: find genes where start <= position
        idx = bisect.bisect_right(genes, (position, float('inf'), '', 0)) - 1

        # Check if this gene contains the position
        while idx >= 0:
            start, stop, direction, gene_id = genes[idx]
            if start <= position < stop:
                return (start, stop, direction, gene_id)
            if stop <= position:
                break
            idx -= 1

        # Also check next gene in case of overlapping genes
        idx = bisect.bisect_right(genes, (position, float('inf'), '', 0))
        if idx < len(genes):
            start, stop, direction, gene_id = genes[idx]
            if start <= position < stop:
                return (start, stop, direction, gene_id)

        return None


    def analyze_mismatch_codon_positions(self, mismatch_positions_contig, query_contig):
        """
        Analyze codon positions of VR/TR mismatches (no SNV data needed).

        Parameters
        ==========
        mismatch_positions_contig : list
            Mismatch positions in contig coordinates
        query_contig : str
            Contig name

        Returns
        =======
        dict
            'mismatch_codon_1': count of mismatches at 1st codon position
            'mismatch_codon_2': count of mismatches at 2nd codon position
            'mismatch_codon_3': count of mismatches at 3rd codon position
            'mismatch_codon_unknown': count of mismatches not in a gene
            'pct_mismatch_codon_3': percentage at 3rd position (of those in genes)
            'vr_gene_id': overlapping gene if found
        """
        counts = {'codon_1': 0, 'codon_2': 0, 'codon_3': 0, 'codon_unknown': 0}
        gene_id = None

        for pos in mismatch_positions_contig:
            gene_info = self.find_overlapping_gene(query_contig, pos)

            if gene_info is None:
                counts['codon_unknown'] += 1
            else:
                gene_start, gene_stop, direction, gid = gene_info
                gene_id = gid  # Store last found gene
                codon_pos = self.get_codon_position(pos, gene_start, gene_stop, direction)
                counts[f'codon_{codon_pos}'] += 1

        # Calculate percentage (only among positions in genes)
        in_gene = counts['codon_1'] + counts['codon_2'] + counts['codon_3']
        pct_codon_3 = (counts['codon_3'] / in_gene * 100) if in_gene > 0 else 0

        return {
            'mismatch_codon_1': counts['codon_1'],
            'mismatch_codon_2': counts['codon_2'],
            'mismatch_codon_3': counts['codon_3'],
            'mismatch_codon_unknown': counts['codon_unknown'],
            'pct_mismatch_codon_3': pct_codon_3,
            'vr_gene_id': gene_id
        }


    def analyze_snv_pattern(self, query_contig, query_start, query_end,
                            mismatch_positions_contig, mutagenesis_base, is_reverse_complement):
        """
        Analyze SNV pattern on a per-sample basis.

        Evaluates each sample individually, sorted by SNV count (descending).
        Returns metrics from the first sample that passes all thresholds,
        or from the sample with most SNVs if none pass.

        Parameters
        ==========
        query_contig : str
            Contig name
        query_start : int
            VR region start position
        query_end : int
            VR region end position
        mismatch_positions_contig : list
            Mismatch positions in contig coordinates
        mutagenesis_base : str
            Expected mutagenesis base ('A' typically)
        is_reverse_complement : bool
            Whether alignment was reverse complemented

        Returns
        =======
        dict or None
            SNV analysis metrics including:
            - All standard fields (n_snvs_total, pct_snvs_explained, etc.)
            - 'snv_supporting_sample': str or None (sample that passed thresholds)
            - 'sample_passed': bool (True if a sample met all thresholds)
        """
        # Fetch full SNV data for this region (including 3rd codon position)
        snv_data = self.get_snvs_for_region(query_contig, query_start, query_end)

        if snv_data is None or len(snv_data['positions']) == 0:
            return None

        positions = snv_data['positions']
        codon_pos = snv_data['codon_pos']
        reference = snv_data['reference']
        sample_ids = snv_data['sample_ids']
        # Nucleotide coverage arrays for diversity calculation
        nuc_A = snv_data['A']
        nuc_C = snv_data['C']
        nuc_G = snv_data['G']
        nuc_T = snv_data['T']

        # Adjust mutagenesis base for reverse complement
        expected_base = mutagenesis_base.upper()
        if is_reverse_complement:
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
            expected_base = complement.get(expected_base, expected_base)

        mismatch_set = set(mismatch_positions_contig)

        # Group SNV indices by sample
        sample_indices = defaultdict(list)
        for i, sample_id in enumerate(sample_ids):
            sample_indices[sample_id].append(i)

        # Sort samples by SNV count (descending) for early exit optimization
        # Samples with more SNVs are checked first; we stop as soon as one passes
        samples_sorted = sorted(sample_indices.keys(),
                               key=lambda s: len(sample_indices[s]),
                               reverse=True)

        best_result = None
        best_sample_snv_count = 0

        for sample_id in samples_sorted:
            indices = sample_indices[sample_id]

            # Get unique positions for this sample (same position can appear multiple times)
            sample_positions = []
            seen = set()
            for i in indices:
                pos = positions[i]
                if pos not in seen:
                    seen.add(pos)
                    sample_positions.append(i)

            n_snvs = len(sample_positions)

            # Skip samples below minimum threshold - not enough data to evaluate
            if n_snvs < self.min_snvs_per_sample:
                continue

            # Compute metrics for this sample
            snv_codon_counts = {1: 0, 2: 0, 3: 0, 0: 0}  # 0 = unknown/not in gene
            n_explained = 0
            n_unexplained = 0
            n_explained_diverse = 0  # explained SNVs with >= 3 nucleotides
            unique_positions_list = []

            for i in sample_positions:
                pos = positions[i]
                unique_positions_list.append(pos)

                # Codon position counting
                cp = codon_pos[i] if i < len(codon_pos) else 0
                if cp in snv_codon_counts:
                    snv_codon_counts[cp] += 1

                # Check nucleotide diversity at this position (>= 3 bases with coverage > 0)
                n_bases_with_coverage = sum(1 for cov in [nuc_A[i], nuc_C[i], nuc_G[i], nuc_T[i]] if cov > 0)
                is_diverse = n_bases_with_coverage >= 3

                # Explained: at mismatch position OR at mutagenesis base site
                # Unexplained: at a position that doesn't match either criterion
                if pos in mismatch_set:
                    n_explained += 1
                    if is_diverse:
                        n_explained_diverse += 1
                elif reference[i] == expected_base:
                    n_explained += 1
                    if is_diverse:
                        n_explained_diverse += 1
                else:
                    n_unexplained += 1

            # Calculate percentages
            pct_explained = (n_explained / n_snvs * 100) if n_snvs > 0 else 100
            in_gene = snv_codon_counts[1] + snv_codon_counts[2] + snv_codon_counts[3]
            pct_snv_codon_3 = (snv_codon_counts[3] / in_gene * 100) if in_gene > 0 else 0

            # Build result dict for this sample
            sample_result = {
                'n_snvs_total': n_snvs,
                'n_snvs_explained': n_explained,
                'n_snvs_unexplained': n_unexplained,
                'n_explained_diverse': n_explained_diverse,
                'pct_snvs_explained': pct_explained,
                'snv_codon_1': snv_codon_counts[1],
                'snv_codon_2': snv_codon_counts[2],
                'snv_codon_3': snv_codon_counts[3],
                'pct_snv_codon_3': pct_snv_codon_3,
                'snv_positions': unique_positions_list,
                'snv_supporting_sample': sample_id,
                'sample_passed': False  # Will be set to True if passes
            }

            # Track best result (sample with most SNVs) as fallback
            if n_snvs > best_sample_snv_count:
                best_sample_snv_count = n_snvs
                best_result = sample_result

            # Check if this sample passes all thresholds
            passes_codon_3 = pct_snv_codon_3 <= self.max_pct_snv_codon_3_per_sample
            passes_explained = pct_explained >= self.min_pct_snvs_explained_per_sample

            if passes_codon_3 and passes_explained:
                # Found a passing sample - return immediately (early exit)
                sample_result['sample_passed'] = True
                return sample_result

        # No sample passed thresholds - return best result (most SNVs) with sample_passed=False
        if best_result is not None:
            best_result['snv_supporting_sample'] = None  # No supporting sample found
            best_result['sample_passed'] = False
            return best_result

        # No samples had enough SNVs to evaluate
        return None


    def compute_dgr_confidence(self, mismatch_analysis, snv_analysis):
        """
        Compute confidence level based on mismatch and SNV analyses.

        Confidence levels:
        - high: Good mismatch distribution AND a sample passed with excellent SNV metrics
        - medium: Good mismatch distribution AND a sample passed with acceptable SNV metrics
        - low: Poor mismatch distribution OR no sample passed SNV thresholds OR no SNVs

        Parameters
        ==========
        mismatch_analysis : dict
            Results from analyze_mismatch_codon_positions()
        snv_analysis : dict or None
            Results from analyze_snv_pattern(), or None if no SNVs

        Returns
        =======
        tuple
            (confidence: str, reasons: list)
            confidence is 'high', 'medium', or 'low'
            reasons is a list of strings explaining any downgrade
        """
        confidence = 'high'
        reasons = []

        # Check mismatch codon distribution (unchanged from before)
        # If >50% of VR/TR mismatches are at 3rd codon position, it's suspicious
        if mismatch_analysis['pct_mismatch_codon_3'] > 50:
            confidence = 'medium'
            reasons.append(f"high_mismatch_codon_3_{mismatch_analysis['pct_mismatch_codon_3']:.1f}%")

        # Check SNV pattern using per-sample logic
        if snv_analysis is None:
            # No SNVs at all in VR region -> low confidence
            # This could mean no coverage or no variation detected
            confidence = 'low'
            reasons.append("no_snvs_in_vr_region")
        elif not snv_analysis.get('sample_passed', False):
            # SNVs exist but no sample passed the thresholds -> low confidence
            # This suggests the SNV pattern doesn't match DGR expectations
            confidence = 'low'
            reasons.append("no_sample_passed_snv_thresholds")
        else:
            # A sample passed - check if it has excellent metrics for high confidence
            pct_codon_3 = snv_analysis['pct_snv_codon_3']
            pct_explained = snv_analysis['pct_snvs_explained']
            n_explained_diverse = snv_analysis.get('n_explained_diverse', 0)

            # Check diversity of explained SNVs (>= 3 nucleotides at position)
            # Low diversity (<=2 diverse SNVs) suggests false positive
            if n_explained_diverse <= 2:
                confidence = 'low'
                reasons.append(f"low_snv_diversity_n_explained_diverse_{n_explained_diverse}")
                return confidence, reasons

            # Excellent: very few SNVs at codon 3 AND almost all SNVs explained
            is_excellent = (pct_codon_3 <= self.high_conf_max_pct_snv_codon_3 and
                           pct_explained >= self.high_conf_min_pct_snvs_explained)

            if is_excellent:
                # Excellent SNV metrics - this strongly supports DGR activity
                # Note: we don't override mismatch-based downgrade, but we note the excellent SNVs
                if confidence == 'medium':
                    reasons.append(f"excellent_snv_support_codon3_{pct_codon_3:.1f}%_explained_{pct_explained:.1f}%")
                # If confidence is still 'high', keep it high (no reason to add)
            else:
                # Passed thresholds but not excellent -> cap at medium confidence
                if confidence == 'high':
                    confidence = 'medium'
                reasons.append(f"snv_passed_not_excellent_codon3_{pct_codon_3:.1f}%_explained_{pct_explained:.1f}%")

        return confidence, reasons


    def find_snv_clusters(self, sample_id_list, contig_sequences):
        """
        Detect clusters of SNVs within contigs, merges overlapping windows,and extracts subsequences as candidate variable regions.

        Parameters
        ==========
        sample_id_list : list of str
            Sample IDs containing SNVs to be processed.

        Returns
        =======
        contig_records : list of SeqRecord
            Extracted subsequences from contigs representing SNV clusters.

        Raises
        ======
        ConfigError
            If no subsequences with sufficient SNV density are found.
        """

        # reset possible windows each run
        self.all_possible_windows = {}

        # For SNV windows, only keep SNVs with >= 3 distinct nucleotides within a sample.
        min_diverse_bases = 3
        diverse_base_counts = (self.snv_panda[nucleotides] > 0).sum(axis=1)
        snv_panda_for_windows = self.snv_panda.loc[diverse_base_counts >= min_diverse_bases]

        # group the DataFrame by 'split_name' and 'sample_id' upfront
        grouped = snv_panda_for_windows.groupby(['split_name', 'sample_id'])

        # now iterate over the grouped data
        for (split, sample), group in grouped:
            # only process groups that match the desired split and sample
            if split in self.split_names_unique and sample in sample_id_list:
                if group.shape[0] == 0:
                    continue

                # extract the contig name and positions for the group
                contig_name = group.contig_name.unique()[0]
                pos_list = group.pos_in_contig.to_list()

                if contig_name not in self.all_possible_windows:
                    # if not, initialize it with an empty list
                    self.all_possible_windows[contig_name] = []

                # Sort positions for efficient binary search when counting SNVs in windows
                sorted_pos = sorted(pos_list)

                if not sorted_pos:
                    continue

                # Determine the range we need to scan with sliding windows
                min_pos = sorted_pos[0]
                max_pos = sorted_pos[-1]
                contig_len = len(contig_sequences[contig_name]['sequence'])

                # Start scanning from before the first SNV (to catch edge cases)
                # End at the last SNV position (windows starting after won't capture it)
                window_start_range = max(0, min_pos - self.snv_window_size)
                window_end_range = min(max_pos, contig_len - self.snv_window_size)

                # Slide the window across the region
                window_pos = window_start_range
                while window_pos <= window_end_range:
                    window_end = window_pos + self.snv_window_size

                    # Use binary search to count SNVs in window [window_pos, window_end)
                    # bisect_left returns the index where window_pos would be inserted
                    # This gives us O(log n) lookup instead of O(n)
                    left_idx = bisect.bisect_left(sorted_pos, window_pos)
                    right_idx = bisect.bisect_left(sorted_pos, window_end)
                    snv_count = right_idx - left_idx

                    # Calculate density as SNVs per window size (fixed denominator now)
                    snv_density = snv_count / self.snv_window_size

                    if snv_density >= self.minimum_snv_density:
                        # Add buffer around the window
                        buffered_start = max(0, window_pos - self.variable_buffer_length)
                        buffered_end = min(contig_len, window_end + self.variable_buffer_length)
                        self.all_possible_windows[contig_name].append((buffered_start, buffered_end))

                    # Move window by step size
                    window_pos += self.snv_window_step

        all_merged_snv_windows = {} # this dictionary will be filled up with the merged window list for each contig

        # merge overlapping windows using standard interval merging algorithm: O(n log n)
        # Key insight: once sorted by start, we only need to check overlap with the LAST merged window
        for contig_name, window_list in self.all_possible_windows.items():
            # sort by start position - O(n log n)
            sorted_windows = sorted(window_list, key=lambda x: x[0])

            # single-pass merge - O(n)
            merged_windows_in_contig = []
            for start, end in sorted_windows:
                if merged_windows_in_contig and start <= merged_windows_in_contig[-1][1]:
                    # current window overlaps with last merged window - extend it
                    prev_start, prev_end = merged_windows_in_contig[-1]
                    merged_windows_in_contig[-1] = (prev_start, max(prev_end, end))
                else:
                    # no overlap - start a new merged window
                    merged_windows_in_contig.append((start, end))

            all_merged_snv_windows[contig_name] = merged_windows_in_contig

        # get short sequences from all_merged_snv_window
        contig_records = {}
        for contig_name in all_merged_snv_windows.keys():
            contig_sequence=contig_sequences[contig_name]['sequence']
            self.positions= all_merged_snv_windows[contig_name]
            for i, (start, end) in enumerate(self.positions):
                section_sequence = contig_sequence[start:end]
                section_id = f"{contig_name}_section_{i}_start_bp{start}_end_bp{end}"
                # add sequence to dictionary of
                contig_records[section_id] = section_sequence
                # contig_records.append(SeqRecord(section_sequence, id=section_id, description=""))

        if len(contig_records) == 0:
            self.run.warning(f"No sequences with SNVs were found with the parameters window size:{self.snv_window_size}bp, "
                            f"step size:{self.snv_window_step}bp, and minimum SNV density:{self.minimum_snv_density}. "
                            f"This means there are no variable region candidates for a BLAST search"
                            , header="NO SEQUENCES WITH SUBSTANTIAL SNVS FOUND")
            raise ConfigError("Therefore, we will exit here because anvi'o has nothing to search for DGRs in, "
                            "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                            "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")

        return contig_records



    def run_blast(self, query_records, target_sequences, mode='activity', bin_name=None):
        """
        Run BLASTn with given query and target sequences.

        This is the unified BLAST execution method used by both activity-based and
        homology-based detection modes. It handles logging and delegates the actual
        BLAST execution to the standalone execute_blast() function.

        Parameters
        ==========
        query_records : dict
            Query sequences as {section_id: sequence}. For activity mode, these are
            SNV cluster regions (VR candidates). For homology mode, these are RT
            window regions (TR candidates).
        target_sequences : dict
            Target sequences as {contig_name: {'sequence': seq}}. Typically all
            contigs or bin-specific contigs.
        mode : str
            Detection mode: 'activity' or 'homology'. Used for logging messages
            and output filename generation.
        bin_name : str or None
            Bin name for collections mode. If provided, creates bin-specific output files.

        Returns
        =======
        blast_output_path : str
            Path to the BLASTn output XML file.

        Raises
        ======
        ConfigError
            If no query or target sequences are provided.
        """
        if not query_records:
            raise ConfigError("No query sequences provided for BLAST search.")

        if not target_sequences:
            raise ConfigError("No target sequences found. Please check your contigs database.")

        # Generate descriptive filenames based on mode and bin
        if bin_name:
            query_fasta_filename = f"blast_query_{mode}_bin_{bin_name}.fasta"
            target_fasta_filename = f"blast_target_{mode}_bin_{bin_name}.fasta"
            blast_output_filename = f"blast_output_{mode}_bin_{bin_name}_wordsize_{self.word_size}.xml"
        else:
            query_fasta_filename = f"blast_query_{mode}.fasta"
            target_fasta_filename = f"blast_target_{mode}.fasta"
            blast_output_filename = f"blast_output_{mode}_wordsize_{self.word_size}.xml"

        # Log query summary with mode-appropriate message
        num_sequences = len(query_records)
        total_bp = sum(len(seq) for seq in query_records.values())
        if mode == 'activity':
            query_desc = "VR candidate region(s) (SNV clusters)"
        else:
            query_desc = "TR candidate region(s) (RT windows)"

        bin_suffix = f" for bin '{bin_name}'" if bin_name else ""
        self.run.info_single(f"Identified {num_sequences} {query_desc} "
                            f"({total_bp:,} bp total){bin_suffix} for BLAST search.", nl_before=1)

        # Execute BLAST using the standalone function
        blast_output_path = execute_blast(
            query_records=query_records,
            target_sequences=target_sequences,
            temp_dir=self.temp_dir,
            word_size=self.word_size,
            num_threads=self.num_threads,
            query_fasta_filename=query_fasta_filename,
            target_fasta_filename=target_fasta_filename,
            blast_output_filename=blast_output_filename
        )

        self.run.info(f'BLAST output ({mode} mode)', blast_output_path)

        return blast_output_path



    @staticmethod
    def process_bin_worker(input_queue, output_queue, config, snv_panda_records, contigs_db_path, profile_db_path):
        """
        Worker function for parallel bin processing. Processes bins from input_queue
        and puts results in output_queue.

        This is a static method to work with multiprocessing - all needed data is passed
        explicitly rather than through self.
        """
        import pandas as pd
        import bisect

        # reconstruct snv_panda from records
        snv_panda = pd.DataFrame.from_records(snv_panda_records)

        while True:
            bin_name, bin_splits_list = input_queue.get(True)

            try:
                # === load_data_and_setup logic ===
                contigs_db = dbops.ContigsDatabase(contigs_db_path, run=run_quiet, progress=progress_quiet)
                contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
                contigs_db.disconnect()

                split_names_unique = bin_splits_list
                bin_contigs = [split.split('_split_')[0] for split in bin_splits_list]
                bin_contig_sequences = {contig: contig_sequences[contig]
                                        for contig in bin_contigs if contig in contig_sequences}

                sample_id_list = list(set(snv_panda.sample_id.unique()))

                # === find_snv_clusters logic ===
                min_diverse_bases = 3
                diverse_base_counts = (snv_panda[nucleotides] > 0).sum(axis=1)
                snv_panda_for_windows = snv_panda.loc[diverse_base_counts >= min_diverse_bases]
                all_possible_windows = {}
                grouped = snv_panda_for_windows.groupby(['split_name', 'sample_id'])

                for (split, sample), group in grouped:
                    if split in split_names_unique and sample in sample_id_list:
                        if group.shape[0] == 0:
                            continue

                        contig_name = group.contig_name.unique()[0]
                        pos_list = group.pos_in_contig.to_list()

                        if contig_name not in all_possible_windows:
                            all_possible_windows[contig_name] = []

                        # Sort positions for efficient binary search
                        sorted_pos = sorted(pos_list)

                        if not sorted_pos:
                            continue

                        # Determine the range to scan with sliding windows
                        min_pos = sorted_pos[0]
                        max_pos = sorted_pos[-1]
                        contig_len = len(bin_contig_sequences[contig_name]['sequence'])

                        window_start_range = max(0, min_pos - config['snv_window_size'])
                        window_end_range = min(max_pos, contig_len - config['snv_window_size'])

                        # Slide the window across the region
                        window_pos = window_start_range
                        while window_pos <= window_end_range:
                            window_end = window_pos + config['snv_window_size']

                            # Use binary search to count SNVs in window
                            left_idx = bisect.bisect_left(sorted_pos, window_pos)
                            right_idx = bisect.bisect_left(sorted_pos, window_end)
                            snv_count = right_idx - left_idx

                            # Calculate density
                            snv_density = snv_count / config['snv_window_size']

                            if snv_density >= config['minimum_snv_density']:
                                buffered_start = max(0, window_pos - config['variable_buffer_length'])
                                buffered_end = min(contig_len, window_end + config['variable_buffer_length'])
                                all_possible_windows[contig_name].append((buffered_start, buffered_end))

                            window_pos += config['snv_window_step']

                # merge overlapping windows
                all_merged_snv_windows = {}
                for contig_name, window_list in all_possible_windows.items():
                    sorted_windows = sorted(window_list, key=lambda x: x[0])
                    merged_windows = []

                    for window in sorted_windows:
                        if merged_windows and window[0] <= merged_windows[-1][1]:
                            merged_windows[-1] = (merged_windows[-1][0], max(merged_windows[-1][1], window[1]))
                        else:
                            merged_windows.append(window)

                    all_merged_snv_windows[contig_name] = merged_windows

                # extract subsequences
                contig_records = {}
                for contig_name in all_merged_snv_windows.keys():
                    contig_sequence = bin_contig_sequences[contig_name]['sequence']
                    for i, (start, end) in enumerate(all_merged_snv_windows[contig_name]):
                        section_sequence = contig_sequence[start:end]
                        section_name = f"{contig_name}_section{i}_start_bp{start}_end_bp{end}"
                        contig_records[section_name] = section_sequence

                if not contig_records:
                    output_queue.put((bin_name, False, None, "No SNV clusters found"))
                    continue

                # === run_blast logic ===
                # Use the standalone execute_blast() function to avoid code duplication
                blast_output_path = execute_blast(
                    query_records=contig_records,
                    target_sequences=bin_contig_sequences,
                    temp_dir=config['temp_dir'],
                    word_size=config['word_size'],
                    num_threads=1,
                    query_fasta_filename=f"bin_{bin_name}_subsequences.fasta",
                    target_fasta_filename=f"bin_{bin_name}_reference_sequences.fasta",
                    blast_output_filename=f"blast_output_for_bin_{bin_name}_wordsize_{config['word_size']}.xml"
                )

                output_queue.put((bin_name, True, blast_output_path, None))

            except Exception as e:
                output_queue.put((bin_name, False, None, str(e)))


    @staticmethod
    def process_bin_worker_homology(input_queue, output_queue, config, rt_windows_list, contigs_db_path):
        """
        Worker function for parallel bin processing in homology mode. Processes bins from input_queue
        and puts results in output_queue.

        This is a static method to work with multiprocessing - all needed data is passed
        explicitly rather than through self.

        Parameters
        ==========
        input_queue : multiprocessing.Queue
            Queue containing (bin_name, bin_contigs_list) tuples.
        output_queue : multiprocessing.Queue
            Queue for results as (bin_name, success, blast_path, error_msg) tuples.
        config : dict
            Configuration dict with 'temp_dir', 'word_size'.
        rt_windows_list : list
            List of RT window dicts (serializable format from rt_windows.values()).
        contigs_db_path : str
            Path to contigs database.
        """
        while True:
            bin_name, bin_contigs_list = input_queue.get(True)

            try:
                # Load contig sequences for this bin
                contigs_db = dbops.ContigsDatabase(contigs_db_path, run=run_quiet, progress=progress_quiet)
                contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
                contigs_db.disconnect()

                # Filter contig sequences to only those in this bin
                bin_contig_sequences = {contig: contig_sequences[contig]
                                        for contig in bin_contigs_list if contig in contig_sequences}

                if not bin_contig_sequences:
                    output_queue.put((bin_name, False, None, "No contig sequences found for bin"))
                    continue

                # Filter RT windows to only those in this bin's contigs
                query_records = {}
                for window_info in rt_windows_list:
                    contig = window_info['contig']
                    if contig not in bin_contigs_list:
                        continue

                    window_start = window_info['window_start']
                    window_end = window_info['window_end']
                    window_id = window_info['window_id']

                    if contig not in bin_contig_sequences:
                        continue

                    # Extract window sequence
                    contig_seq = bin_contig_sequences[contig]['sequence']
                    window_seq = contig_seq[window_start:window_end + 1]

                    # Create section_id compatible with existing parsing
                    section_id = f"{contig}_section_{window_id}_start_bp{window_start}_end_bp{window_end}"
                    query_records[section_id] = window_seq

                if not query_records:
                    output_queue.put((bin_name, False, None, "No RT windows found in bin"))
                    continue

                # Run BLAST using the standalone execute_blast() function
                blast_output_path = execute_blast(
                    query_records=query_records,
                    target_sequences=bin_contig_sequences,
                    temp_dir=config['temp_dir'],
                    word_size=config['word_size'],
                    num_threads=1,
                    query_fasta_filename=f"bin_{bin_name}_homology_query.fasta",
                    target_fasta_filename=f"bin_{bin_name}_homology_target.fasta",
                    blast_output_filename=f"blast_output_homology_for_bin_{bin_name}_wordsize_{config['word_size']}.xml"
                )

                output_queue.put((bin_name, True, blast_output_path, None))

            except Exception as e:
                output_queue.put((bin_name, False, None, str(e)))


    def process_collections_mode(self):
        """
        Process contigs in collections mode by running BLASTn separately
        for each bin in the collection.

        If num_bins > num_threads: process bins in parallel (each BLAST uses 1 thread)
        If num_bins <= num_threads: process bins sequentially (each BLAST uses all threads)

        Parameters
        ==========
        None

        Returns
        =======
        blast_output : str
            Path to the BLASTn output XML file from the last successfully processed bin.

        Raises
        ======
        ConfigError
            If no valid SNV clusters are found for a bin.
        """

        # get collections data
        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        where_clause = f'''collection_name == "{self.collections_given}"'''
        split_collections_dict = profile_db.db.get_some_rows_from_table_as_dict(
            t.collections_splits_table_name, where_clause=where_clause, error_if_no_data=True)
        profile_db.disconnect()

        # group splits by bin
        bin_splits_dict = {}
        for key, value in split_collections_dict.items():
            bin_name = value['bin_name']
            split = value['split']
            if bin_name not in bin_splits_dict:
                bin_splits_dict[bin_name] = []
            bin_splits_dict[bin_name].append(split)

        self.bin_names_list = list(bin_splits_dict.keys())
        num_bins = len(bin_splits_dict)

        # Initialize blast_output to None
        self.blast_output = None
        successful_bins = []
        skipped_bins = {}  # bin_name -> error message

        # === DECIDE: PARALLEL vs SEQUENTIAL ===
        use_parallel = num_bins > self.num_threads and self.num_threads > 1

        if use_parallel:
            # === PARALLEL BIN PROCESSING ===
            self.run.info_single(f"Processing {num_bins} bins in parallel using {self.num_threads} threads "
                                f"(each BLAST uses 1 thread).", nl_before=1)

            # prepare config dict for workers
            config = {
                'temp_dir': self.temp_dir,
                'word_size': self.word_size,
                'snv_window_size': self.snv_window_size,
                'snv_window_step': self.snv_window_step,
                'minimum_snv_density': self.minimum_snv_density,
                'variable_buffer_length': self.variable_buffer_length,
            }

            # convert snv_panda to records for pickling
            snv_panda_records = self.snv_panda.to_dict('records')

            # setup queues
            manager = multiprocessing.Manager()
            input_queue = manager.Queue()
            output_queue = manager.Queue()

            # put all bins in input queue
            for bin_name, bin_splits_list in bin_splits_dict.items():
                input_queue.put((bin_name, bin_splits_list))

            # start workers
            workers = []
            for i in range(self.num_threads):
                worker = multiprocessing.Process(
                    target=DGR_Finder.process_bin_worker,
                    args=(input_queue, output_queue, config, snv_panda_records,
                          self.contigs_db_path, self.profile_db_path))
                workers.append(worker)
                worker.start()

            # monitor progress
            self.progress.new('Processing bins for BLAST (parallel)', progress_total_items=num_bins)
            self.progress.update(f"Processing {num_bins} bins across {self.num_threads} workers...")

            bins_processed = 0
            while bins_processed < num_bins:
                try:
                    bin_name, success, blast_path, error_msg = output_queue.get()
                    bins_processed += 1
                    self.progress.increment(increment_to=bins_processed)

                    if success:
                        successful_bins.append(bin_name)
                        self.blast_output = blast_path
                    else:
                        skipped_bins[bin_name] = error_msg

                    if bins_processed < num_bins:
                        self.progress.update(f"Processed {bins_processed}/{num_bins} bins...")
                    else:
                        self.progress.update("All done!")

                except KeyboardInterrupt:
                    self.run.info_single("Received kill signal, terminating workers...", nl_before=1)
                    break

            self.progress.end()

            # terminate workers
            for worker in workers:
                worker.terminate()

        else:
            # === SEQUENTIAL BIN PROCESSING (original behavior) ===
            if self.num_threads > 1:
                self.run.info_single(f"Processing {num_bins} bins sequentially "
                                    f"(each BLAST uses {self.num_threads} threads).", nl_before=1)

            self.progress.new('Processing bins for BLAST', progress_total_items=num_bins)

            for bin_name, bin_splits_list in bin_splits_dict.items():
                self.progress.update(f"{bin_name}", increment=True)

                try:
                    sample_id_list, bin_contig_sequences = self.load_data_and_setup(bin_splits_list)

                    contig_records = self.find_snv_clusters(sample_id_list, bin_contig_sequences)
                    blast_output_path = self.run_blast(contig_records, bin_contig_sequences,
                                                       mode='activity', bin_name=bin_name)
                    self.blast_output = blast_output_path

                    successful_bins.append(bin_name)

                except ConfigError as e:
                    skipped_bins[bin_name] = str(e)
                    continue
                except Exception as e:
                    skipped_bins[bin_name] = str(e)
                    continue

            self.progress.end()

        # Check if any bins were successfully processed
        if not successful_bins:
            raise ConfigError(f"No bins in collection '{self.collections_given}' could be processed successfully. "
                            f"All {len(skipped_bins)} bins were skipped. "
                            f"Common causes: sequences too short for word_size={self.word_size}, "
                            f"insufficient SNV density, or BLAST failures.")

        # summary message
        self.run.info_single(f"Processed {PL('bin', num_bins)}: {len(successful_bins)} successful, "
                            f"{len(skipped_bins)} skipped.", nl_before=1)

        # report skipped bins if any
        if skipped_bins:
            self.run.warning(f"{PL('bin', len(skipped_bins))} skipped during BLAST: {', '.join(skipped_bins.keys())}. "
                           "Use '--debug' to see detailed error messages.",
                           header="BINS SKIPPED")
            if anvio.DEBUG:
                for bin_name, error_msg in skipped_bins.items():
                    self.run.info_single(f"  {bin_name}: {error_msg}")

        # return the last successful blast output (maintains original behavior)
        return self.blast_output


    def process_collections_mode_homology(self):
        """
        Process contigs in collections mode for homology-based detection by running BLASTn
        separately for each bin in the collection.

        Similar to process_collections_mode() but uses RT windows as queries instead of
        SNV clusters. This is used when detection_mode is 'homology' or 'both' and
        collections_mode is True.

        If num_bins > num_threads: process bins in parallel (each BLAST uses 1 thread)
        If num_bins <= num_threads: process bins sequentially (each BLAST uses all threads)

        Returns
        =======
        blast_outputs : dict
            Dictionary mapping bin_name to blast_output_path for successful bins.
        """
        self.run.info_single("Running homology-based detection in collections mode...", nl_before=1)

        # Get collections data
        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        where_clause = f'''collection_name == "{self.collections_given}"'''
        split_collections_dict = profile_db.db.get_some_rows_from_table_as_dict(
            t.collections_splits_table_name, where_clause=where_clause, error_if_no_data=True)
        profile_db.disconnect()

        # Group splits by bin and get contigs
        bin_contigs_dict = {}
        for key, value in split_collections_dict.items():
            bin_name = value['bin_name']
            split = value['split']
            contig = split.split('_split_')[0]
            if bin_name not in bin_contigs_dict:
                bin_contigs_dict[bin_name] = set()
            bin_contigs_dict[bin_name].add(contig)

        # Convert sets to lists
        for bin_name in bin_contigs_dict:
            bin_contigs_dict[bin_name] = list(bin_contigs_dict[bin_name])

        num_bins = len(bin_contigs_dict)

        # Pre-compute RT windows (shared across all bins)
        if not hasattr(self, 'rt_windows') or not self.rt_windows:
            self.get_rt_windows()

        if not self.rt_windows:
            self.run.warning("No RT windows found. Homology-based detection cannot proceed in collections mode.",
                           header="NO RT WINDOWS")
            return {}

        # Initialize results tracking
        blast_outputs = {}
        successful_bins = []
        skipped_bins = {}

        # === DECIDE: PARALLEL vs SEQUENTIAL ===
        use_parallel = num_bins > self.num_threads and self.num_threads > 1

        if use_parallel:
            # === PARALLEL BIN PROCESSING ===
            self.run.info_single(f"Processing {num_bins} bins in parallel using {self.num_threads} threads "
                                f"(homology mode, each BLAST uses 1 thread).", nl_before=1)

            # Prepare config dict for workers
            config = {
                'temp_dir': self.temp_dir,
                'word_size': self.word_size,
            }

            # Convert rt_windows to serializable list format with window_id
            rt_windows_list = []
            for window_id, window_info in self.rt_windows.items():
                window_dict = dict(window_info)
                window_dict['window_id'] = window_id
                rt_windows_list.append(window_dict)

            # Setup queues
            manager = multiprocessing.Manager()
            input_queue = manager.Queue()
            output_queue = manager.Queue()

            # Put all bins in input queue
            for bin_name, bin_contigs_list in bin_contigs_dict.items():
                input_queue.put((bin_name, bin_contigs_list))

            # Start workers
            workers = []
            for i in range(self.num_threads):
                worker = multiprocessing.Process(
                    target=DGR_Finder.process_bin_worker_homology,
                    args=(input_queue, output_queue, config, rt_windows_list, self.contigs_db_path))
                workers.append(worker)
                worker.start()

            # Monitor progress
            self.progress.new('Processing bins for homology BLAST (parallel)', progress_total_items=num_bins)
            self.progress.update(f"Processing {num_bins} bins across {self.num_threads} workers...")

            bins_processed = 0
            while bins_processed < num_bins:
                try:
                    bin_name, success, blast_path, error_msg = output_queue.get()
                    bins_processed += 1
                    self.progress.increment(increment_to=bins_processed)

                    if success:
                        successful_bins.append(bin_name)
                        blast_outputs[bin_name] = blast_path
                    else:
                        skipped_bins[bin_name] = error_msg

                    if bins_processed < num_bins:
                        self.progress.update(f"Processed {bins_processed}/{num_bins} bins...")
                    else:
                        self.progress.update("All done!")

                except KeyboardInterrupt:
                    self.run.info_single("Received kill signal, terminating workers...", nl_before=1)
                    break

            self.progress.end()

            # Terminate workers
            for worker in workers:
                worker.terminate()

        else:
            # === SEQUENTIAL BIN PROCESSING ===
            if self.num_threads > 1:
                self.run.info_single(f"Processing {num_bins} bins sequentially "
                                    f"(homology mode, each BLAST uses {self.num_threads} threads).", nl_before=1)

            self.progress.new('Processing bins for homology BLAST', progress_total_items=num_bins)

            # Load contig sequences once
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
            contigs_db.disconnect()

            for bin_name, bin_contigs_list in bin_contigs_dict.items():
                self.progress.update(f"{bin_name}", increment=True)

                try:
                    # Filter contig sequences to this bin
                    bin_contig_sequences = {contig: contig_sequences[contig]
                                           for contig in bin_contigs_list if contig in contig_sequences}

                    # Run homology mode BLAST for this bin
                    blast_output_path = self.run_blast_homology_mode(
                        contig_sequences=bin_contig_sequences,
                        bin_name=bin_name,
                        bin_contigs=bin_contigs_list
                    )

                    if blast_output_path:
                        blast_outputs[bin_name] = blast_output_path
                        successful_bins.append(bin_name)
                    else:
                        skipped_bins[bin_name] = "No RT windows in bin"

                except ConfigError as e:
                    skipped_bins[bin_name] = str(e)
                    continue
                except Exception as e:
                    skipped_bins[bin_name] = str(e)
                    continue

            self.progress.end()

        # Summary message
        self.run.info_single(f"Homology mode processed {PL('bin', num_bins)}: {len(successful_bins)} successful, "
                            f"{len(skipped_bins)} skipped.", nl_before=1)

        # Report skipped bins if any
        if skipped_bins:
            self.run.warning(f"{PL('bin', len(skipped_bins))} skipped during homology BLAST: {', '.join(skipped_bins.keys())}. "
                           "Use '--debug' to see detailed error messages.",
                           header="BINS SKIPPED (HOMOLOGY)")
            if anvio.DEBUG:
                for bin_name, error_msg in skipped_bins.items():
                    self.run.info_single(f"  {bin_name}: {error_msg}")

        return blast_outputs


    def process_standard_mode(self):
        """
        Process all contigs in standard mode by running BLASTn once across all SNV clusters.

        Parameters
        ==========
        None

        Returns
        =======
        blast_output : str
            Path to the BLASTn output XML file.

        Raises
        ======
        ConfigError
            If no valid SNV clusters are found.
        """

        sample_id_list, contig_sequences = self.load_data_and_setup()
        contig_records = self.find_snv_clusters(sample_id_list, contig_sequences)
        self.blast_output = self.run_blast(contig_records, contig_sequences, mode='activity')

        return self.blast_output



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

        # extract all starts and stops
        all_start = []
        all_end = []
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



    def process_blast_results(self, max_percent_identity=100, vr_in_query=True, apply_snv_filters=True, mode='activity'):
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
        vr_in_query : bool
            If True (activity mode), query=VR, subject=TR.
            If False (homology mode), query=TR, subject=VR.
        apply_snv_filters : bool
            If True, perform SNV-based filtering. If False, skip SNV analysis.
        mode : str
            Detection mode ('activity' or 'homology'). Affects BLAST output file naming
            pattern in collections mode.

        Returns
        =======
        mismatch_hits : dict
            A dictionary of all of the BLASTn hits that are less than 100%

        Raises
        =======
        ConfigError
            If no  BLAST output then exit
        """

        # if collections_mode is enabled, process multiple bins separately
        if self.collections_mode:
            self.merged_mismatch_hits = defaultdict(lambda: defaultdict(dict))

            tmp_directory_path = self.temp_dir

            for bin_name in self.bin_names_list:

                # reset mismatch hits for each bin
                self.mismatch_hits = defaultdict(lambda: defaultdict(dict))

                # Use mode-specific file naming pattern
                if mode == 'homology':
                    blast_file = os.path.join(
                        tmp_directory_path,
                        f"blast_output_homology_for_bin_{bin_name}_wordsize_{self.word_size}.xml"
                    )
                else:
                    blast_file = os.path.join(
                        tmp_directory_path,
                        f"blast_output_for_bin_{bin_name}_wordsize_{self.word_size}.xml"
                    )

                if not os.path.exists(blast_file):
                    self.run.warning(f"Warning: BLAST output file for {bin_name} ({mode} mode) not found. Skipping...")
                    continue

                if os.stat(blast_file).st_size == 0:
                    self.run.warning(f"No DGR like sequences are being found via BLAST ({mode} mode).", header="NO DGRS FOUND")
                    raise ConfigError("Therefore, we will exit here because anvi'o has found no DGRs in your data, "
                                    "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                                    "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")

                # process mismatches from this bin
                self.parse_and_process_blast_results(blast_file, bin_name, max_percent_identity,
                                                     vr_in_query=vr_in_query, apply_snv_filters=apply_snv_filters)

                # merge results into `merged_mismatch_hits`
                for section_id, hits_dict in self.mismatch_hits.items():
                    for hit_id, hit_data in hits_dict.items():
                        self.merged_mismatch_hits[section_id][hit_id] = hit_data

            self.run.info_single(f"Total unique mismatches: {len(self.merged_mismatch_hits)}", nl_before=1)

            return self.merged_mismatch_hits

        # run in normal none collections mode
        else:
            # single BLAST output processing
            self.mismatch_hits = {}

            if os.stat(f"{self.blast_output}").st_size == 0:
                self.run.warning("No DGR like sequences are being found via BLAST.", header="NO DGRS FOUND")
                raise ConfigError("Therefore, we will exit here because anvi'o has found no DGRs in your data, "
                                "nada, nowt, nothin'! However, you can go back and tinker with the parameters "
                                "of this tool if you believe this should not be the case. Anvi'o wishes you a nice day :)")

            # process the BLAST output normally
            self.parse_and_process_blast_results(self.blast_output, bin_name=None, max_percent_identity=max_percent_identity,
                                                 vr_in_query=vr_in_query, apply_snv_filters=apply_snv_filters)

            return self.mismatch_hits



    def has_repeat(self, seq: str, qseq: str, hseq: str) -> bool:
        """
        Check whether a sequence contains tandem repeats or STRs.
        Returns True if a repeat is found, otherwise False.

        Both the sequence and its reverse complement are checked, because
        depending on the dominant base (A vs T), the output may be reverse
        complemented. We need to ensure neither orientation exceeds the
        repeat threshold.
        """

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

        # only build a new string if needed
        if "-" in seq:
            seq = seq.replace("-", "")

        # guard against empty sequence after removing dashes
        if not seq:
            return False

        # using pytantan (https://doi.org/10.1093/nar/gkq1212) as a repeat finder
        # historically we used pytrf (https://doi.org/10.1186/s12859-025-06168-3)
        # citation notice is shown once at the start of parse_and_process_blast_results()

        # Check both the sequence and its reverse complement, since the output
        # may be reverse complemented if the dominant base is T
        seq_rc = utils.rev_comp(seq)

        for seq_to_check in [seq, seq_rc]:
            masked_seq = pytantan.mask_repeats(seq_to_check)
            # Count lowercase (masked) characters - str.count() is C-optimized
            num_masked = sum(masked_seq.count(c) for c in 'acgtn')
            frac_masked = num_masked / len(masked_seq)
            if frac_masked > self.repeat_threshold:
                if anvio.DEBUG and self.verbose:
                    self.run.warning(f"Removing the candidate DGR with its VR on this contig: {qseq} and TR on this contig: {hseq}. "
                                     f"This is because of repeats found in the sequence by pytantan ({frac_masked:.1%} masked).")
                return True

        # for atr1 in pytrf.ATRFinder('name', seq):
        #     if int(atr1.repeat) > self.numb_imperfect_tandem_repeats:
        #         #has_repeat = True
        #         if anvio.DEBUG and self.verbose:
        #             self.run.warning(f"Removing the DGR with its VR on this contig: {qseq} and TR on this contig: {hseq}. Found imperfect tandem repeat in the query sequence {seq} with repeat count {atr1.repeat} and motif {atr1.motif}")
        #         return True

        # #look for tandem homopolymers and 2 base short tandem repeats that are over 4 (so occur 5 times) times in the sequence
        # for ssr in pytrf.STRFinder('name', seq, mono=10, di=5):
        #     #if ((len(ssr.motif) == 1) or (len(ssr.motif) == 2)) and ssr.repeat > 6:
        #     if ssr:
        #         #has_repeat = True
        #         if anvio.DEBUG and self.verbose:
        #             self.run.warning(f"Removing the DGR with its VR on this contig: {qseq} and TR on this contig: {hseq}. Found tandem repeat in the query sequence {seq} with repeat count {ssr.repeat} and motif {ssr.motif}")
        #         return True

        # #look for approximate tandem repeats that in the VR, using a coverage value of the motif length times by the number of repeats
        # # divided by the sequence length
        # for atr in pytrf.ATRFinder('name', seq, min_motif=4, max_motif=10, min_seedrep=2, min_identity=70):
        #     coverage = (len(atr.motif)*atr.repeat) / len(seq)
        #     if coverage > self.repeat_motif_coverage:
        #         if anvio.DEBUG and self.verbose:
        #             self.run.warning(f"Removing the DGR with its VR on this contig: {qseq} and TR on this contig: {hseq}.Found approximate tandem repeat in the query sequence {seq} with repeat count {atr.repeat}, motif {atr.motif} and coverage {coverage}")
        #         return True

        # # look for approximate tandem repeats that in the VR, using a coverage value of the motif length times by the number of repeats
        # # divided by the sequence length
        # # need to do this for shorter motifs too because pytrf misses them otherwise
        # for atr2 in pytrf.ATRFinder('name', seq, min_motif=3, max_motif=10, min_seedrep=2, min_seedlen=9, min_identity=70):
        #     coverage = (len(atr2.motif)*atr2.repeat) / len(seq)
        #     if coverage > self.repeat_motif_coverage:
        #         if anvio.DEBUG and self.verbose:
        #             self.run.warning(f"Removing the DGR with its VR on this contig: {qseq} and TR on this contig: {hseq}. Found approximate tandem repeat in the query sequence {seq} with repeat count {atr2.repeat}, motif {atr2.motif} and coverage {coverage}")
        #         return True

        return False

    def find_optimal_mismatch_window(self, vr_seq, tr_seq, chars_to_skip, min_length, min_mismatches, max_non_dominant=1, max_gaps=0):
        """
        Find the longest contiguous window where:
        1. First and last mismatches are to the dominant base
        2. At most max_non_dominant mismatches are to non-dominant bases

        The window extends through matching bases, stopping just before
        any mismatch outside the valid window.

        Parameters
        ==========
        vr_seq : str
            Variable region sequence from BLAST alignment.
        tr_seq : str
            Template region sequence from BLAST alignment (used for dominant base detection).
        chars_to_skip : set
            Characters to ignore when counting mismatches (e.g., {'N', '-'}).
        min_length : int
            Minimum length of the trimmed alignment.
        min_mismatches : int
            Minimum number of mismatches required in the trimmed window.
        max_non_dominant : int
            Maximum number of non-dominant mismatches allowed in window (default 1).
        max_gaps : int
            Maximum number of gap characters ('-') allowed in the trimmed alignment window.

        Returns
        =======
        tuple or None
            (trim_start, trim_end, dominant_base) if valid window found, None otherwise.
        """
        # Collect all mismatch positions and their TR base
        mismatch_positions = []
        for idx in range(len(vr_seq)):
            if vr_seq[idx] in chars_to_skip or tr_seq[idx] in chars_to_skip:
                continue
            if vr_seq[idx] != tr_seq[idx]:
                mismatch_positions.append((idx, tr_seq[idx]))

        if len(mismatch_positions) < min_mismatches:
            return None

        # Find the dominant base (most frequent on TR side)
        base_counts = defaultdict(int)
        for _, tr_base in mismatch_positions:
            base_counts[tr_base] += 1

        dominant_base = max(base_counts, key=base_counts.get)

        # Mark each mismatch as dominant or not
        is_dominant = [tr_base == dominant_base for _, tr_base in mismatch_positions]

        n_mismatches = len(mismatch_positions)

        # Find longest valid window [i, j) where:
        # - mismatch i is dominant (first in window)
        # - mismatch j-1 is dominant (last in window)
        # - at most max_non_dominant non-dominant mismatches in window
        best_start = None
        best_end = None
        best_length = 0

        for i in range(n_mismatches):
            # First mismatch in window must be dominant
            if not is_dominant[i]:
                continue

            for j in range(i + 1, n_mismatches + 1):
                # Last mismatch in window must be dominant
                if not is_dominant[j - 1]:
                    continue

                # Count non-dominant mismatches in this window
                non_dominant_count = sum(1 for k in range(i, j) if not is_dominant[k])

                if non_dominant_count > max_non_dominant:
                    continue

                # Valid window - compute alignment coordinates with extension

                # Extend backward: start at beginning or after the mismatch before i
                if i > 0:
                    align_start = mismatch_positions[i - 1][0] + 1
                else:
                    align_start = 0

                # Extend forward: end at sequence end or before the mismatch after j-1
                if j < n_mismatches:
                    align_end = mismatch_positions[j][0]
                else:
                    align_end = len(vr_seq)

                # Check gap count in the extended window
                gap_count = vr_seq[align_start:align_end].count('-') + tr_seq[align_start:align_end].count('-')
                if gap_count > max_gaps:
                    continue

                window_length = align_end - align_start
                mismatch_count = j - i

                # Check minimum requirements and track best
                if window_length >= min_length and mismatch_count >= min_mismatches and window_length > best_length:
                    best_start = align_start
                    best_end = align_end
                    best_length = window_length

        if best_start is not None:
            return (best_start, best_end, dominant_base)
        return None


    def parse_and_process_blast_results(self, xml_file_path, bin_name, max_percent_identity,
                                         vr_in_query=True, apply_snv_filters=True):
        """
        Parse and process BLAST XML results for a single bin or dataset.

        This function filters BLAST HSPs for TR/VR alignment characteristics:
        - Mismatched hits below 100% identity
        - Filters out highly repeated sequences
        - Requires mismatches to be biased toward one base type (typically A on TR side)
        - Checks sequence diversity requirements

        The function supports two modes based on vr_in_query:
        - Activity mode (vr_in_query=True): Query=VR candidates, Subject=genome (contains TR)
        - Homology mode (vr_in_query=False): Query=TR candidates (RT windows), Subject=genome (contains VR)

        Parameters
        ==========
        xml_file_path : str
            Path to BLAST XML output file.
        bin_name : str or None
            Name of the bin associated with the BLAST search. If None,
            results are treated as coming from a single dataset.
        max_percent_identity : float
            Maximum percent identity allowed for considering a hit as mismatched.
            Hits with higher identity are ignored.
        vr_in_query : bool
            If True (default, activity mode), query sequences are VR candidates and
            subjects are potential TRs. If False (homology mode), query sequences are
            TR candidates (RT windows) and subjects are potential VRs.
        apply_snv_filters : bool
            If True (default), perform SNV-based analysis and confidence scoring.
            If False (homology mode), skip SNV analysis and set confidence to "homology-based".

        Returns
        =======
        None
            Updates self.mismatch_hits dictionary with filtered results.
        """

        hit_id_counter = 0
        current_section_id = None
        current_subject_contig = None

        if not hasattr(self, 'mismatch_hits') or not isinstance(self.mismatch_hits, defaultdict):
            self.mismatch_hits = defaultdict(lambda: defaultdict(dict))

        # Show pytantan citation once before processing
        self.run.warning("Anvi'o will now review the candidate DGRs for repeated sequences with the "
                        "python wrapped tantan repeat finder. DOI: https://doi.org/10.1093/nar/gkq1212",
                        lc='green', header="CITATION")

        # === PRE-INDEX SNV DATA BY CONTIG AS SORTED NUMPY ARRAYS ===
        # Only needed for activity-based detection (when apply_snv_filters=True)
        snv_index = {}
        if apply_snv_filters and hasattr(self, 'snv_panda') and self.snv_panda is not None:
            # Using sorted arrays enables O(log n) binary search for range queries
            # instead of O(n) pandas filtering per HSP
            for contig, group in self.snv_panda.groupby('contig_name'):
                sorted_group = group.sort_values('pos_in_contig')
                snv_index[contig] = {
                    'positions': sorted_group['pos_in_contig'].values,
                    'codon_pos': sorted_group['base_pos_in_codon'].values,
                    'reference': sorted_group['reference'].values
                }

        chars_to_skip = []
        if self.skip_Ns:
            chars_to_skip.append('N')
        if self.skip_dashes:
            chars_to_skip.append('-')

        # iterate over XML HSPs one by one
        # use iterparse with start/end events for better context tracking
        context = ET.iterparse(xml_file_path, events=("start", "end"))
        try:
            for event,elem in context:
                if event == "end":
        #for event, elem in ET.iterparse(self.blast_output, events=("end",)):

                    # track the current iteration (query)
                    if elem.tag == "Iteration_query-def":
                        current_section_id = elem.text
                        elem.clear()
                        continue

                    # track the current hit (subject)
                    if elem.tag == "Hit_def":
                        current_subject_contig = elem.text
                        elem.clear()
                        continue

                    if elem.tag != "Hsp":
                        continue

                    if elem.tag == "Hsp":
                        identical_positions = int(elem.find('Hsp_identity').text)
                        alignment_length = int(elem.find('Hsp_align-len').text)
                        percentage_identity = (identical_positions / alignment_length) * 100

                        if percentage_identity >= max_percent_identity:
                            elem.clear()
                            continue

                        # get parent Iteration and Hit elements
                        section_id = current_section_id

                        # extract start position from section_id (using pre-compiled regex)
                        match = SECTION_ID_PATTERN.search(section_id)
                        query_start_position = int(match.group(1)) if match else 0

                        # Stage 2: Parse sequences for mismatch analysis
                        qseq = elem.find('Hsp_qseq').text
                        hseq = elem.find('Hsp_hseq').text

                        query_mismatch_positions = []
                        # use defaultdict for lazy initialization - only creates entries for actual mismatches
                        query_mismatch_counts = defaultdict(int)
                        subject_mismatch_counts = defaultdict(int)

                        for idx in range(len(qseq)):
                            if qseq[idx] in chars_to_skip or hseq[idx] in chars_to_skip:
                                continue
                            if qseq[idx] != hseq[idx]:
                                query_mismatch_counts[qseq[idx]] += 1
                                query_mismatch_positions.append(idx)
                                subject_mismatch_counts[hseq[idx]] += 1

                        # get number of mismatches
                        mismatch_length_bp = len(query_mismatch_positions)

                        # if num of mismatches = 0, skip DGR search sanity check
                        if mismatch_length_bp == 0:
                            continue

                        # ========== TWO-STAGE FILTERING ==========
                        # Stage 1: Permissive initial filter
                        # Check if any base has >= initial_threshold (60%) of mismatches
                        # IMPORTANT: We check the TR side for dominant base (A)
                        # In activity mode: TR is subject (hseq), so use subject_mismatch_counts
                        # In homology mode: TR is query (qseq), so use query_mismatch_counts
                        tr_mismatch_counts = subject_mismatch_counts if vr_in_query else query_mismatch_counts
                        passes_initial_filter = False
                        for letter, count in tr_mismatch_counts.items():
                            initial_percentage = count / mismatch_length_bp
                            if initial_percentage >= self.initial_mismatch_bias_threshold:
                                passes_initial_filter = True
                                break

                        if not passes_initial_filter or mismatch_length_bp < self.number_of_mismatches:
                            continue

                        # Stage 2: Optimal window trimming
                        # Find the longest region where >= 95% of mismatches are to dominant base
                        # The function expects (vr_seq, tr_seq) so we swap order for homology mode
                        if vr_in_query:
                            # Activity mode: query=VR, subject=TR
                            trim_result = self.find_optimal_mismatch_window(
                                qseq, hseq, chars_to_skip,
                                self.minimum_vr_length,
                                self.number_of_mismatches,
                                self.max_non_dominant,
                                self.max_alignment_gaps
                            )
                        else:
                            # Homology mode: query=TR, subject=VR - swap order
                            trim_result = self.find_optimal_mismatch_window(
                                hseq, qseq, chars_to_skip,
                                self.minimum_vr_length,
                                self.number_of_mismatches,
                                self.max_non_dominant,
                                self.max_alignment_gaps
                            )

                        if trim_result is None:
                            # No valid trimmed window found
                            continue

                        trim_start, trim_end, dominant_base = trim_result

                        # Trim sequences to optimal window
                        trimmed_qseq = qseq[trim_start:trim_end]
                        trimmed_hseq = hseq[trim_start:trim_end]

                        # Stage 4: Repeat check on trimmed sequences only
                        if self.has_repeat(trimmed_qseq, trimmed_qseq, trimmed_hseq) or \
                           self.has_repeat(trimmed_hseq, trimmed_qseq, trimmed_hseq):
                            elem.clear()
                            continue

                        # Recalculate mismatch counts for trimmed region
                        query_mismatch_counts = defaultdict(int)
                        subject_mismatch_counts = defaultdict(int)
                        query_mismatch_positions = []

                        for idx in range(len(trimmed_qseq)):
                            if trimmed_qseq[idx] in chars_to_skip or trimmed_hseq[idx] in chars_to_skip:
                                continue
                            if trimmed_qseq[idx] != trimmed_hseq[idx]:
                                query_mismatch_counts[trimmed_qseq[idx]] += 1
                                query_mismatch_positions.append(idx)
                                subject_mismatch_counts[trimmed_hseq[idx]] += 1

                        mismatch_length_bp = len(query_mismatch_positions)

                        # Use the dominant base from trimming (handles T -> A conversion later)
                        letter = dominant_base

                        # Recalculate percentage for the trimmed region
                        # The dominant base should be on the TR side
                        # Activity mode: TR is subject; Homology mode: TR is query
                        tr_counts = subject_mismatch_counts if vr_in_query else query_mismatch_counts
                        percentage_of_mismatches = tr_counts.get(letter, 0) / mismatch_length_bp if mismatch_length_bp > 0 else 0

                        # To test for VR diversity of base types in the sequence
                        # VR should have multiple different bases at mismatch positions
                        # Activity mode: VR is query; Homology mode: VR is subject
                        vr_counts = query_mismatch_counts if vr_in_query else subject_mismatch_counts
                        non_zero_bases = sum(1 for count in vr_counts.values() if count > 0)
                        if not non_zero_bases >= self.min_mismatching_base_types_vr:
                            continue

                        # Stage 6-7: Parse remaining XML elements (deferred until after initial filters pass)
                        hit_from = int(elem.find('Hsp_hit-from').text)
                        hit_to = int(elem.find('Hsp_hit-to').text)
                        query_from = int(elem.find('Hsp_query-from').text)
                        query_to = int(elem.find('Hsp_query-to').text)
                        query_frame = int(elem.find('Hsp_query-frame').text)
                        subject_frame = int(elem.find('Hsp_hit-frame').text)
                        midline = elem.find('Hsp_midline').text

                        # Determine strand orientation from BLAST coordinates
                        # When hit_from > hit_to, the subject hit is on the reverse strand
                        subject_on_reverse_strand = hit_from > hit_to
                        query_on_reverse_strand = query_from > query_to

                        # Compute base genome positions (0-based)
                        # For query: always forward strand relative to our input sequences
                        query_genome_start_position = query_start_position + min(query_from - 1, query_to - 1)
                        query_genome_end_position = query_start_position + max(query_from - 1, query_to - 1)

                        # For subject: need to track strand for proper coordinate handling
                        subject_genome_start_position = min(hit_from - 1, hit_to - 1)
                        subject_genome_end_position = max(hit_from - 1, hit_to - 1)

                        # Update coordinates to reflect trimming
                        # trim_start/trim_end are alignment-relative positions
                        alignment_length = trim_end - trim_start

                        # Query trimming (query is always forward strand in our usage)
                        query_genome_start_position = query_genome_start_position + trim_start
                        query_genome_end_position = query_genome_start_position + alignment_length - 1

                        # Subject trimming depends on strand
                        if subject_on_reverse_strand:
                            # Reverse strand: alignment position 0 = genomic end
                            # Trimming from alignment start removes from genomic end
                            subject_genome_end_position = subject_genome_end_position - trim_start
                            subject_genome_start_position = subject_genome_end_position - alignment_length + 1
                        else:
                            # Forward strand: alignment position 0 = genomic start
                            subject_genome_start_position = subject_genome_start_position + trim_start
                            subject_genome_end_position = subject_genome_start_position + alignment_length - 1

                        # Trim midline to match trimmed sequences
                        trimmed_midline = midline[trim_start:trim_end]

                        # Now proceed with the trimmed sequences
                        original_query_frame = query_frame
                        original_subject_frame = subject_frame
                        original_midline = trimmed_midline

                        # if the letter is T, then we assume that it is an A base and we reverse complement EVERYTHING
                        if letter == 'T':
                            hit_sequence = str(Seq(trimmed_hseq).reverse_complement())
                            query_sequence = str(Seq(trimmed_qseq).reverse_complement())
                            midline = ''.join(reversed(original_midline))
                            base = 'A'
                            is_reverse_complement = True
                            # Flip frames for reverse complement
                            # (actual VR/TR assignment happens in semantic swap below)
                            new_query_frame = original_query_frame * -1
                            new_subject_frame = original_subject_frame * -1
                        else:
                            hit_sequence = str(trimmed_hseq)
                            query_sequence = str(trimmed_qseq)
                            midline = original_midline
                            base = letter
                            is_reverse_complement = False
                            new_query_frame = original_query_frame
                            new_subject_frame = original_subject_frame

                        # to test for VR diversity of base types in the sequence
                        # count the distinct base types in the sequence
                        # ignore gaps and ambiguous bases if needed
                        # Activity mode: VR is query_sequence; Homology mode: VR is hit_sequence
                        vr_sequence_for_diversity = query_sequence if vr_in_query else hit_sequence
                        vr_unique_bases = set(vr_sequence_for_diversity) - {"-", "N"}

                        # ensure the sequence has at least the required number of distinct base types (default 2)
                        if len(vr_unique_bases) <= self.min_base_types_vr:
                            continue

                        # By default (conservative), filter for A-base mutagenesis only.
                        if not self.allow_any_base and base != "A":
                            continue

                        # Extract contig names from BLAST query/subject
                        query_contig = section_id.split('_section', 1)[0]
                        subject_contig = current_subject_contig

                        # ========== ASSIGN VR/TR BASED ON DETECTION MODE ==========
                        # In activity mode (vr_in_query=True): query=VR, subject=TR
                        # In homology mode (vr_in_query=False): query=TR (RT window), subject=VR
                        if vr_in_query:
                            # Activity mode: query is VR, subject is TR
                            vr_contig = query_contig
                            tr_contig = subject_contig
                            vr_sequence = query_sequence
                            tr_sequence = hit_sequence
                            vr_start = query_genome_start_position
                            vr_end = query_genome_end_position
                            tr_start = subject_genome_start_position
                            tr_end = subject_genome_end_position
                            vr_frame = new_query_frame
                            tr_frame = new_subject_frame
                            vr_mismatch_counts = query_mismatch_counts
                            tr_mismatch_counts = subject_mismatch_counts
                        else:
                            # Homology mode: query is TR (RT window), subject is VR
                            vr_contig = subject_contig
                            tr_contig = query_contig
                            vr_sequence = hit_sequence
                            tr_sequence = query_sequence
                            vr_start = subject_genome_start_position
                            vr_end = subject_genome_end_position
                            tr_start = query_genome_start_position
                            tr_end = query_genome_end_position
                            vr_frame = new_subject_frame
                            tr_frame = new_query_frame
                            vr_mismatch_counts = subject_mismatch_counts
                            tr_mismatch_counts = query_mismatch_counts

                        # Mismatch positions relative to VR contig
                        # query_mismatch_positions are alignment-relative (0-indexed within trimmed alignment)
                        # Need to convert to genomic coordinates, accounting for strand
                        if vr_in_query:
                            # Activity mode: VR is query (always forward strand in our usage)
                            mismatch_pos_contig_relative = [x + vr_start for x in query_mismatch_positions]
                        else:
                            # Homology mode: VR is subject (can be forward or reverse strand)
                            if subject_on_reverse_strand:
                                # Reverse strand: alignment position 0 = genomic end (vr_end)
                                # Position X in alignment = genomic position (vr_end - X)
                                mismatch_pos_contig_relative = [vr_end - x for x in query_mismatch_positions]
                            else:
                                # Forward strand: alignment position 0 = genomic start (vr_start)
                                mismatch_pos_contig_relative = [x + vr_start for x in query_mismatch_positions]

                        # ========== SNV ANALYSIS (only for activity mode) ==========
                        if apply_snv_filters:
                            # Activity mode: perform full SNV analysis
                            # subset snv by VR contig and VR range using binary search
                            contig_data = snv_index.get(vr_contig)
                            if contig_data is not None:
                                positions = contig_data['positions']
                                left_idx = bisect.bisect_left(positions, vr_start)
                                right_idx = bisect.bisect_right(positions, vr_end)
                                snv_positions = positions[left_idx:right_idx]
                            else:
                                snv_positions = np.array([], dtype=int)

                            snv_VR_positions = list(dict.fromkeys(snv_positions))

                            # Analyze mismatch codon positions (uses gene position data)
                            mismatch_analysis = self.analyze_mismatch_codon_positions(
                                mismatch_pos_contig_relative,
                                vr_contig
                            )

                            # Analyze SNV pattern
                            snv_analysis = self.analyze_snv_pattern(
                                vr_contig,
                                vr_start,
                                vr_end,
                                mismatch_pos_contig_relative,
                                base,
                                is_reverse_complement
                            )

                            # Compute confidence level
                            confidence, confidence_reasons = self.compute_dgr_confidence(
                                mismatch_analysis,
                                snv_analysis
                            )

                            numb_of_SNVs = snv_analysis['n_snvs_total'] if snv_analysis else 0
                        else:
                            # Homology mode: skip SNV analysis, use default values
                            snv_VR_positions = []

                            # Still analyze mismatch codon positions (doesn't require SNV data)
                            mismatch_analysis = self.analyze_mismatch_codon_positions(
                                mismatch_pos_contig_relative,
                                vr_contig
                            )

                            snv_analysis = None
                            confidence = "homology-based"
                            confidence_reasons = []
                            numb_of_SNVs = 0

                        # Basic counts
                        numb_of_mismatches = len(query_mismatch_positions)

                        # ========== BUILD HIT DATA ==========
                        # Use semantic VR/TR naming regardless of BLAST query/subject direction
                        hit_data = {
                            'bin': bin_name if bin_name else "N/A",
                            'query_section': section_id,
                            # Semantic VR/TR fields
                            'query_seq': vr_sequence,  # VR sequence (for backwards compat, query_seq = VR)
                            'hit_seq': tr_sequence,    # TR sequence (for backwards compat, hit_seq = TR)
                            'midline': midline,
                            'query_contig': vr_contig,  # VR contig
                            'subject_contig': tr_contig,  # TR contig
                            'subject_genome_start_position': tr_start,  # TR start
                            'subject_genome_end_position': tr_end,  # TR end
                            'query_mismatch_counts': vr_mismatch_counts,
                            'subject_mismatch_counts': tr_mismatch_counts,
                            'position': query_mismatch_positions,
                            'alignment_length': alignment_length,
                            'query_genome_start_position': vr_start,  # VR start
                            'query_genome_end_position': vr_end,  # VR end
                            'query_frame': vr_frame,
                            'subject_frame': tr_frame,
                            'base': base,
                            'is_reverse_complement': is_reverse_complement,
                            'numb_of_mismatches': numb_of_mismatches,
                            'numb_of_SNVs': numb_of_SNVs,
                            'percentage_of_mismatches': percentage_of_mismatches,
                            'mismatch_pos_contig_relative': mismatch_pos_contig_relative,
                            'snv_VR_positions': snv_analysis['snv_positions'] if snv_analysis else snv_VR_positions,
                            # Mismatch codon analysis
                            'mismatch_codon_1': mismatch_analysis['mismatch_codon_1'],
                            'mismatch_codon_2': mismatch_analysis['mismatch_codon_2'],
                            'mismatch_codon_3': mismatch_analysis['mismatch_codon_3'],
                            'pct_mismatch_codon_3': mismatch_analysis['pct_mismatch_codon_3'],
                            'vr_gene_id': mismatch_analysis['vr_gene_id'],
                            # SNV analysis (defaults for homology mode)
                            'n_snvs_total': snv_analysis['n_snvs_total'] if snv_analysis else 0,
                            'n_snvs_explained': snv_analysis['n_snvs_explained'] if snv_analysis else 0,
                            'n_snvs_unexplained': snv_analysis['n_snvs_unexplained'] if snv_analysis else 0,
                            'n_explained_diverse': snv_analysis['n_explained_diverse'] if snv_analysis else 0,
                            'pct_snvs_explained': snv_analysis['pct_snvs_explained'] if snv_analysis else 100,
                            'snv_codon_1': snv_analysis['snv_codon_1'] if snv_analysis else 0,
                            'snv_codon_2': snv_analysis['snv_codon_2'] if snv_analysis else 0,
                            'snv_codon_3': snv_analysis['snv_codon_3'] if snv_analysis else 0,
                            'pct_snv_codon_3': snv_analysis['pct_snv_codon_3'] if snv_analysis else 0,
                            # Per-sample SNV analysis results
                            'snv_supporting_sample': snv_analysis.get('snv_supporting_sample') if snv_analysis else None,
                            'sample_passed': snv_analysis.get('sample_passed', False) if snv_analysis else False,
                            # Confidence
                            'confidence': confidence,
                            'confidence_reasons': confidence_reasons,
                            # Backwards compatibility fields
                            'numb_snvs_in_3rd_codon_pos': snv_analysis['snv_codon_3'] if snv_analysis else 0,
                            'numb_of_snv_in_matches_not_mutagen_base': snv_analysis['n_snvs_unexplained'] if snv_analysis else 0,
                            'DGR_looks_false': confidence == 'low' if apply_snv_filters else False,
                            'snv_at_3_codon_over_a_third': mismatch_analysis['pct_mismatch_codon_3'] > 33,
                            # Detection method tracking
                            'detection_method': 'activity' if apply_snv_filters else 'homology',
                        }

                        # Filter out low confidence hits in activity mode
                        # We only keep medium and high confidence for activity-based detection
                        if apply_snv_filters and confidence == 'low':
                            elem.clear()
                            continue

                        # Stage 7: Only now increment counter and create unique ID for hits that passed all filters
                        hit_id_counter += 1
                        hit_identity_unique = f"{section_id}_count_{hit_id_counter}"

                        # store grouped by query_section (VR region)
                        self.mismatch_hits[section_id][hit_identity_unique] = hit_data

                        # free memory for this element
                        elem.clear()

                    # clear other elements we don't need
                    elif elem.tag not in ("Iteration", "Hit", "BlastOutput"):
                        elem.clear()

        finally:
            # ensure cleanup
            del context


    def filter_for_best_VR_TR(self):
        #NOTE: need to check this code to clean up and maybe put into one function to remove redundancy - Iva offered to help :)
        # We did this on 21.05.2024, created functions add_new_DGR and update_existing_DGR
        #NOTE: the arguments for add_new_DGR and update_existing_DGR that are based on the query and subject genome are positional, so they need
        # to be in the same place as those function deal with whether these are TR or VR. BUT the frame argument is not positional and you need
        # to change that based on whether it is in the query or subject. YES it is gross, but if it works, it works.

        """
        This function takes the hits of the BLASTn and chooses one singular hit that qualifies as a TR VR pair based on the filters for template and variable regions. So that each VR has one TR.

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
        bins_with_dgrs = set()

        # possible DGR dictionary
        self.DGRs_found_dict = {}

        # Display warning if running in permissive mode
        if self.allow_any_base:
                self.run.warning("You are using --allow-any-base, which reports DGRs with any dominant base (not just A). "
                                "This may include non-canonical DGRs or false positives. Review results carefully.",
                                header="PERMISSIVE MODE: Allowing any dominant base")

        if self.collections_mode:
            hits_item = self.merged_mismatch_hits
        else:
            hits_item = self.mismatch_hits

        for section_id, hits_dict in hits_item.items():
            # ensure hits_dict is a mapping of hit_id -> hit_data
            if isinstance(hits_dict, list):
                hits_dict = {f"hit_{i+1}": h for i, h in enumerate(hits_dict)}

            # create grouping dict for this section
            hits_by_query = defaultdict(list)

            # populate hits_by_query
            for hit_id, hit_data in hits_dict.items():
                # defensive: skip non-dicts
                if not isinstance(hit_data, dict):
                    continue
                query_section = hit_data.get('query_section')
                if query_section is None:
                    # if no query_section, skip or decide appropriate fallback
                    continue
                hits_by_query[query_section].append(hit_data)

            # now filter per query_section
            for query_section, query_hits in hits_by_query.items():

                # if only one hit, keep it
                if len(query_hits) == 1:
                    best_hit = query_hits[0]
                    best_amongst_multiple_TRs_for_one_VR = False
                else:
                    best_amongst_multiple_TRs_for_one_VR = True
                    # Default: prefer A-based hits when multiple TRs match one VR
                    if not self.allow_any_base:
                        a_hits = [h for h in query_hits if h.get('base') == 'A']
                        if a_hits:
                            query_hits = a_hits

                    # sort hits by your ranking criteria
                    # (negative for longest alignment first, others ascending)
                    query_hits.sort(
                        key=lambda h: (
                            -h.get('alignment_length', 0),
                            h.get('numb_of_SNVs', float('inf')),
                            h.get('numb_of_snv_in_matches_not_mutagen_base', float('inf')),
                            h.get('numb_snvs_in_3rd_codon_pos', float('inf')),
                        )
                    )

                best_hit = query_hits[0]

                # unpack dict
                query_section = best_hit['query_section']
                bin = best_hit['bin']
                subject_genome_start_position = best_hit['subject_genome_start_position']
                subject_genome_end_position = best_hit['subject_genome_end_position']
                TR_sequence = Seq(best_hit['hit_seq'])
                midline = best_hit['midline']
                VR_sequence = Seq(best_hit['query_seq'])
                query_genome_start_position = best_hit['query_genome_start_position']
                query_genome_end_position = best_hit['query_genome_end_position']
                VR_frame = int(best_hit['query_frame'])
                TR_frame = int(best_hit['subject_frame'])
                query_contig = best_hit['query_contig']
                subject_contig = best_hit['subject_contig']
                is_reverse_complement = best_hit['is_reverse_complement']
                base = best_hit['base']
                numb_of_snv_in_matches_not_mutagen_base= best_hit['numb_of_snv_in_matches_not_mutagen_base']
                numb_of_SNVs= best_hit['numb_of_SNVs']
                DGR_looks_snv_false = best_hit['DGR_looks_false']
                snv_at_3_codon_over_a_third = best_hit['snv_at_3_codon_over_a_third']
                numb_of_mismatches = best_hit['numb_of_mismatches']
                percentage_of_mismatches = best_hit['percentage_of_mismatches']
                mismatch_pos_contig_relative = best_hit['mismatch_pos_contig_relative']
                snv_VR_positions = best_hit['snv_VR_positions']
                # New codon/SNV analysis fields
                mismatch_codon_1 = best_hit.get('mismatch_codon_1', 0)
                mismatch_codon_2 = best_hit.get('mismatch_codon_2', 0)
                mismatch_codon_3 = best_hit.get('mismatch_codon_3', 0)
                pct_mismatch_codon_3 = best_hit.get('pct_mismatch_codon_3', 0)
                vr_gene_id = best_hit.get('vr_gene_id', None)
                n_snvs_total = best_hit.get('n_snvs_total', 0)
                n_snvs_explained = best_hit.get('n_snvs_explained', 0)
                n_snvs_unexplained = best_hit.get('n_snvs_unexplained', 0)
                n_explained_diverse = best_hit.get('n_explained_diverse', 0)
                pct_snvs_explained = best_hit.get('pct_snvs_explained', 100)
                snv_codon_1 = best_hit.get('snv_codon_1', 0)
                snv_codon_2 = best_hit.get('snv_codon_2', 0)
                snv_codon_3 = best_hit.get('snv_codon_3', 0)
                pct_snv_codon_3 = best_hit.get('pct_snv_codon_3', 0)
                snv_supporting_sample = best_hit.get('snv_supporting_sample', None)
                confidence = best_hit.get('confidence', 'N/A')
                confidence_reasons = best_hit.get('confidence_reasons', [])
                detection_method = best_hit.get('detection_method', 'activity')

                # need to check if the new TR you're looping through exists in the DGR_found_dict, see if position overlap
                if not self.DGRs_found_dict:
                    # add first DGR
                    num_DGR += 1
                    if anvio.DEBUG:
                        self.run.warning(f"Adding new DGR {num_DGR} in the bin: {bin}, the VR is on this contig: {query_contig}", header="NEW DGR", lc='yellow')
                    bins_with_dgrs.add(bin)
                    self.add_new_DGR(num_DGR, bin, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig,
                                base, is_reverse_complement, TR_frame, VR_sequence, VR_frame, subject_genome_start_position, subject_genome_end_position,
                                subject_contig, midline, percentage_of_mismatches, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative,
                                snv_VR_positions, numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs, best_amongst_multiple_TRs_for_one_VR,
                                mismatch_codon_1, mismatch_codon_2, mismatch_codon_3, pct_mismatch_codon_3, vr_gene_id,
                                n_snvs_total, n_snvs_explained, n_snvs_unexplained, n_explained_diverse, pct_snvs_explained,
                                snv_codon_1, snv_codon_2, snv_codon_3, pct_snv_codon_3, snv_supporting_sample, confidence, confidence_reasons,
                                detection_method)
                else:
                    was_added = False
                    for dgr in self.DGRs_found_dict:
                        if self.DGRs_found_dict[dgr]['TR_contig'] == subject_contig and self.range_overlapping(subject_genome_start_position,
                                                                                                        subject_genome_end_position,
                                                                                                        self.DGRs_found_dict[dgr]['TR_start_position'],
                                                                                                        self.DGRs_found_dict[dgr]['TR_end_position']):
                            was_added = True
                            #TODO can rename consensus_TR
                            self.update_existing_DGR(dgr, bin, TR_frame, VR_sequence, VR_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement, query_genome_start_position,
                                            query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                                            subject_contig, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative, snv_VR_positions,
                                            numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs, best_amongst_multiple_TRs_for_one_VR,
                                            mismatch_codon_1, mismatch_codon_2, mismatch_codon_3, pct_mismatch_codon_3, vr_gene_id,
                                            n_snvs_total, n_snvs_explained, n_snvs_unexplained, n_explained_diverse, pct_snvs_explained,
                                            snv_codon_1, snv_codon_2, snv_codon_3, pct_snv_codon_3, snv_supporting_sample, confidence, confidence_reasons,
                                            detection_method)
                            break
                    if not was_added:
                        # add new TR and its first VR
                        num_DGR += 1
                        if anvio.DEBUG:
                            self.run.warning(f"Adding new DGR {num_DGR} in the bin: {bin}, the VR is on this contig: {query_contig}", header="NEW DGR", lc='yellow')
                        bins_with_dgrs.add(bin)
                        self.add_new_DGR(num_DGR,
                                        bin,
                                        TR_sequence,
                                        query_genome_start_position,
                                        query_genome_end_position,
                                        query_contig,
                                        base,
                                        is_reverse_complement,
                                        TR_frame,
                                        VR_sequence,
                                        VR_frame,
                                        subject_genome_start_position, subject_genome_end_position,
                                        subject_contig,
                                        midline,
                                        percentage_of_mismatches,
                                        DGR_looks_snv_false,
                                        snv_at_3_codon_over_a_third, mismatch_pos_contig_relative,
                                        snv_VR_positions, numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches,
                                        numb_of_SNVs, best_amongst_multiple_TRs_for_one_VR,
                                        mismatch_codon_1, mismatch_codon_2, mismatch_codon_3, pct_mismatch_codon_3, vr_gene_id,
                                        n_snvs_total, n_snvs_explained, n_snvs_unexplained, n_explained_diverse, pct_snvs_explained,
                                        snv_codon_1, snv_codon_2, snv_codon_3, pct_snv_codon_3, snv_supporting_sample, confidence, confidence_reasons,
                                        detection_method)

        # summary of DGRs found
        if num_DGR > 0:
            if self.collections_mode:
                self.run.warning(f"Anvi'o found {PL('DGR', num_DGR)} across {PL('bin', len(bins_with_dgrs))}. "
                                "This is very exciting and you should celebrate.",
                                header="DGRs FOUND 🎉", lc='green')
            else:
                self.run.warning(f"Anvi'o found {PL('DGR', num_DGR)}. This is very exciting and you should celebrate.",
                                header="DGRs FOUND 🎉", lc='green')

        if anvio.DEBUG:
            self.run.warning(f"The temp directory, '{self.temp_dir}', is kept. Don't forget to clean it up later!", header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (use `--debug` to keep it for testing purposes)", nl_before=1)
            shutil.rmtree(self.temp_dir)
        return



    def add_new_DGR(self, DGR_number, bin, TR_sequence, query_genome_start_position, query_genome_end_position, query_contig, base,
                    is_reverse_complement, TR_frame, VR_sequence, VR_frame, subject_genome_start_position, subject_genome_end_position, subject_contig, midline,
                    percentage_of_mismatches, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative, snv_VR_positions,
                    numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs, best_amongst_multiple_TRs_for_one_VR,
                    mismatch_codon_1=0, mismatch_codon_2=0, mismatch_codon_3=0, pct_mismatch_codon_3=0, vr_gene_id=None,
                    n_snvs_total=0, n_snvs_explained=0, n_snvs_unexplained=0, n_explained_diverse=0, pct_snvs_explained=100,
                    snv_codon_1=0, snv_codon_2=0, snv_codon_3=0, pct_snv_codon_3=0,
                    snv_supporting_sample=None, confidence='N/A', confidence_reasons=None,
                    detection_method='activity'):
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

        Returns
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
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['best_amongst_multiple_TRs_for_one_VR'] = best_amongst_multiple_TRs_for_one_VR

        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_start_position'] = subject_genome_start_position
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_end_position'] = subject_genome_end_position
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_sequence'] = TR_sequence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_reverse_complement'] = is_reverse_complement
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['TR_frame'] = TR_frame
        # New codon/SNV analysis fields
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['mismatch_codon_1'] = mismatch_codon_1
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['mismatch_codon_2'] = mismatch_codon_2
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['mismatch_codon_3'] = mismatch_codon_3
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['pct_mismatch_codon_3'] = pct_mismatch_codon_3
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['vr_gene_id'] = vr_gene_id
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['n_snvs_total'] = n_snvs_total
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['n_snvs_explained'] = n_snvs_explained
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['n_snvs_unexplained'] = n_snvs_unexplained
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['n_explained_diverse'] = n_explained_diverse
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['pct_snvs_explained'] = pct_snvs_explained
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['snv_codon_1'] = snv_codon_1
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['snv_codon_2'] = snv_codon_2
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['snv_codon_3'] = snv_codon_3
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['pct_snv_codon_3'] = pct_snv_codon_3
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['snv_supporting_sample'] = snv_supporting_sample
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['confidence'] = confidence
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['confidence_reasons'] = confidence_reasons if confidence_reasons else []
        self.DGRs_found_dict[DGR_key]['VRs']['VR_001']['detection_method'] = detection_method


    def update_existing_DGR(self, existing_DGR_key, bin, TR_frame, VR_sequence, VR_frame, TR_sequence, midline, percentage_of_mismatches, is_reverse_complement,
                            query_genome_start_position, query_genome_end_position, query_contig, subject_genome_start_position, subject_genome_end_position,
                            subject_contig, DGR_looks_snv_false, snv_at_3_codon_over_a_third, mismatch_pos_contig_relative, snv_VR_positions,
                            numb_of_snv_in_matches_not_mutagen_base, numb_of_mismatches, numb_of_SNVs, best_amongst_multiple_TRs_for_one_VR,
                            mismatch_codon_1=0, mismatch_codon_2=0, mismatch_codon_3=0, pct_mismatch_codon_3=0, vr_gene_id=None,
                            n_snvs_total=0, n_snvs_explained=0, n_snvs_unexplained=0, n_explained_diverse=0, pct_snvs_explained=100,
                            snv_codon_1=0, snv_codon_2=0, snv_codon_3=0, pct_snv_codon_3=0,
                            snv_supporting_sample=None, confidence='N/A', confidence_reasons=None,
                            detection_method='activity'):
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

        Returns
        =======
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions and the corresponding info for those regions
        """

        if existing_DGR_key not in self.DGRs_found_dict:
            raise KeyError(f"Existing DGR key {existing_DGR_key} not found in DGRs_found_dict")

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
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['best_amongst_multiple_TRs_for_one_VR'] = best_amongst_multiple_TRs_for_one_VR
        # New codon/SNV analysis fields
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['mismatch_codon_1'] = mismatch_codon_1
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['mismatch_codon_2'] = mismatch_codon_2
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['mismatch_codon_3'] = mismatch_codon_3
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['pct_mismatch_codon_3'] = pct_mismatch_codon_3
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['vr_gene_id'] = vr_gene_id
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['n_snvs_total'] = n_snvs_total
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['n_snvs_explained'] = n_snvs_explained
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['n_snvs_unexplained'] = n_snvs_unexplained
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['n_explained_diverse'] = n_explained_diverse
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['pct_snvs_explained'] = pct_snvs_explained
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['snv_codon_1'] = snv_codon_1
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['snv_codon_2'] = snv_codon_2
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['snv_codon_3'] = snv_codon_3
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['pct_snv_codon_3'] = pct_snv_codon_3
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['snv_supporting_sample'] = snv_supporting_sample
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['confidence'] = confidence
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['confidence_reasons'] = confidence_reasons if confidence_reasons else []
        self.DGRs_found_dict[existing_DGR_key]['VRs'][new_VR_key]['detection_method'] = detection_method

        self.DGRs_found_dict[existing_DGR_key]['TR_start_position'] = min(subject_genome_start_position, self.DGRs_found_dict[existing_DGR_key]['TR_start_position'])
        self.DGRs_found_dict[existing_DGR_key]['TR_end_position'] = max(subject_genome_end_position, self.DGRs_found_dict[existing_DGR_key]['TR_end_position'])



    def get_gene_info(self):
        """
        This function collects information on the genes in which the variable regions act in.

        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======
        """

        if len(self.DGRs_found_dict) == 0:
            return

        self.run.info_single(f"Computing the Genes the Variable Regions occur in and creating a '{self.output_directory}_DGR_genes_found.tsv'.", nl_before=1)

        # initiate a dictionary for the gene where we find VR
        self.vr_gene_info = {}
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)

        # are there genes?
        if not contigs_db.meta['genes_are_called']:
            self.run.warning("There are no gene calls in your contigs database, therefore there is context to "
                            "learn about :/ Your reports will not include a file to study the genomic context "
                            "that surrounds the DGRs")

            contigs_db.disconnect()
            return

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

                        # while we are here, let's add more info about the gene
                        # DNA sequence:
                        where_clause = f'''contig == "{contig_name}"'''
                        contig_sequence = contigs_db.db.get_some_rows_from_table_as_dict(t.contig_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                        dna_sequence = contig_sequence[contig_name]['sequence'][gene_call['start']:gene_call['stop']]
                        if gene_call['direction'] == 'f':
                            gene_call['DNA_sequence'] = dna_sequence
                        else:
                            gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)

                        # add AA sequence
                        where_clause = f'''gene_callers_id == "{gene_callers_id}"'''
                        aa_sequence = contigs_db.db.get_some_rows_from_table_as_dict(t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                        gene_call['AA_sequence'] = aa_sequence[gene_callers_id]['sequence']

                        # gene length
                        gene_call['length'] = gene_call['stop'] - gene_call['start']

                        self.vr_gene_info[dgr][vr] = gene_call
                        break
        contigs_db.disconnect()

        # write the results to TSV file
        self.write_dgr_genes_found_tsv()

        return



    def write_dgr_genes_found_tsv(self):
        """
        Write the DGR gene information to a TSV file.

        Parameters
        ==========

        Returns
        =======
        {self.output_directory}_DGR_genes_found.tsv: tsv
            A TSV file containing the DGR gene information, including gene IDs, contig names etc
        """

        if not hasattr(self, 'vr_gene_info') or not self.vr_gene_info:
            self.run.warning("No gene information available to write to TSV file.")
            return

        output_directory_path = self.output_directory
        output_prefix = os.path.basename(self.output_directory)
        output_path_for_genes_found = os.path.join(output_directory_path, f"{output_prefix}_DGR_genes_found.tsv")

        # define the header for the TSV file
        csv_header = ['DGR_ID', 'VR_ID', 'Contig', 'Start', 'Stop', 'Direction', 'Partial', 'Call_Type', 'Gene_Caller_Source', 'Version', 'Gene_Caller_ID', 'DNA_Sequence', 'AA_Sequence', 'Length', 'Gene_Functions', 'Gene_Function_Source', 'Gene_Function_Accession']

        # open the TSV file in write mode
        with open(output_path_for_genes_found, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(csv_header)  # Write the header row
            # iterate through the dictionary and write each gene's information to the TSV file
            for dgr_id, vr_data in self.vr_gene_info.items():
                for vr_id, gene_info in vr_data.items():
                    if not gene_info:
                        continue

                    # convert list of functions to string
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

            self.run.info_single(f"DGR genes information successfully written to '{self.output_directory}'" , nl_before=1)
        return



    def get_rt_windows(self):
        """
        Find all RT HMM hits and create search windows around each for homology-based DGR detection.

        This method queries the contigs database for HMM hits matching the specified RT HMM sources,
        then creates genomic windows extending rt_window_size bp on each side of each RT gene.
        Windows are clipped to contig boundaries to avoid invalid coordinates.

        Returns
        =======
        dict
            Dictionary mapping window_id to window metadata:
            {
                'RT_window_0': {
                    'contig': str,           # Contig name
                    'window_start': int,     # Window start position (0-based)
                    'window_end': int,       # Window end position (0-based, inclusive)
                    'rt_gene_id': int,       # Gene callers ID of the RT
                    'rt_start': int,         # RT gene start position
                    'rt_end': int,           # RT gene end position
                    'rt_direction': str,     # RT gene direction ('f' or 'r')
                    'rt_gene_name': str,     # HMM gene name (e.g., 'Clean_DGR_clade_3')
                    'hmm_source': str,       # HMM source name
                    'e_value': float,        # HMM hit e-value
                },
                ...
            }
        """
        self.run.info_single("Finding RT HMM hits and creating search windows for homology-based detection...", nl_before=1)

        # Open contigs database and get required tables
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        hmm_hits_dict = contigs_db.db.get_table_as_dict(t.hmm_hits_table_name)
        genes_in_contigs = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
        contigs_db.disconnect()

        # Build contig length lookup
        contig_lengths = {contig: len(data['sequence']) for contig, data in contig_sequences.items()}

        # Find all RT HMM hits matching specified sources
        # For each gene_callers_id, keep the hit with lowest e-value
        rt_hits = {}
        for index, entry in hmm_hits_dict.items():
            if entry['source'] in self.hmm:
                gene_callers_id = entry['gene_callers_id']
                gene_name = entry['gene_name']
                e_value = entry['e_value']
                hmm_source = entry['source']

                if gene_callers_id not in rt_hits:
                    rt_hits[gene_callers_id] = {
                        'gene_name': gene_name,
                        'e_value': e_value,
                        'hmm_source': hmm_source
                    }
                elif e_value < rt_hits[gene_callers_id]['e_value']:
                    rt_hits[gene_callers_id]['gene_name'] = gene_name
                    rt_hits[gene_callers_id]['e_value'] = e_value
                    rt_hits[gene_callers_id]['hmm_source'] = hmm_source

        if not rt_hits:
            self.run.warning("No RT HMM hits were found in the contigs database. Homology-based detection "
                           "cannot proceed without RT genes. Please ensure you have run "
                           "`anvi-run-hmms -I Reverse_Transcriptase` on your contigs database.",
                           header="NO RT HMM HITS FOUND")
            return {}

        # Create windows around each RT hit
        rt_windows = {}
        window_counter = 0

        for gene_callers_id, hit_info in rt_hits.items():
            # Get gene position information
            if gene_callers_id not in genes_in_contigs:
                continue

            gene_info = genes_in_contigs[gene_callers_id]
            contig = gene_info['contig']
            rt_start = gene_info['start']
            rt_end = gene_info['stop']
            rt_direction = gene_info['direction']

            # Get contig length for boundary clipping
            contig_length = contig_lengths.get(contig, 0)
            if contig_length == 0:
                continue

            # Calculate window boundaries (clip to contig edges)
            window_start = max(0, rt_start - self.rt_window_size)
            window_end = min(contig_length - 1, rt_end + self.rt_window_size)

            # Create window entry
            window_id = f"RT_window_{window_counter}"
            rt_windows[window_id] = {
                'contig': contig,
                'window_start': window_start,
                'window_end': window_end,
                'rt_gene_id': gene_callers_id,
                'rt_start': rt_start,
                'rt_end': rt_end,
                'rt_direction': rt_direction,
                'rt_gene_name': hit_info['gene_name'],
                'hmm_source': hit_info['hmm_source'],
                'e_value': hit_info['e_value'],
            }
            window_counter += 1

        self.run.info('RT windows created', len(rt_windows))

        # Store for later use
        self.rt_windows = rt_windows

        return rt_windows



    def get_rt_window_query_records(self, contig_sequences, bin_contigs=None):
        """
        Extract RT window sequences and create query records for BLAST.

        Parameters
        ==========
        contig_sequences : dict
            Contig sequences as {contig_name: {'sequence': seq}}.
        bin_contigs : list or None
            If provided, only include RT windows from these contigs (for collections mode).

        Returns
        =======
        dict
            Query records as {section_id: sequence} for BLAST input.
        """
        # Get RT windows if not already computed
        if not hasattr(self, 'rt_windows') or not self.rt_windows:
            self.get_rt_windows()

        if not self.rt_windows:
            return {}

        # Extract window sequences and create query records
        # Format section_id to be compatible with existing parsing code:
        # {contig_name}_section_{window_id}_start_bp{start}_end_bp{end}
        query_records = {}
        for window_id, window_info in self.rt_windows.items():
            contig = window_info['contig']
            window_start = window_info['window_start']
            window_end = window_info['window_end']

            # Skip if not in target contigs (for collections mode)
            if bin_contigs is not None and contig not in bin_contigs:
                continue

            if contig not in contig_sequences:
                continue

            # Extract window sequence (window_end is inclusive, so add 1 for Python slicing)
            contig_seq = contig_sequences[contig]['sequence']
            window_seq = contig_seq[window_start:window_end + 1]

            # Create section_id compatible with existing parsing
            section_id = f"{contig}_section_{window_id}_start_bp{window_start}_end_bp{window_end}"
            query_records[section_id] = window_seq

        return query_records


    def run_blast_homology_mode(self, contig_sequences=None, bin_name=None, bin_contigs=None):
        """
        Run BLASTn using RT neighborhood windows as queries against contig sequences.

        This method is used for homology-based DGR detection. It:
        1. Gets RT windows from get_rt_windows()
        2. Extracts genomic sequences for each window
        3. BLASTs these TR-candidate regions against the genome
        4. Returns the path to the BLAST output for parsing

        Parameters
        ==========
        contig_sequences : dict or None
            Contig sequences as {contig_name: {'sequence': seq}}.
            If None, loads all contigs from contigs.db.
        bin_name : str or None
            Bin name for collections mode output naming.
        bin_contigs : list or None
            If provided, only include RT windows from these contigs.

        Returns
        =======
        str
            Path to the BLASTn output XML file.

        Raises
        ======
        ConfigError
            If no RT windows are found or no contig sequences are available.
        """
        self.run.info_single("Running BLAST for homology-based DGR detection...", nl_before=1)

        # Load contig sequences if not provided
        if contig_sequences is None:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
            contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
            contigs_db.disconnect()

        if not contig_sequences:
            raise ConfigError("No contig sequences found in the contigs database.")

        # Get query records (RT window sequences)
        query_records = self.get_rt_window_query_records(contig_sequences, bin_contigs)

        if not query_records:
            if bin_name:
                self.run.warning(f"No RT windows found in bin '{bin_name}'. Skipping homology-based detection for this bin.")
                return None
            else:
                raise ConfigError("No RT windows were found. Homology-based detection requires RT HMM hits "
                                "in the contigs database. Please run `anvi-run-hmms -I Reverse_Transcriptase`.")

        # Use unified run_blast function
        blast_output_path = self.run_blast(query_records, contig_sequences,
                                           mode='homology', bin_name=bin_name)

        # Store for later use
        self.homology_blast_output = blast_output_path

        return blast_output_path



    def process_homology_mode(self):
        """
        Run the complete homology-based DGR detection pipeline.

        This method orchestrates homology-based detection:
        1. Runs BLAST with RT windows as queries
        2. Parses BLAST results with vr_in_query=False (query=TR, subject=VR)
        3. Skips SNV-based filtering (apply_snv_filters=False)

        Returns
        =======
        dict
            The mismatch_hits dictionary containing homology-based DGR candidates.
            Each hit has detection_method='homology' and confidence='homology-based'.
        """
        self.run.info_single("Starting homology-based DGR detection pipeline...", nl_before=1, nl_after=1)

        # Step 1: Run BLAST with RT windows as queries
        blast_output = self.run_blast_homology_mode()

        # Check if BLAST was run and produced output
        if blast_output is None or os.stat(blast_output).st_size == 0:
            self.run.warning("No BLAST hits were found in homology mode. This could mean there are no "
                           "TR/VR pairs near the RT genes, or the RT windows don't contain template regions.",
                           header="NO HOMOLOGY-MODE HITS")
            self.mismatch_hits = defaultdict(lambda: defaultdict(dict))
            return self.mismatch_hits

        # Step 2: Parse BLAST results
        # In homology mode: query=TR (RT window), subject=VR
        # So we set vr_in_query=False
        self.mismatch_hits = defaultdict(lambda: defaultdict(dict))
        self.parse_and_process_blast_results(
            blast_output,
            bin_name=None,
            max_percent_identity=100,
            vr_in_query=False,  # Query contains TR (RT windows), subject contains VR
            apply_snv_filters=False  # No SNV data in homology mode
        )

        num_hits = sum(len(hits) for hits in self.mismatch_hits.values())
        self.run.info('Homology mode candidates found', num_hits)

        return self.mismatch_hits



    def merge_detection_results(self, activity_hits, homology_hits):
        """
        Merge results from activity-based and homology-based detection modes.

        When the same VR/TR pair is found by both methods, they are merged into a single
        entry with detection_method='both'. The activity-based confidence score is kept
        since it's more informative (based on SNV analysis).

        Parameters
        ==========
        activity_hits : dict
            Mismatch hits from activity-based detection (query=VR, subject=TR).
            Each hit has detection_method='activity'.
        homology_hits : dict
            Mismatch hits from homology-based detection (query=TR, subject=VR).
            Each hit has detection_method='homology'.

        Returns
        =======
        dict
            Merged mismatch_hits dictionary with detection_method indicating
            'activity', 'homology', or 'both' for each hit.
        """
        self.run.info_single("Merging results from activity and homology detection modes...", nl_before=1)

        merged_hits = defaultdict(lambda: defaultdict(dict))

        # Track statistics
        activity_only = 0
        homology_only = 0
        found_by_both = 0

        # First, add all activity hits to merged results
        for section_id, hits_dict in activity_hits.items():
            for hit_id, hit_data in hits_dict.items():
                merged_hits[section_id][hit_id] = hit_data.copy()
                merged_hits[section_id][hit_id]['detection_method'] = 'activity'

        # Now process homology hits, checking for overlaps
        for section_id, hits_dict in homology_hits.items():
            for hit_id, hit_data in hits_dict.items():
                # Extract key coordinates for matching
                # In homology mode, hit_data already has semantic VR/TR naming from parsing
                h_vr_contig = hit_data.get('query_contig')  # VR contig
                h_vr_start = hit_data.get('query_genome_start_position')  # VR start
                h_vr_end = hit_data.get('query_genome_end_position')  # VR end
                h_tr_contig = hit_data.get('subject_contig')  # TR contig
                h_tr_start = hit_data.get('subject_genome_start_position')  # TR start
                h_tr_end = hit_data.get('subject_genome_end_position')  # TR end

                # Look for matching entry in activity hits
                match_found = False
                for a_section_id, a_hits_dict in activity_hits.items():
                    for a_hit_id, a_hit_data in a_hits_dict.items():
                        a_vr_contig = a_hit_data.get('query_contig')
                        a_vr_start = a_hit_data.get('query_genome_start_position')
                        a_vr_end = a_hit_data.get('query_genome_end_position')
                        a_tr_contig = a_hit_data.get('subject_contig')
                        a_tr_start = a_hit_data.get('subject_genome_start_position')
                        a_tr_end = a_hit_data.get('subject_genome_end_position')

                        # Check if VR and TR coordinates overlap significantly
                        vr_overlap = (h_vr_contig == a_vr_contig and
                                     self.range_overlapping(h_vr_start, h_vr_end, a_vr_start, a_vr_end))
                        tr_overlap = (h_tr_contig == a_tr_contig and
                                     self.range_overlapping(h_tr_start, h_tr_end, a_tr_start, a_tr_end))

                        if vr_overlap and tr_overlap:
                            # Found a match - update detection_method to 'both'
                            # Keep the activity hit data (has SNV-based confidence)
                            merged_hits[a_section_id][a_hit_id]['detection_method'] = 'both'
                            match_found = True
                            found_by_both += 1
                            break
                    if match_found:
                        break

                if not match_found:
                    # No match found - add as homology-only hit
                    new_hit_id = f"{hit_id}_homology"
                    merged_hits[section_id][new_hit_id] = hit_data.copy()
                    merged_hits[section_id][new_hit_id]['detection_method'] = 'homology'
                    homology_only += 1

        # Count activity-only hits (those not upgraded to 'both')
        for section_id, hits_dict in merged_hits.items():
            for hit_id, hit_data in hits_dict.items():
                if hit_data.get('detection_method') == 'activity':
                    activity_only += 1

        # Report merge statistics
        total_merged = activity_only + homology_only + found_by_both
        self.run.info('Total merged candidates', total_merged)
        self.run.info('  Activity-only', activity_only)
        self.run.info('  Homology-only', homology_only)
        self.run.info('  Found by both', found_by_both)

        return merged_hits



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
        : tsv
            A tsv tabular file containing the template and variable regions
        """
        dgrs_dict = self.DGRs_found_dict

        if len(dgrs_dict) == 0:
            return

        self.run.info_single("Computing the closest HMMs to the Template Regions and printing them in your output tsv.", nl_before=1)

        dgrs_without_rt = set()

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
                if gene_callers_id not in found_HMMS_dict.keys():
                    found_HMMS_dict[gene_callers_id] = {'gene_name':gene_name,
                                                            'e_value':e_value,
                                                            'HMM_source':HMM_source}
                elif e_value < found_HMMS_dict[gene_callers_id]['e_value']:
                    found_HMMS_dict[gene_callers_id]['gene_name'] = gene_name
                    found_HMMS_dict[gene_callers_id]['e_value'] = e_value
                    found_HMMS_dict[gene_callers_id]['HMM_source'] = HMM_source

        # check if the gene_caller_id exists in genes_in_contigs
        for gene_callers_id, hmm_dict in found_HMMS_dict.items():
            # retrieve the 'start' and 'stop' values
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

        # look at general consensus TR in the level up so all the TRs have the same HMM if in the same DGR.
        for DGR_id, DGR_info in dgrs_dict.items():
            TR_start_position = DGR_info['TR_start_position']
            TR_end_position = DGR_info['TR_end_position']
            TR_middle_position = (TR_start_position + TR_end_position) / 2

            # initialize closest_distances dictionary inside the loop
            closest_distances = {}  # initialize an empty dictionary to store closest distances
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
                dgrs_without_rt.add(DGR_id)
                if anvio.DEBUG:
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

        # summary of DGRs without RT
        if dgrs_without_rt:
            self.run.warning(f"{PL('DGR', len(dgrs_without_rt))} had no Reverse Transcriptase HMM on the same contig "
                            f"as the Template Region: {', '.join(sorted(dgrs_without_rt))}. This does not affect the "
                            "outcome of the tool, but you may want to manually inspect your data.",
                            header="NO REVERSE TRANSCRIPTASE HMMs")

        return



    def create_found_tr_vr_tsv(self):
        """
        This function creates a tsv tabular format of the template and variable regions that are found from this tool.

        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======
        : tsv
            A tsv tabular file containing the template and variable regions

        """
        output_directory_path = self.output_directory
        output_prefix = os.path.basename(self.output_directory)

        dgrs_dict = self.DGRs_found_dict
        output_path_dgrs = os.path.join(output_directory_path, f"{output_prefix}_DGRs_found.tsv")
        headers = [
            "DGR", "VR", "VR_contig", "VR_frame_reported", "VR_sequence", "Midline",
            "VR_start_position", "VR_end_position", "VR_bin", "Mismatch %",
            "TR_contig", "TR_frame_Reported", "TR_sequence", "Base", "Reverse Complemented_from_BLAST",
            "TR_start_position", "TR_end_position", "TR_bin", "TR_in_gene", "HMM_source",
            "distance_to_HMM", "HMM_gene_name", "HMM_direction", "HMM_start",
            "HMM_stop", "HMM_gene_callers_id", "numb_of_mismatches", "numb_of_SNVs",
            "VR_TR_mismatch_positions", "snv_VR_positions", "best_amongst_multiple_TRs_for_one_VR",
            # New codon/SNV analysis columns
            "mismatch_codon_1", "mismatch_codon_2", "mismatch_codon_3", "pct_mismatch_codon_3",
            "vr_gene_id", "n_snvs_total", "n_snvs_explained", "n_snvs_unexplained", "n_explained_diverse", "pct_snvs_explained",
            "snv_codon_1", "snv_codon_2", "snv_codon_3", "pct_snv_codon_3",
            "snv_supporting_sample", "confidence", "confidence_reasons",
            "detection_method"]

        # check if either dictionary is empty or lacks meaningful keys
        if not any(dgrs_dict.values()):
            self.run.warning("No DGRS were found so no output file will be written :( However, you can go "
                            "back and tinker with the parameters of this tool if you believe this should not "
                            "be the case. Anvi'o wishes you a nice day :)", header="NO DIVERSITY-GENERATING RETROELEMENTS")
            return

        # open the TSV file and write headers and rows
        with open(output_path_dgrs, 'w', encoding='utf-8', newline='\n') as f:
            # Write header line
            f.write('\t'.join(headers) + '\n')

            for dgr, tr in dgrs_dict.items():
                for vr, vr_data in tr['VRs'].items():
                    # populate tsv_row for DGRs_found.tsv format
                    csv_row = [
                        dgr, vr, vr_data['VR_contig'], vr_data.get('VR_frame', 'N/A'),
                        vr_data['VR_sequence'], vr_data['midline'], vr_data['VR_start_position'],
                        vr_data['VR_end_position'], vr_data.get('VR_bin'),
                        vr_data['percentage_of_mismatches'], tr['TR_contig'],
                        vr_data.get('TR_frame', 'N/A'), vr_data['TR_sequence'],
                        tr['base'], vr_data['TR_reverse_complement'],
                        vr_data['TR_start_position'], vr_data['TR_end_position'],
                        tr.get('TR_bin', 'N/A'), tr.get('TR_in_gene', 'N/A'), tr['HMM_source'], tr["distance_to_HMM"],
                        tr["HMM_gene_name"], tr["HMM_direction"], tr["HMM_start"],
                        tr["HMM_stop"], tr["HMM_gene_callers_id"],
                        vr_data["numb_of_mismatches"], vr_data["numb_of_SNVs"], vr_data["VR_TR_mismatch_positions"],
                        vr_data["snv_VR_positions"], vr_data["best_amongst_multiple_TRs_for_one_VR"],
                        # New codon/SNV analysis columns
                        vr_data.get("mismatch_codon_1", 0), vr_data.get("mismatch_codon_2", 0),
                        vr_data.get("mismatch_codon_3", 0), vr_data.get("pct_mismatch_codon_3", 0),
                        vr_data.get("vr_gene_id", "N/A"),
                        vr_data.get("n_snvs_total", 0), vr_data.get("n_snvs_explained", 0),
                        vr_data.get("n_snvs_unexplained", 0), vr_data.get("n_explained_diverse", 0), vr_data.get("pct_snvs_explained", 100),
                        vr_data.get("snv_codon_1", 0), vr_data.get("snv_codon_2", 0),
                        vr_data.get("snv_codon_3", 0), vr_data.get("pct_snv_codon_3", 0),
                        vr_data.get("snv_supporting_sample", "N/A"),
                        vr_data.get("confidence", "N/A"), vr_data.get("confidence_reasons", []),
                        vr_data.get("detection_method", "activity")
                    ]
                    f.write('\t'.join(map(str, csv_row)) + '\n')
        return



    def recover_genomic_context_surrounding_dgrs(self):
        """
        Learn about what surrounds the variable region sites of each found DGR

        Returns
        =======
        None
            Updates the instance attribute `self.genomic_context_surrounding_dgrs` with a
            nested dictionary structure
        """

        # in which we will store the genomic context that surrounds dgrs for downstream fun
        self.genomic_context_surrounding_dgrs = {}

        dgrs_dict = self.DGRs_found_dict

        if len(dgrs_dict) == 0:
            return

        # we know when we are not wanted
        if self.skip_recovering_genomic_context:
            self.run.info_single('Skipping genomic context recovery due to self.skip_recovering_genomic_context being True', nl_before=1)
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

            # initialize TR context
            tr_context_genes = {}
            self.genomic_context_surrounding_dgrs[dgr_id] = {}

            is_tr_in_gene = False

            if TR_contig_name not in gene_calls_per_TR_contig:
                where_clause = f'''contig="{TR_contig_name}" and source="{self.gene_caller_to_consider_in_context}"'''
                gene_calls_per_TR_contig[TR_contig_name] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

            gene_calls_in_TR_contig = gene_calls_per_TR_contig[TR_contig_name]

            if not gene_calls_in_TR_contig:
                # case: TR has no genes
                trs_with_no_gene_calls_around.add(dgr_id)
                dgrs_dict[dgr_id]['TR_in_gene'] = False
                self.genomic_context_surrounding_dgrs[dgr_id]["TR"] = {}  # empty TR context
            else:
                # case: TR has genes → process them
                min_distance_to_TR_start, min_distance_to_TR_end = float('inf'), float('inf')
                closest_gene_call_to_TR_start, closest_gene_call_to_TR_end = None, None

                # process TR gene calls
                for gene_callers_id, gene_call in gene_calls_in_TR_contig.items():
                    # cache distance calculations to avoid redundant abs() calls
                    dist_to_start = abs(gene_call['start'] - TR_start)
                    dist_to_end = abs(gene_call['start'] - TR_end)

                    if dist_to_start < min_distance_to_TR_start:
                        closest_gene_call_to_TR_start = gene_callers_id
                        min_distance_to_TR_start = dist_to_start

                    if dist_to_end < min_distance_to_TR_end:
                        closest_gene_call_to_TR_end = gene_callers_id
                        min_distance_to_TR_end = dist_to_end

                TR_range = range(closest_gene_call_to_TR_start - self.num_genes_to_consider_in_context,
                                closest_gene_call_to_TR_end + self.num_genes_to_consider_in_context)
                tr_gene_caller_ids_of_interest = [c for c in TR_range if c in gene_calls_in_TR_contig]

                # === BATCH QUERIES BEFORE THE LOOP (query once, use many times) ===
                # 1. Fetch all functions for genes of interest in one query
                tr_functions_by_gene = defaultdict(list)
                if function_sources_found and tr_gene_caller_ids_of_interest:
                    where_clause = '''gene_callers_id IN (%s)''' % (', '.join([str(g) for g in tr_gene_caller_ids_of_interest]))
                    all_tr_function_hits = list(contigs_db.db.get_some_rows_from_table_as_dict(
                        t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
                    for hit in all_tr_function_hits:
                        tr_functions_by_gene[hit['gene_callers_id']].append(hit)

                # 2. Fetch contig sequence once
                where_clause = f'''contig == "{TR_contig_name}"'''
                tr_contig_seq_data = contigs_db.db.get_some_rows_from_table_as_dict(
                    t.contig_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                tr_contig_sequence = tr_contig_seq_data[TR_contig_name]['sequence'] if TR_contig_name in tr_contig_seq_data else ''

                # 3. Fetch all AA sequences in one query
                tr_aa_sequences = {}
                if tr_gene_caller_ids_of_interest:
                    where_clause = '''gene_callers_id IN (%s)''' % (', '.join([str(g) for g in tr_gene_caller_ids_of_interest]))
                    tr_aa_sequences = contigs_db.db.get_some_rows_from_table_as_dict(
                        t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)

                # === NOW LOOP - only in-memory operations ===
                for gene_callers_id in tr_gene_caller_ids_of_interest:
                    gene_call = gene_calls_in_TR_contig[gene_callers_id]
                    gene_call['gene_callers_id'] = gene_callers_id

                    # O(1) dict lookup instead of database query
                    gene_call['functions'] = tr_functions_by_gene.get(gene_callers_id, [])

                    # Use pre-fetched contig sequence
                    dna_sequence = tr_contig_sequence[gene_call['start']:gene_call['stop']]

                    rev_compd = None
                    if gene_call['direction'] == 'f':
                        gene_call['DNA_sequence'] = dna_sequence
                        rev_compd = False
                    else:
                        gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)
                        rev_compd = True

                    # O(1) dict lookup instead of database query
                    gene_call['AA_sequence'] = tr_aa_sequences.get(gene_callers_id, {}).get('sequence', '')

                    gene_call['length'] = gene_call['stop'] - gene_call['start']

                    header = '|'.join([f"contig:{gene_call['contig']}",
                                    f"start:{gene_call['start']}",
                                    f"stop:{gene_call['stop']}",
                                    f"direction:{gene_call['direction']}",
                                    f"rev_compd:{rev_compd}",
                                    f"length:{gene_call['length']}"])
                    gene_call['header'] = ' '.join([str(gene_callers_id), header])

                    # store the gene call in the dictionary using gene_callers_id as the key
                    tr_context_genes[gene_callers_id] = gene_call

                    # check if the TR is inside this gene
                    if TR_start >= gene_call['start'] and TR_end <= gene_call['stop']:
                        is_tr_in_gene = True

            # after processing all the gene calls for the TR, add the result to dgrs_dict
            dgrs_dict[dgr_id]['TR_in_gene'] = is_tr_in_gene

            self.genomic_context_surrounding_dgrs[dgr_id] = copy.deepcopy(tr_context_genes)

            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key

                self.progress.update(f"{dgr_id} / {vr_id}")

                VR_contig = vr_data.get('VR_contig')
                VR_start = vr_data.get('VR_start_position')
                VR_end = vr_data.get('VR_end_position')

                # initialize VR context
                vr_context_genes = {}

                if VR_contig not in gene_calls_per_VR_contig:
                    where_clause = f'''contig="{VR_contig}" and source="{self.gene_caller_to_consider_in_context}"'''
                    gene_calls_per_VR_contig[VR_contig] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

                gene_calls_in_VR_contig = gene_calls_per_VR_contig[VR_contig]

                if not len(gene_calls_in_VR_contig):
                    vrs_with_no_gene_calls_around.add(vr_id)
                    continue

                min_distance_to_VR_start, min_distance_to_VR_end = float('inf'), float('inf')
                closest_gene_call_to_VR_start, closest_gene_call_to_VR_end = None, None
                for gene_callers_id, gene_call in gene_calls_in_VR_contig.items():
                    # cache distance calculations to avoid redundant abs() calls
                    dist_to_start = abs(gene_call['start'] - VR_start)
                    dist_to_end = abs(gene_call['start'] - VR_end)

                    if dist_to_start < min_distance_to_VR_start:
                        closest_gene_call_to_VR_start = gene_callers_id
                        min_distance_to_VR_start = dist_to_start

                    if dist_to_end < min_distance_to_VR_end:
                        closest_gene_call_to_VR_end = gene_callers_id
                        min_distance_to_VR_end = dist_to_end

                VR_range = range(closest_gene_call_to_VR_start - self.num_genes_to_consider_in_context,
                                closest_gene_call_to_VR_end + self.num_genes_to_consider_in_context)
                vr_gene_caller_ids_of_interest = [c for c in VR_range if c in gene_calls_in_VR_contig]

                # === BATCH QUERIES BEFORE THE LOOP (query once, use many times) ===
                # 1. Fetch all functions for genes of interest in one query
                vr_functions_by_gene = defaultdict(list)
                if function_sources_found and vr_gene_caller_ids_of_interest:
                    where_clause = '''gene_callers_id IN (%s)''' % (', '.join([str(g) for g in vr_gene_caller_ids_of_interest]))
                    all_vr_function_hits = list(contigs_db.db.get_some_rows_from_table_as_dict(
                        t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
                    for hit in all_vr_function_hits:
                        vr_functions_by_gene[hit['gene_callers_id']].append(hit)

                # 2. Fetch contig sequence once
                where_clause = f'''contig == "{VR_contig}"'''
                vr_contig_seq_data = contigs_db.db.get_some_rows_from_table_as_dict(
                    t.contig_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
                vr_contig_sequence = vr_contig_seq_data[VR_contig]['sequence'] if VR_contig in vr_contig_seq_data else ''

                # 3. Fetch all AA sequences in one query
                vr_aa_sequences = {}
                if vr_gene_caller_ids_of_interest:
                    where_clause = '''gene_callers_id IN (%s)''' % (', '.join([str(g) for g in vr_gene_caller_ids_of_interest]))
                    vr_aa_sequences = contigs_db.db.get_some_rows_from_table_as_dict(
                        t.gene_amino_acid_sequences_table_name, where_clause=where_clause, error_if_no_data=False)

                # === NOW LOOP - only in-memory operations ===
                for gene_callers_id in vr_gene_caller_ids_of_interest:
                    gene_call = gene_calls_in_VR_contig[gene_callers_id]
                    gene_call['gene_callers_id'] = gene_callers_id

                    # O(1) dict lookup instead of database query
                    gene_call['functions'] = vr_functions_by_gene.get(gene_callers_id, [])

                    # Use pre-fetched contig sequence
                    dna_sequence = vr_contig_sequence[gene_call['start']:gene_call['stop']]

                    rev_compd = None
                    if gene_call['direction'] == 'f':
                        gene_call['DNA_sequence'] = dna_sequence
                        rev_compd = False
                    else:
                        gene_call['DNA_sequence'] = utils.rev_comp(dna_sequence)
                        rev_compd = True

                    # O(1) dict lookup instead of database query
                    gene_call['AA_sequence'] = vr_aa_sequences.get(gene_callers_id, {}).get('sequence', '')

                    gene_call['length'] = gene_call['stop'] - gene_call['start']

                    header = '|'.join([f"contig:{gene_call['contig']}",
                                    f"start:{gene_call['start']}",
                                    f"stop:{gene_call['stop']}",
                                    f"direction:{gene_call['direction']}",
                                    f"rev_compd:{rev_compd}",
                                    f"length:{gene_call['length']}"])
                    gene_call['header'] = ' '.join([str(gene_callers_id), header])

                    # store the gene call in the dictionary using gene_callers_id as the key
                    vr_context_genes[gene_callers_id] = gene_call

                self.genomic_context_surrounding_dgrs[dgr_id][vr_id] = copy.deepcopy(vr_context_genes)

        contigs_db.disconnect()
        self.progress.end()
        self.run.info_single('Completed recovering genomic context surrounding the DGRs', nl_before=1)

        if len(self.genomic_context_surrounding_dgrs) > 0:
            # count total VRs with context recovered across all DGRs
            total_vrs_with_context = sum(
                len([k for k in dgr_context.keys() if isinstance(k, str)])
                for dgr_context in self.genomic_context_surrounding_dgrs.values()
            )
            self.run.info(f"[Genomic Context] Searched {PL('DGR', len(dgrs_dict))}",
                        f"recovered TR context for {PL('DGR', len(self.genomic_context_surrounding_dgrs))}",
                        f"and VR context for {PL('VR', total_vrs_with_context)}",
                        nl_before=1,
                        lc="yellow")

        if len(trs_with_no_gene_calls_around):
            self.run.warning(
                'No gene calls around the following TRs. '
                f"Here is the list in case you would like to track them down: {', '.join(trs_with_no_gene_calls_around)}."
            )

        if len(vrs_with_no_gene_calls_around):
            self.run.warning(
                'No gene calls around the following VRs. '
                f"Here is the list in case you would like to track them down: {', '.join(vrs_with_no_gene_calls_around)}."
            )

        if not len(self.genomic_context_surrounding_dgrs):
            self.run.warning(f"Even though the tool went through all {PL('DGR', len(dgrs_dict))} "
                            "it was unable to recover any genomic context for any of them. So your final reports will "
                            "not include any insights into the surrounding genomic context of all your DGRs ( a little bit sad "
                            "but otherwise you will be fine).")



    def report_genomic_context_surrounding_dgrs(self):
        """
        Reports two long-format output files for genes and functions around each DGR element
        STOLEN (modified) FROM INVERSIONS CODE (line 1925)
        Generate per-DGR reports of genes and annotated functions surrounding TR and VRs.

        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======
        None
            Writes report files to the output directory and logs their paths.
        """

        if self.skip_recovering_genomic_context:
            self.run.info_single('Skipping reporting genomic context due to self.skip_recovering_genomic_context being True', nl_before=1)
            return

        if not len(self.genomic_context_surrounding_dgrs):
            self.run.warning("No genomic context data available to report on any of the DGRs :(")
            return

        dgrs_dict = self.DGRs_found_dict

        # we are in business
        genes_output_headers = ["gene_callers_id", "start", "stop", "direction", "partial", "call_type", "source", "version", "contig"]
        functions_output_headers = ["gene_callers_id", "source", 'accession', 'function']

        # tracking for summary
        num_tr_reports = 0
        num_vr_reports = 0
        vrs_without_genes = set()
        unexpected_gene_calls = []

        # process each DGR and its VRs
        for dgr_key, dgr_data in dgrs_dict.items():
            # assuming dgr_key itself is the dgr_id or a dictionary containing it
            dgr_id = dgr_key  # If dgr_key is the dgr_id itself

            # create output directory for DGR
            dgr_directory = os.path.join(self.output_directory, "PER_DGR", dgr_id)
            filesnpaths.gen_output_directory(dgr_directory, delete_if_exists=False)

            # TR output paths
            tr_genes_output_path = os.path.join(dgr_directory, 'TR_SURROUNDING-GENES.txt')
            tr_functions_output_path = os.path.join(dgr_directory, 'TR_SURROUNDING-FUNCTIONS.txt')

            with open(tr_genes_output_path, 'w') as tr_genes_output, open(tr_functions_output_path, 'w') as tr_functions_output:
                tr_genes_output.write("dgr_id\tentry_type\t%s\n" % '\t'.join(genes_output_headers))
                tr_functions_output.write("dgr_id\t%s\n" % '\t'.join(functions_output_headers))

                # create fake gene call entries for TR and VRs:
                d = dict([(h, '') for h in genes_output_headers])

                # fill in non-empty data for the TR in the DGR and insert it:
                d['contig'] = dgr_data.get('TR_contig')  # Use dgr_data to access TR info
                d['start'] = dgr_data.get('TR_start_position')
                d['stop'] = dgr_data.get('TR_end_position')
                tr_genes_output.write(f"{dgr_id}_TR\tTEMPLATE_REGION\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                # check if there are surrounding genes for the TR and write them
                if dgr_id in self.genomic_context_surrounding_dgrs:
                    # get the context for this DGR
                    dgr_context = self.genomic_context_surrounding_dgrs[dgr_id]

                    # filter for TR genes (integer keys only, not VR string keys)
                    tr_genes = {k: v for k, v in dgr_context.items() if isinstance(k, int)}

                    # iterate over the gene_callers_id and the associated gene_call dictionary
                    for gene_callers_id, gene_call in tr_genes.items():
                        # ensure the gene_call is a dictionary
                        if isinstance(gene_call, dict):
                            # write gene information to the output, ensuring we access the correct fields
                            tr_genes_output.write(f"{dgr_id}_TR\tGENE\t%s\n" % '\t'.join([f"{gene_call.get(h, '')}" for h in genes_output_headers]))

                            # if 'functions' is present in the gene_call, write functions to the output
                            if 'functions' in gene_call and gene_call['functions']:
                                for hit in gene_call['functions']:
                                    tr_functions_output.write(f"{dgr_id}_TR\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                            else:
                                # write placeholder if no functions are found
                                tr_functions_output.write(f"{dgr_id}_TR\t{gene_call.get('gene_callers_id', '')}\t\t\t\n")
                        else:
                            # track unexpected gene_call types (should be dict)
                            unexpected_gene_calls.append((dgr_id, gene_callers_id, type(gene_call).__name__))
                            if len(unexpected_gene_calls) == 1:
                                self.run.warning(f"Unexpected data type encountered for gene_call (expected dict, got {type(gene_call).__name__}). "
                                                "This is highly unusual and may indicate data corruption. "
                                                "Will continue processing but please check your data.",
                                                header="UNEXPECTED DATA TYPE")

                # log information about the reporting files
                num_tr_reports += 1
                if anvio.DEBUG:
                    self.run.info(f"Reporting file on gene context for {dgr_id} TR", tr_genes_output_path)
                    self.run.info(f"Reporting file on functional context for {dgr_id} TR", tr_functions_output_path, nl_after=1)

            # fill in non-empty data for each VR in the DGR and insert it:
            for vr_key, vr_data in dgr_data['VRs'].items():
                vr_id = vr_key
                vr_directory = os.path.join(dgr_directory, f"VR_{vr_id}")
                filesnpaths.gen_output_directory(vr_directory, delete_if_exists=False)

                vr_genes_output_path = os.path.join(vr_directory, 'VR_SURROUNDING-GENES.txt')
                vr_functions_output_path = os.path.join(vr_directory, 'VR_SURROUNDING-FUNCTIONS.txt')

                with open(vr_genes_output_path, 'w') as vr_genes_output, open(vr_functions_output_path, 'w') as vr_functions_output:
                    vr_genes_output.write("dgr_id\tentry_type\t%s\n" % '\t'.join(genes_output_headers))
                    vr_functions_output.write("dgr_id\t%s\n" % '\t'.join(functions_output_headers))

                    # create fake gene call entry for VR:
                    d['contig'] = vr_data.get('VR_contig')
                    d['start'] = vr_data.get('VR_start_position')
                    d['stop'] = vr_data.get('VR_end_position')
                    vr_genes_output.write(f"{dgr_id} {vr_id}\tVARIABLE_REGION\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                    # check if there are surrounding genes for the VR and write them
                    # VR genes are stored at: self.genomic_context_surrounding_dgrs[dgr_id][vr_id]
                    if dgr_id in self.genomic_context_surrounding_dgrs and vr_id in self.genomic_context_surrounding_dgrs[dgr_id]:
                        vr_genes = self.genomic_context_surrounding_dgrs[dgr_id][vr_id]

                        # Iterate over VR genes
                        for gene_callers_id, gene_call in vr_genes.items():
                            if isinstance(gene_call, dict):
                                vr_genes_output.write(f"{dgr_id}_{vr_id}\tGENE\t%s\n" % '\t'.join([f"{gene_call.get(h, '')}" for h in genes_output_headers]))

                            if 'functions' in gene_call:
                                for hit in gene_call['functions']:
                                    vr_functions_output.write(f"{dgr_id} {vr_id}\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                            else:
                                vr_functions_output.write(f"{dgr_id} {vr_id}\t{gene_call['gene_callers_id']}\t\t\t\n")
                    else:
                        vrs_without_genes.add(f"{dgr_id}_{vr_id}")
                        if anvio.DEBUG:
                            self.run.info_single(f"No VR genes found for {dgr_id} {vr_id}", nl_before=1)

                    num_vr_reports += 1
                    if anvio.DEBUG:
                        self.run.info(f'    Reporting file on gene context for {dgr_id} {vr_id}', vr_genes_output_path)
                        self.run.info(f'    Reporting file on functional context for {dgr_id} {vr_id}', vr_functions_output_path, nl_after=1)

        # summary of genomic context reporting
        self.run.info_single(f"Genomic context reports generated for {PL('TR', num_tr_reports)} and {PL('VR', num_vr_reports)} "
                             f"in '{self.output_directory}/PER_DGR/'", nl_before=1)
        if vrs_without_genes:
            self.run.warning(f"{PL('VR', len(vrs_without_genes))} had no surrounding genes: {', '.join(sorted(vrs_without_genes))}",
                            header="VRs WITHOUT SURROUNDING GENES")
        if len(unexpected_gene_calls) > 1:
            self.run.warning(f"Encountered {len(unexpected_gene_calls)} gene calls with unexpected data types in total.",
                            header="UNEXPECTED DATA TYPE SUMMARY")


    # function to get the consensus base
    @staticmethod
    def get_consensus_base(row):
        """
        Determine the consensus base from a row of nucleotide frequencies.
        """
        for nucleotide in nucleotides:
            if row[nucleotide] > 0.5:  # assuming a threshold for consensus, adjust as necessary
                return nucleotide
        return None



    @staticmethod
    def compute_dgr_variability_profiling_per_vr(input_queue, output_queue, samples_artifact, primers_dict, output_directory_path, run=run_quiet, progress=progress_quiet, sample_primers_dict=None, use_sample_primers=False):
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
            if sample_name is None:
                break

            # extract sample-specific primers from the nested structure
            # convert from {primer_name: {sample_name: {data}}}
            # to {primer_name: {data}} for this specific sample
            primers_for_sample = {}
            for primer_name, samples_data in primers_dict.items():
                if sample_name in samples_data:
                    primers_for_sample[primer_name] = samples_data[sample_name]

            samples_dict = samples_artifact.as_dict()
            samples_dict_for_sample = {sample_name: samples_dict[sample_name]}
            sample_artifact = SamplesTxt.from_dict(samples_dict_for_sample)

            # setup the args object
            args = argparse.Namespace(samples_artifact= sample_artifact,
                                    primers_dict= primers_for_sample,
                                    output_dir=output_directory_path,
                                    only_report_primer_matches = True
                                    )

            s = PrimerSearch(args, run=run, progress=progress)
            sample_dict, primer_hits = s.process(return_dicts = True)

            output_queue.put(sample_name)



    def init_vr_contigs(self):
        """
        """
        vr_contigs = ()
        for dgr_id, dgr_data in self.DGRs_found_dict.items():
            for vr_id, vr_data in dgr_data['VRs'].items():
                vr_contig = vr_data['VR_contig']
                vr_contigs += (vr_contig,)

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        where_clause = '''contig IN (%s)''' % (', '.join([f'"{str(g)}"' for g in vr_contigs]))
        self.vr_contig_sequences = contigs_db.db.get_some_rows_from_table_as_dict(t.contig_sequences_table_name, where_clause=where_clause, error_if_no_data=False)
        contigs_db.disconnect()
        return



    def generate_primers_for_vrs(self, dgrs_dict):
        """
        A function to generate primers for each and every VR. These are composed of not only an initial primer sequence before the VR
        but also of anchor points in the VR. These anchor points are the bases in the TR that are not A bases and only at the places the TR and VR match.
        The primer_sequence is then composed of the initial primer sequence, the masked primer sequence, and a flag stating if the primer is a L_to_R or not (left to right or right to left) because they are all on the same strand. These primers don't take the total primer length into account yet.

                    TR:GCTAACTGACATAATT
        masked_primer :GCT..C.G.C.T..TT
                    VR:GCTCACGGACTTCATT

        Parameters
        ==========
        DGRs_found_dict : dict
            A dictionary containing the template and variable regions

        Returns
        =======
        primers_dict : dict
            A dictionary containing the various primers for each sample including the compositional sections of each primer (i.e. masked primer and the initial primer)
        """

        # create primers dictionary
        primers_dict = {}

        # initialise the vr contig sequences
        self.init_vr_contigs()
        sample_names = set(self.samples_artifact.samples())

        if self.pre_computed_dgrs_path:
            self.init_snv_table()

        # Sanity check for samples (move this outside the main loops)
        sample_names_in_snv_table = set(self.snv_panda['sample_id'])
        samples_missing_in_snv_table = sample_names.difference(sample_names_in_snv_table)

        if anvio.DEBUG:
            self.run.info("Samples given", ", ".join(list(sample_names)))
            self.run.info("Samples in profile.db's nucleotide variability table", ", ".join(list(sample_names_in_snv_table)))
            self.run.info("Missing samples from profile.db's nucleotide variability table", ", ".join(list(samples_missing_in_snv_table)))

        if sample_names == samples_missing_in_snv_table:
            raise ConfigError(f"Anvi'o is not angry, just disappointed :/ You gave 'anvi-report-dgrs' these samples ({list(sample_names)}), but you have none of them in your profile.db.")

        for dgr_id, dgr_data in dgrs_dict.items():
            for vr_id, vr_data in dgr_data['VRs'].items():
                vr_contig = vr_data['VR_contig']
                contig_sequence = self.vr_contig_sequences[vr_contig]['sequence']
                contig_length = len(contig_sequence)

                # get TR VR sequences
                VR_sequence = vr_data['VR_sequence']
                TR_sequence = vr_data['TR_sequence']
                vr_start = vr_data.get('VR_start_position')
                vr_end = vr_data.get('VR_end_position')

                # keep the frame info
                VR_frame = vr_data['VR_frame']

                # decide where initial primer comes from
                skip_initial_primer = False
                initial_primer_right = False

                # sanity check that the VR is short enough to have an initial primer
                if (len(VR_sequence) + self.initial_primer_length) > contig_length:
                    skip_initial_primer = True
                    self.run.warning(
                        f"The VR is too close to the start and the end of the contig {vr_contig} "
                        f"for DGR {dgr_id} VR {vr_id}. The initial primer will therefore be non-existent.")

                # which way do we get the primer or skip an initial primer if VR at start/end
                if (vr_start < self.initial_primer_length and VR_frame == 1) or (VR_frame == -1 and not vr_end > (contig_length - self.initial_primer_length)):
                    initial_primer_right = True

                # store the flag in vr_data so we can reuse it later
                vr_data['initial_primer_right'] = initial_primer_right

                # get initial primer region if allowed
                vr_initial_primer_region = None
                if not skip_initial_primer:
                    if initial_primer_right:
                        vr_initial_primer_region = contig_sequence[vr_end : vr_end + self.initial_primer_length]
                    else:
                        vr_initial_primer_region = contig_sequence[vr_start - self.initial_primer_length : vr_start]

                # base primer based on TR VR mismatches and mutagenesis base in the TR sequence
                base_vr_masked_primer_list = []
                for i, (tr_base, vr_base) in enumerate(zip(TR_sequence, VR_sequence)):
                    if tr_base == dgr_data['base']:
                        base_vr_masked_primer_list.append('.')
                    elif tr_base != vr_base:
                        base_vr_masked_primer_list.append('.')
                    else:
                        base_vr_masked_primer_list.append(tr_base)

                base_vr_masked_primer = ''.join(base_vr_masked_primer_list)

                # initialize dict entry for this VR
                dgr_vr_key = f"{dgr_id}_{vr_id}_Primer"
                primers_dict[dgr_vr_key] = {}

                # track if we warned about missing initial primer (warn only once)
                warned_about_initial_primer = False

                # now create sample-specific primers if we have SNV info
                for sample_name in sample_names:
                    if anvio.DEBUG:
                        self.run.info_single(f"Processing sample {sample_name} for DGR {dgr_id} VR {vr_id}", nl_before=1)

                    # get SNVs for this sample in the VR region
                    sample_snvs = self.snv_panda[
                        (self.snv_panda['sample_id'] == sample_name) &
                        (self.snv_panda['contig_name'] == vr_contig) &
                        (self.snv_panda['pos_in_contig'] >= vr_start) &
                        (self.snv_panda['pos_in_contig'] < vr_end)
                    ]
                    snv_positions = set(sample_snvs['pos_in_contig'])

                    # apply SNV-based masking
                    sample_primer_list = list(base_vr_masked_primer)
                    for pos in snv_positions:
                        idx = pos - vr_start  # position relative to VR start
                        if 0 <= idx < len(sample_primer_list):
                            sample_primer_list[idx] = '.'

                    vr_masked_primer = ''.join(sample_primer_list)

                    # determine if this sample used the original (no-SNV) primer
                    used_original_primer = (len(snv_positions) == 0)

                    # reverse complement the masked primer if the VR is in reverse frame
                    if VR_frame == -1:
                        vr_masked_primer = utils.rev_comp(vr_masked_primer)

                    # combine initial + masked primer
                    if not skip_initial_primer:
                        if initial_primer_right:
                            primer_sequence = vr_masked_primer + vr_initial_primer_region
                        else:
                            primer_sequence = vr_initial_primer_region + vr_masked_primer
                    else:
                        primer_sequence = vr_masked_primer
                        if not warned_about_initial_primer:
                            self.run.warning(
                                f"Creating primer for DGR {dgr_id} VR {vr_id} without an initial primer region. "
                                f"Only a masked primer is present, more risky."
                            )
                            warned_about_initial_primer = True

                    # cut to correct total length
                    if initial_primer_right:
                        # cut from the LEFT (keep the right side)
                        primer_sequence = primer_sequence[-self.whole_primer_length:]
                    else:
                        # cut from the RIGHT (keep the left side)
                        primer_sequence = primer_sequence[:self.whole_primer_length]

                    # store sample-specific primer data
                    primers_dict[dgr_vr_key][sample_name] = {
                        'used_original_primer': used_original_primer,
                        'initial_primer_sequence': vr_initial_primer_region if not skip_initial_primer else 'N/A',
                        'vr_masked_primer': vr_masked_primer,
                        'primer_sequence': primer_sequence
                    }

        return primers_dict



    def populate_dgrs_dict_from_input_file(self):
        """
        Get the DGRs from a previously generated output file and populate the self.DGRs_found_dict with proper formatting.

        Parameters
        ==========
        DGRs_found_tsv : tsv
            A tsv generated from `anvi-report-dgrs` containing all associated information about the template and variable regions

        Returns
        =======
        dgrs_dict : nested dict
            A self.DGRs_found_dict that has all the information associated with the DGRs
        """

        dgrs_dict = utils.get_TAB_delimited_file_as_dictionary(
            self.pre_computed_dgrs_path,
            ignore_duplicated_keys=True
        )

        self.DGRs_found_dict = {}

        for dgr_id, row in dgrs_dict.items():

            # skip empty rows
            if all(v in [None, '', 'N/A'] for v in row.values()):
                continue

            # Ensure DGR slot exists
            if dgr_id not in self.DGRs_found_dict:
                self.DGRs_found_dict[dgr_id] = {'VRs': {}}

            # TR fields
            tr_fields = [
                'TR_sequence', 'Base', 'Reverse Complemented_from_BLAST',
                'TR_start_position', 'TR_end_position', 'TR_bin',
                'TR_contig', 'HMM_gene_callers_id', 'distance_to_HMM',
                'HMM_start', 'HMM_stop', 'HMM_direction', 'HMM_source',
                'HMM_gene_name'
            ]

            # Add TR fields that have valid values
            for field, converter in self.essential_keys_to_describe_dgrs:
                if field in tr_fields:
                    raw_value = row.get(field)
                    if raw_value not in [None, '', 'N/A']:   # skip empty fields
                        converted = self.safe_convert(
                            raw_value, converter,
                            field_name=field,
                            dgr_id=dgr_id
                        )
                        if converted is not None:
                            if field == 'Base':
                                target_field = 'base'
                            else:
                                target_field = field
                            self.DGRs_found_dict[dgr_id][target_field] = converted

            # VR identification
            vr_id = row.get('VR', 'VR_001')
            if vr_id not in self.DGRs_found_dict[dgr_id]['VRs']:
                self.DGRs_found_dict[dgr_id]['VRs'][vr_id] = {}

            # VR fields
            vr_fields = [
                'VR_sequence', 'Midline', 'Mismatch %', 'VR_start_position',
                'VR_end_position', 'VR_bin', 'VR_contig', 'VR_frame_reported',
                'TR_start_position', 'TR_end_position', 'TR_sequence',
                'TR_frame_Reported', 'VR_TR_mismatch_positions',
                'snv_VR_positions', 'numb_of_mismatches', 'numb_of_SNVs',
                'best_amongst_multiple_TRs_for_one_VR', 'mismatch_codon_1',
                'mismatch_codon_2', 'mismatch_codon_3', 'pct_mismatch_codon_3',
                'vr_gene_id', 'n_snvs_total', 'n_snvs_explained',
                'n_snvs_unexplained', 'n_explained_diverse', 'pct_snvs_explained', 'snv_codon_1',
                'snv_codon_2', 'snv_codon_3', 'pct_snv_codon_3', 'confidence',
                'confidence_reasons'
            ]

            # Add VR fields that have valid values
            for field, converter in self.essential_keys_to_describe_dgrs:
                if field in vr_fields:
                    raw = row.get(field)
                    if raw not in [None, '', 'N/A']:
                        converted = self.safe_convert(raw, converter, field_name=field, dgr_id=dgr_id, vr_id=vr_id)

                        if converted is not None:
                            # rename VR_frame_reported → VR_frame
                            if field == 'VR_frame_reported':
                                target_field = 'VR_frame'
                            else:
                                target_field = field

                            self.DGRs_found_dict[dgr_id]['VRs'][vr_id][target_field] = converted



    def safe_convert(self, value, converter, field_name='', dgr_id='', vr_id=''):
        """
        A very small helper function for empty fields in tab-delimited
        """
        if value in [None, '', 'N/A']:
            return None

        try:
            return converter(value)
        except Exception as e:
            vr_info = f", VR '{vr_id}'" if vr_id else ""
            raise ConfigError(f"Could not convert field '{field_name}' for DGR '{dgr_id}'{vr_info}: {e}. "
                             f"Your pre-computed DGR input file '{self.pre_computed_dgrs_path}' appears to have "
                             "malformed data. Please check the file format and try again.")



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

        if self.skip_compute_DGR_variability_profiling:
            return

        # define defaults
        dgrs_dict = self.DGRs_found_dict

        if not len(dgrs_dict):
            raise ConfigError("Compute DGR variability profile function speaking: There are no DGRs to compute in-sample variability :/")

        sample_names = self.samples_artifact.samples()
        num_samples = len(sample_names)

        # let the user know what is going on
        msg = (f"Now anvi'o will compute in-sample activity of {PL('DGR VR', len(self.DGRs_found_dict))} "
            f"across {PL('sample', num_samples)}. Brace yourself and please note that this can "
            "take a very long time since for each sample, anvi'o will go through each short read to search for ever variable region "
            "sequence/s per DGR. You can always skip this step and search for individual primers "
            "listed in the output file using the program `anvi-search-primers` with the parameter "
            "`--min-remainder-length` set to the length of the consensus template region of the DGRs' VR you are interested in "
            "and the flag `--only-report-remainders` to explore the variable region activity of that one variable region "
            "manually. Maybe this makes no sense? See the documentation for `anvi-report-dgrs` (and hope for "
            "the best)")
        self.run.warning(None, header="PERFORMANCE NOTE", lc="yellow")

        if num_samples > self.num_threads:
            self.run.info_single(f"You have {PL('sample', num_samples)} but {PL('thread', self.num_threads)}. Therefore, not all samples will be processed "
                                f"in parallel. Just an FYI. {msg}.", level=0, nl_before=1)
        elif self.num_threads > num_samples:
            self.run.info_single(f"You have {PL('sample', num_samples)} but {PL('thread', self.num_threads)}. Since only samples are run in "
                                f"parallel, the additional {PL('thread', self.num_threads - num_samples)} you have there are not really "
                                f"useful for anything. Just an FYI. {msg}.", level=0, nl_before=1)
            self.num_threads = num_samples
        else:
            self.run.info_single(f"{msg}.", level=0, nl_before=1)

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

        # output the final primers dictionary
        self.run.info_single("Computing the Variable Regions Primers and creating a 'DGR_Primers_used_for_VR_diversity.tsv' file.", nl_before=1)
        self.print_primers_dict_to_csv(primers_dict)

        ##################
        # MULTITHREADING #
        ##################

        # setup the input/output queues
        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        self.dgr_activity = []
        # create directory for Primer matches
        primer_folder= os.path.join(self.output_directory, "PRIMER_MATCHES")

        # put all the sample names in our input queue
        for sample_name in sample_names:
            input_queue.put(sample_name)

        # engage the proletariat, our hard-working wage-earner class
        workers = []
        for i in range(self.num_threads):
            worker = multiprocessing.Process(target=DGR_Finder.compute_dgr_variability_profiling_per_vr,
                                            args=(input_queue,
                                                output_queue,
                                                self.samples_artifact,
                                                primers_dict,
                                                primer_folder),
                                            kwargs=({'progress': self.progress if self.num_threads == 1 else progress_quiet
                                                }))
            workers.append(worker)
            worker.start()

        # monitor progress
        self.progress.new('DGR variability profile', progress_total_items=num_samples)
        if self.num_threads > 1:
            self.progress.update(f"Processing {PL('sample', num_samples)} and {PL('primer', len(primers_dict))} in {PL('thread', self.num_threads)}.")

        num_samples_processed = 0
        while num_samples_processed < num_samples:
            try:
                sample_finished_processing = output_queue.get()
                if anvio.DEBUG:
                    self.progress.reset()
                    self.run.info_single(f"Sample {sample_finished_processing} has finished processing.", nl_before=1)
                num_samples_processed += 1
                self.progress.increment(increment_to=num_samples_processed)

                if self.num_threads > 1:
                    if num_samples_processed < num_samples:
                        self.progress.update(f"Samples processed: {num_samples_processed} of {num_samples}. Still working ...")
                    else:
                        self.progress.update("All done!")
            except KeyboardInterrupt:
                self.run.info_single("Received kill signal, terminating all processes... Don't believe anything you see "
                                    "below this and destroy all the output files with fire.", nl_before=1)
                break

        if self.num_threads > 1:
            self.progress.end()

        # always double-tap?
        for worker in workers:
            worker.terminate()

        ######################
        # END MULTITHREADING #
        ######################



    def print_primers_dict_to_csv(self, primers_dict):
        """
        Turn the primers dictionary into a tsv file for users to visualise post analysis.
        This function will create a TSV file named 'DGR_Primers_used_for_VR_diversity.tsv' in the output directory.

        Parameters
        ==========
        primers_dict : dict
            Dictionary of primer metadata
        Returns
        =======
        None
            Creates a TSV file with primer details for downstream review.
        """

        self.run.info_single("Storing the primers used to calculate VR diversity in 'DGR_Primers_used_for_VR_diversity.tsv'.", nl_before=1)

        output_directory_path = self.output_directory
        output_path= os.path.join(output_directory_path, "DGR_Primers_used_for_VR_diversity.tsv")

        # define the header for the TSV file
        csv_header = ['Primer_ID', 'Sample_ID', 'No_SNV_Primer', 'Initial_Primer', 'Masked_Primer', 'Whole_Primer']

        # Open the TSV file in write mode
        with open(output_path, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(csv_header)  # Write the header row

            # iterate through the dictionary and write each primer's information to the TSV file
            for primer_name, samples in primers_dict.items():
                for sample_id, primer_info in samples.items():
                    used_original_primer = primer_info['used_original_primer']
                    initial_primer = primer_info['initial_primer_sequence']
                    masked_primer = primer_info['vr_masked_primer']
                    whole_primer = primer_info['primer_sequence']

                    writer.writerow([
                    primer_name,
                    sample_id,
                    used_original_primer,
                    initial_primer,
                    masked_primer,
                    whole_primer
                    ])
        return



    def process_dgr_data_for_HTML_summary(self):
        """
        Take everything that is known, turn them into data that can be used from Django templates.

        A lot of ugly/frightening things happening here to prepare coordinates for SVG objects to be displayed
        or store boolean variables to use the Django template engine effectively. IF YOU DON'T LIKE
        IT DON'T LOOK AT IT. IT MIGHT MAKE YOU CRY

        Returns
        =======
        None
            Generates a directory (`summary_html_output`) containing
            the HTML summary and supporting files.

        """
        if not len(self.DGRs_found_dict):
            self.run.warning("No DGRs were found so no HTML file will be written :(", header="No HTML OUTPUT")
            return

        # in which we will store all the static HTML output related stuff
        self.summary = {}

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.summary['meta'] = {'summary_type': 'dgrs',
                                'num_dgrs': len(self.DGRs_found_dict),
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

        # track DGRs without TR gene context
        dgrs_without_tr_genes = []

        for dgr_key, dgr_data in dgrs_dict.items():
            # assuming dgr_key itself is the dgr_id or a dictionary containing it
            dgr_id = dgr_key

            self.summary['dgrs'][dgr_id] = {'dgr_data': copy.deepcopy(dgr_data), 'tr_genes': {}, 'vr_genes': {}}

            # handle genomic context recovery
            if not self.skip_recovering_genomic_context:
                # Get the genomic context for this DGR (might be empty if no genes found)
                dgr_context = self.genomic_context_surrounding_dgrs.get(dgr_id, {})

                # deep copy the genomic context for TR and VR genes
                tr_genes = {gene_id: gene_info
                    for gene_id, gene_info in self.genomic_context_surrounding_dgrs.get(dgr_id, {}).items()
                    if isinstance(gene_id, int)}

                # then we will learn about these so we can transform the coordinates of anything we wish
                # to display in the output
                # if we have genes we can display them else not
                if tr_genes:
                    # get the start and stop positions of the first and last genes in tr_genes
                    genomic_context_start_tr = min(gene['start'] for gene in tr_genes.values()) - 100
                    genomic_context_end_tr = max(gene['stop'] for gene in tr_genes.values()) + 100
                else:
                    genomic_context_start_tr = int(dgr_data['TR_start_position']) - 500
                    genomic_context_end_tr = int(dgr_data['TR_end_position']) + 500
                    dgrs_without_tr_genes.append(dgr_id)
                    if anvio.DEBUG:
                        self.run.warning(
                            f"No TR gene calls found around {dgr_id}. "
                            f"Using TR coordinates {genomic_context_start_tr}-{genomic_context_end_tr} instead."
                            )

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
                    # transform start and stop coordinates for the gene
                    gene['start_tr_g'] = (gene['start'] - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length
                    gene['stop_tr_g'] = (gene['stop'] - genomic_context_start_tr) / (genomic_context_end_tr - genomic_context_start_tr) * new_context_length

                    if (gene['stop_tr_g'] - gene['start_tr_g']) < default_gene_arrow_width:
                        # if the transformed length of the gene is smaller than the default arrow width
                        gene_arrow_width = gene['stop_tr_g'] - gene['start_tr_g']
                        gene['stop_tr_g'] = gene['start_tr_g']
                        gene['TRW'] = 0
                    else:
                        gene_arrow_width = default_gene_arrow_width
                        gene['TRW'] = (gene['stop_tr_g'] - gene['start_tr_g']) - gene_arrow_width

                    # determine color based on gene functions first
                    if gene['functions']:
                        gene['has_functions'] = True
                        gene['COLOR'] = '#008000'  # Green for genes with functions
                    else:
                        gene['has_functions'] = False
                        gene['COLOR'] = '#c3c3c3'  # Grey for genes without functions

                    # check for hmm_id and gene_id match if provided in the dgr_data
                    try:
                        hmm_id = int(dgr_data.get('HMM_gene_callers_id'))
                        current_gene_id = int(gene.get('gene_callers_id'))
                    except (TypeError, ValueError):
                        hmm_id = None
                        current_gene_id = None

                    if hmm_id is not None and current_gene_id is not None and hmm_id == current_gene_id:
                        gene['COLOR'] = '#c366e8'  # Purple if there's a match based on hmm_id and gene_id

                    # additional transformations for coordinates
                    gene['TRX'] = gene['start_tr_g']
                    gene['TCX'] = (gene['start_tr_g'] + (gene['stop_tr_g'] - gene['start_tr_g']) / 2)
                    gene['TGY'] = gene['TRX'] + gene['TRW'] + gene_arrow_width
                    gene['TGTRANS'] = gene['TRX'] + gene['TRX'] + gene['TRW'] + gene_arrow_width
                    gene['TRX_TRW'] = gene['TRX'] + gene['TRW'] - 0.5

                    # append the gene to the summary
                    self.summary['dgrs'][dgr_id]['tr_genes'][gene_id] = gene

                for vr_key, vr_data in dgr_data.get('VRs', {}).items():
                    vr_id = vr_key
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'] = self.summary['dgrs'][dgr_id]['dgr_data'].get('VRs', {})

                    # extract VR genes (string keys like 'VR_001')
                    vr_genes = dgr_context.get(vr_key, {})

                    # calculate genomic context for VR genes
                    if vr_genes:
                        genomic_context_start_vr = min(
                            vr_genes[gene]['start']
                            for gene in vr_genes
                        ) - 100

                        genomic_context_end_vr = max(
                            vr_genes[gene]['stop']
                            for gene in vr_genes
                        ) + 100
                    else:
                        genomic_context_start_vr = int(vr_data['VR_start_position']) - 500
                        genomic_context_end_vr = int(vr_data['VR_end_position']) + 500

                    # transform coordinates for VR data
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

                        # determine color based on gene functions first
                        if gene['functions']:
                            gene['has_functions'] = True
                            gene['COLOR'] = '#008000'  # green for genes with functions
                        else:
                            gene['has_functions'] = False
                            gene['COLOR'] = '#c3c3c3'  # grey for genes without functions

                        # check for hmm_id and gene_id match if provided in the dgr_data
                        try:
                            hmm_id = int(dgr_data.get('HMM_gene_callers_id'))
                            current_gene_id = int(gene.get('gene_callers_id'))
                        except (TypeError, ValueError):
                            hmm_id = None
                            current_gene_id = None

                        if hmm_id is not None and current_gene_id is not None and hmm_id == current_gene_id:
                            gene['COLOR'] = '#c366e8'  # purple if there's a match based on hmm_id and gene_id

                        gene['VRX'] = gene['start_vr_g']
                        gene['VCX'] = (gene['start_vr_g'] + (gene['stop_vr_g'] - gene['start_vr_g']) / 2)
                        gene['VGY'] = gene['VRX'] + gene['VRW'] + gene_arrow_width
                        gene['VGTRANS'] = gene['VRX'] + gene['VRX'] + gene['VRW'] + gene_arrow_width
                        gene['VRX_VRW'] = gene['VRX'] + gene['VRW'] - 0.5 # <-- minus 0.5 makes the arrow nicely cover the rest of the gene

                    # append each transformed VR gene individually to the summary
                    self.summary['dgrs'][dgr_id]['dgr_data']['VRs'][vr_id]['vr_genes'] = vr_genes

                    # Store individual VR genes in vr_genes dict (for compatibility)
                    for gene_id, gene in vr_genes.items():
                        self.summary['dgrs'][dgr_id]['vr_genes'][gene_id] = gene

                    # Store output file paths for DGRs and VRs
                    if dgr_id not in self.summary['files']:
                        self.summary['files'][dgr_id] = {}

                    # Add TR files if TR genes exist
                    if tr_genes:
                        self.summary['files'][dgr_id]['tr_genes'] = os.path.join('PER_DGR', dgr_id, 'SURROUNDING-GENES.txt')
                        self.summary['files'][dgr_id]['functions'] = os.path.join('PER_DGR', dgr_id, 'SURROUNDING-FUNCTIONS.txt')

                    # Add VR files
                    self.summary['files'][dgr_id][f'{vr_id}_vr_genes'] = os.path.join('PER_DGR', vr_id, 'SURROUNDING-GENES.txt')
                    self.summary['files'][dgr_id][f'{vr_id}_vr_functions'] = os.path.join('PER_DGR', vr_id, 'SURROUNDING-FUNCTIONS.txt')


        # summary of DGRs without TR gene context
        if dgrs_without_tr_genes:
            self.run.warning(f"{PL('DGR', len(dgrs_without_tr_genes))} had no gene calls around their Template Region "
                            f"(using TR coordinates as fallback): {', '.join(dgrs_without_tr_genes)}. "
                            "Use '--debug' to see individual messages.",
                            header="DGRs WITHOUT TR GENE CONTEXT")

        # ensure the destination directory does not exist before generating the summary HTML
        destination_dir = 'summary_html_output'
        if os.path.exists(destination_dir):
            shutil.rmtree(destination_dir)

        summary_html_output = SummaryHTMLOutput(self.summary)
        summary_html_output.generate('summary_html_output')

        return



    def parameter_output_sheet(self):
        """
        This function creates a tsv tabular format of all the parameters the user input in the current run.

        Returns
        =======
        : tsv
            A tsv tabular file containing the template and variable regions

        """
        output_path_parameters = os.path.join(self.output_directory, "Parameters_used_in_DGRs_found.tsv")
        with open(output_path_parameters, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter='\t')
            headers = ["Parameter", "Value"]
            csv_writer.writerow(headers)

            parameters = [
                ("Contig.db", self.contigs_db_path),
                ("Profile.db", self.profile_db_path),
                ("Temporary Directory", self.temp_dir if self.temp_dir else None),
                ("Output Directory", self.output_directory if self.output_directory else "default"),
                ("Word Size of BLAST", self.word_size if self.word_size else "8"),
                ("Number of Threads for BLAST", self.num_threads if self.num_threads else "1"),
                ("Skip 'Ns'", self.skip_Ns if self.skip_Ns else "FALSE"),
                ("Skip '-'", self.skip_dashes if self.skip_dashes else "FALSE"),
                ("Number of Mismatches", self.number_of_mismatches if self.number_of_mismatches else "7"),
                ("Initial Mismatch Bias Threshold", self.initial_mismatch_bias_threshold if self.initial_mismatch_bias_threshold else "0.6"),
                ("Max Non-Dominant Base Threshold", self.max_non_dominant if self.max_non_dominant else "1"),
                ("Minimum VR Length", self.minimum_vr_length if self.minimum_vr_length else "50"),
                ("Allow Any Base", self.allow_any_base or "FALSE"),
                ("SNV Window Size", self.snv_window_size if self.snv_window_size else "50"),
                ("SNV Window Step", self.snv_window_step if self.snv_window_step else "10"),
                ("Minimum Mismatching Base Types in VR", self.min_base_types_vr if self.min_base_types_vr else "2"),
                ("Minimum Base Types in VR", self.min_base_types_vr if self.min_base_types_vr else "2"),
                ("Minimum Base Types in TR", self.min_base_types_tr if self.min_base_types_tr else "2"),
                ("Number of imperfect tandem repeats", self.numb_imperfect_tandem_repeats or "10"),
                ("Repeat motif coverage", self.repeat_motif_coverage or "0.8"),
                ("Variable Buffer Length", self.variable_buffer_length if self.variable_buffer_length else "35"),
                ("SNV Matching Proportion", self.snv_matching_proportion or "0.25/0.30"),
                ("SNV Codon Position", self.snv_codon_position or "0.33"),
                ("SNV Density", self.minimum_snv_density or "0.1"),
                ("Gene caller", self.gene_caller_to_consider_in_context if self.gene_caller_to_consider_in_context else "prodigal"),
                ("Number of genes to consider in context", self.num_genes_to_consider_in_context or "3"),
                ("Skip Recovering Genomic Context", self.skip_recovering_genomic_context or "FALSE"),
                ("HMMs Provided to Search through", self.hmm or "Reverse_Transcriptase"),
                ("Discovery mode", self.discovery_mode or "FALSE"),
                ("Collections Mode", self.collections_mode or "FALSE"),
                ("Collections Given", self.collections_given or "None"),
                ("Parameter Outputs", self.parameter_outputs or "FALSE"),
                ("Just Do It", self.just_do_it or "FALSE"),
                ("Verbose", self.verbose or "FALSE"),
                ("Samples.txt", self.samples_txt or "None"),
                ("Initial Primer Length", self.initial_primer_length or "12"),
                ("Variable Region Primer Length", self.whole_primer_length or "65"),
                ("Skip Compute DGR Variability Profiling", self.skip_compute_DGR_variability_profiling or "FALSE"),
            ]

            csv_writer.writerows(parameters)
        return



    def process(self,args):
        """
        Here we process all of the functions in our class and call upon different functions depending on the inputs used,
        to run the full DGR analysis workflow.

        Parameters
        =======
        args : argparse.Namespace
            Command-line arguments controlling optional outputs.
        """
        self.sanity_check()

        # do we have a previously computed list of dgrs to focus for dgr activity calculations?
        if self.pre_computed_dgrs_path:
            self.run.warning("Anvi'o is taking a shortcut to calculate DGR activity using the DGRs you have "
                            "provided in the 'pre computed dgrs' file. There is nothing for you to be concerned about "
                            "-- except the fact that some very fancy coding is at play here and catastrophic failures "
                            "are not a remote possibility. Still, anvi'o prints this message in green for positive "
                            "vibes 🥲", header="🌈 PRE-COMPUTED DGRS FOUND 🌈", lc="green")

            self.populate_dgrs_dict_from_input_file()
            self.compute_dgr_variability_profiling()

            # WE'RE DOnE HERE. DoNe.
            return

        # Run detection based on the selected mode
        self.run.info('Detection mode', self.detection_mode)

        activity_hits = None
        homology_hits = None

        # Activity-based detection (requires profile.db with SNV data)
        if self.detection_mode in ('activity', 'both'):
            self.run.info_single("Running activity-based DGR detection...", nl_before=1, mc='green')
            self.get_blast_results()
            self.process_blast_results(vr_in_query=True, apply_snv_filters=True)

            # Save activity hits for potential merging
            if self.collections_mode:
                activity_hits = copy.deepcopy(self.merged_mismatch_hits)
            else:
                activity_hits = copy.deepcopy(self.mismatch_hits)

        # Homology-based detection (uses RT HMM windows)
        if self.detection_mode in ('homology', 'both'):
            self.run.info_single("Running homology-based DGR detection...", nl_before=1, mc='green')

            if self.collections_mode:
                # Run homology BLAST for each bin in the collection
                homology_blast_outputs = self.process_collections_mode_homology()

                if homology_blast_outputs:
                    # Store bin names list for process_blast_results to iterate over
                    # (only bins that had successful BLAST runs)
                    self.bin_names_list = list(homology_blast_outputs.keys())

                    # Process BLAST results for all bins (homology mode)
                    self.process_blast_results(vr_in_query=False, apply_snv_filters=False, mode='homology')
                    homology_hits = copy.deepcopy(self.merged_mismatch_hits)
                else:
                    self.run.warning("No bins had successful homology-based BLAST runs.",
                                   header="NO HOMOLOGY RESULTS")
                    homology_hits = defaultdict(lambda: defaultdict(dict))
            else:
                # Standard (non-collections) homology mode
                homology_hits = self.process_homology_mode()

        # Merge results if running both modes
        if self.detection_mode == 'both':
            if activity_hits and homology_hits:
                merged = self.merge_detection_results(activity_hits, homology_hits)
                if self.collections_mode:
                    self.merged_mismatch_hits = merged
                else:
                    self.mismatch_hits = merged
            elif activity_hits:
                # Only activity hits found
                if self.collections_mode:
                    self.merged_mismatch_hits = activity_hits
                else:
                    self.mismatch_hits = activity_hits
            elif homology_hits:
                # Only homology hits found
                if self.collections_mode:
                    self.merged_mismatch_hits = homology_hits
                else:
                    self.mismatch_hits = homology_hits
            else:
                # No hits from either mode
                self.run.warning("No DGR candidates were found by either activity or homology detection modes.",
                               header="NO DGRs FOUND")
                if self.collections_mode:
                    self.merged_mismatch_hits = defaultdict(lambda: defaultdict(dict))
                else:
                    self.mismatch_hits = defaultdict(lambda: defaultdict(dict))

        self.filter_for_best_VR_TR()
        if args.parameter_output:
            self.run.info_single("Writing to Parameters used file.", nl_before=1)
            self.parameter_output_sheet()
        self.get_gene_info()
        self.get_hmm_info()
        self.recover_genomic_context_surrounding_dgrs()
        self.report_genomic_context_surrounding_dgrs()
        self.create_found_tr_vr_tsv()
        self.compute_dgr_variability_profiling()
        self.process_dgr_data_for_HTML_summary()

        return
