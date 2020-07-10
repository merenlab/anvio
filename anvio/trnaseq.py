# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for tRNA-seq dataset operations.

    anvi-trnaseq is the default client using this module.
"""

import itertools
import multiprocessing
import os
import pysam

from collections import OrderedDict

import anvio
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.trnaidentifier as trnaidentifier
import anvio.utils as utils

from anvio.constants_package.trnaseq import THREEPRIME_VARIANTS, TRNA_FEATURE_NAMES
from anvio.dbops_package.trnaseq import TRNASeqDatabase
from anvio.errors import ConfigError, FilesNPathsError
from anvio.sequence import SequenceDereplicator

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"

class TRNASeqDataset:
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_fasta_path = A('fasta_file')
        self.project_name = A('project_name')
        self.output_dir = A('output_dir')
        self.write_uniquing_report = A('uniquing_report')
        self.write_normalization_reports = A('normalization_reports')
        self.max_possible_alignments = A('max_possible_alignments')
        self.description_file_path = A('description')
        self.skip_fasta_check = A('skip_fasta_check')
        self.num_threads = A('num_threads')
        self.write_buffer_size = A('write_buffer_size')

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using --fasta-file or -f.")
        if not self.project_name:
            raise ConfigError("Please set a project name using --project-name or -n.")
        if not self.output_dir:
            raise ConfigError("Please provide a non-existent output directory using --output-dir or -o.")

        self.input_seq_count = 0

        self.trnaseq_db_path = os.path.join(self.output_dir, self.project_name + "-TRNASEQ.db")
        self.nontrna_report_path = os.path.join(self.output_dir, self.project_name + "-unique_nontrna.txt")
        self.NONTRNA_REPORT_HEADER = ["representative_id", "member_ids", "member_count", "sequence"]

        if self.write_uniquing_report:
            self.uniquing_report_path = os.path.join(self.output_dir, self.project_name + "-uniquing_report.txt")
        else:
            self.uniquing_report_path = None
        self.UNIQUING_REPORT_HEADER = ["representative_id", "member_ids"]

        self.profiled_unique_seq_count = 0
        self.trna_count = 0
        self.unique_trna_count = 0
        self.trna_containing_anticodon_count = 0
        self.mature_trna_count = 0
        self.trna_with_one_to_five_extra_fiveprime_bases_count = 0
        self.trna_with_more_than_five_extra_fiveprime_bases_count = 0
        self.trna_with_extrapolated_fiveprime_feature_count = 0
        self.trna_with_threeprime_cca_count = 0
        self.trna_with_threeprime_cc_count = 0
        self.trna_with_threeprime_c_count = 0
        self.trna_with_threeprime_nca_cna_ccn_count = 0
        self.trna_with_threeprime_ccan_ccann_count = 0
        self.trna_with_extra_threeprime_bases_count = 0

        self.derep_mapping_temp_dir = os.path.join(self.output_dir, "DEREP_MAPPING_TEMP")

        if self.write_normalization_reports:
            self.fiveprime_normalization_report_path = os.path.join(self.output_dir, self.project_name + "-fiveprime_normalization_report.txt")
            self.threeprime_normalization_report_path = os.path.join(self.output_dir, self.project_name + "-threeprime_normalization_report.txt")
            self.subsequence_report_path = os.path.join(self.output_dir, self.project_name + "-normalized_subsequence_report.txt")
        else:
            self.fiveprime_normalization_report_path = None
            self.threeprime_normalization_report_path = None
            self.subsequence_report_path = None
        self.FIVEPRIME_NORMALIZATION_REPORT_HEADER = ["representative_id", "member_id", "fiveprime_sequence", "threeprime_sequence", "input_seq_count"]
        self.THREEPRIME_NORMALIZATION_REPORT_HEADER = ["representative_id", "member_id", "threeprime_sequence", "input_seq_count"]
        self.SUBSEQUENCE_REPORT_HEADER = ["representative_id", "threeprime_member_ids", "interior_member_ids"]

        self.tiny_threeprime_trna_ids_path = os.path.join(self.output_dir, self.project_name + "-tiny_threeprime_trna_ids.txt")

        self.seed_count = 0

        self.log_path = os.path.join(self.output_dir, "log.txt")


    def check_programs(self):
        utils.is_program_exists('bowtie2')


    def check_output_dir(self):
        if os.path.exists(self.output_dir):
            raise ConfigError("The directory that was specified by --output-dir or -o, %s, already exists. "
                              "Please try again with a non-existent directory." % self.output_dir)
        os.mkdir(self.output_dir)
        filesnpaths.is_output_dir_writable(self.output_dir)
        self.output_dir = os.path.abspath(self.output_dir)


    def get_output_file_path(self, file_name, delete_if_exists=False):
        output_file_path = os.path.join(self.output_dir, file_name)

        if delete_if_exists:
            if os.path.exists(output_file_path):
                os.remove(output_file_path)

        return output_file_path


    def sanity_check(self):

        self.check_programs()

        self.check_output_dir()

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path)).read()
        else:
            self.description = None
        self.run.info("Description", self.description_file_path if self.description_file_path else "No description given")

        self.run.log_path = self.log_path

        if self.max_possible_alignments < 2:
            raise ConfigError("The maximum number of possible alignments in Bowtie2 clustering "
                                "of uniqued, normalized sequences to themselves cannot be less than 2, "
                                "as one of these alignments will be the sequence mapped to itself. Try again!")

        if not 1 < self.num_threads < multiprocessing.cpu_count():
            ConfigError("The number of threads to use must be a positive integer "
                        "less than or equal to %d. Try again!" % multiprocessing.cpu_count())

        if self.write_buffer_size < 1:
            ConfigError("The write buffer size must be a positive integer. Try again!")

        self.run.info("Input FASTA file", self.input_fasta_path)
        utils.check_fasta_id_uniqueness(self.input_fasta_path)
        if not self.skip_fasta_check:
            self.progress.new("Checking input FASTA defline format")
            utils.check_fasta_id_formatting(self.input_fasta_path)
            self.progress.end()
            self.run.info_single("FASTA deflines were found to be anvi'o-compliant", mc='green')


    def generate_trnaseq_db(self):
        meta_values = {'project_name': self.project_name,
                       'description': self.description if self.description else '_No description is provided_',
                       'max_possible_alignments': self.max_possible_alignments}
        TRNASeqDatabase(self.trnaseq_db_path, quiet=False).create(meta_values)


    def profile_trna(self):

        def get_unique_input_seqs():
            fasta = fastalib.SequenceSource(self.input_fasta_path)
            ids = []
            seqs = []
            while next(fasta):
                ids.append(fasta.id)
                seqs.append(fasta.seq)
            self.input_seq_count = len(ids)
            clusters = SequenceDereplicator(ids, seqs).full_length_dereplicate()

            unique_member_ids = [t[1] for t in clusters]
            unique_representative_ids = [l[0] for l in unique_member_ids]
            unique_input_seqs = [t[0] for t in clusters]

            if self.write_uniquing_report:
                with open(self.uniquing_report_path, 'w') as uniquing_report:
                    uniquing_report.write("\t".join(self.UNIQUING_REPORT_HEADER) + "\n")
                    for rep_id, member_ids in zip(unique_representative_ids, unique_member_ids):
                        uniquing_report.write(rep_id + "\t")
                        uniquing_report.write(",".join(member_ids) + "\n")

            return unique_representative_ids, unique_member_ids, unique_input_seqs


        def write_results():

            # List of entries for each table
            trnaseq_sequences_table_entries = []
            trnaseq_info_table_entries = []
            trnaseq_features_table_entries = []
            trnaseq_unconserved_table_entries = []
            trnaseq_unpaired_table_entries = []

            retrieved_profile_count = 0

            while retrieved_profile_count < self.profiled_unique_seq_count:
                # Retrieve profiles from the output queue and write tRNA profiles to the database.
                trna_profile = output_queue.get()
                retrieved_profile_count += 1

                output_id = trna_profile.name
                output_seq = trna_profile.input_seq

                if len(trna_profile.features) < num_features_through_t_loop:
                    replicate_ids = replicate_dict[output_id]
                    nontrna_report.write(output_id + "\t"
                                         + ",".join(replicate_ids) + "\t"
                                         + str(len(replicate_ids)) + "\t"
                                         + output_seq + "\n")
                    continue

                self.unique_trna_count += 1
                output_seq_length = len(output_seq)

                replicate_ids = replicate_dict[output_id]
                num_replicates = len(replicate_ids)
                self.trna_count += num_replicates

                # Recover nucleotides that did not fit expectation,
                # either by not being the expected nucleotide or type of nucleotide
                # or by not base pairing in a stem.
                unconserved_info = trna_profile.get_unconserved_positions()
                unpaired_info = trna_profile.get_unpaired_positions()

                if trna_profile.anticodon_seq:
                    self.trna_containing_anticodon_count += num_replicates
                if trna_profile.has_complete_feature_set:
                    self.mature_trna_count += num_replicates
                if trna_profile.num_in_extrapolated_fiveprime_feature > 0:
                    self.trna_with_extrapolated_fiveprime_feature_count += num_replicates

                if trna_profile.num_extra_threeprime > 0:
                    self.trna_with_threeprime_ccan_ccann_count += num_replicates
                elif trna_profile.acceptor.string == 'CCA':
                    self.trna_with_threeprime_cca_count += num_replicates
                elif trna_profile.acceptor.string == 'CC':
                    self.trna_with_threeprime_cc_count += num_replicates
                elif trna_profile.acceptor.string == 'C':
                    self.trna_with_threeprime_c_count += num_replicates
                else:
                    self.trna_with_threeprime_nca_cna_ccn_count += num_replicates

                if trna_profile.num_extra_fiveprime > 5:
                    self.trna_with_more_than_five_extra_fiveprime_bases_count += num_replicates
                elif trna_profile.num_extra_fiveprime > 0:
                    self.trna_with_one_to_five_extra_fiveprime_bases_count += num_replicates

                trnaseq_sequences_table_entries.append((output_id, num_replicates, output_seq))

                # The alpha and beta regions of the D loop vary in length.
                # Record their start and stop positions in the sequence if they were profiled.
                # These are included in the info table rather than the feature table,
                # because they are subfeatures of the D loop feature,
                # and the positions of the features in the feature table are not overlapping.
                alpha_start = trna_profile.alpha_start if trna_profile.alpha_start else '??'
                alpha_stop = trna_profile.alpha_stop - 1 if trna_profile.alpha_stop else '??'
                beta_start = trna_profile.beta_start if trna_profile.beta_start else '??'
                beta_stop = trna_profile.beta_stop - 1 if trna_profile.beta_stop else '??'

                trnaseq_info_table_entries.append((
                    output_id,
                    trna_profile.has_complete_feature_set,
                    trna_profile.anticodon_seq,
                    trna_profile.anticodon_aa,
                    output_seq_length,
                    # Zero-based start index of identified tRNA features within the input sequence.
                    output_seq_length - len(trna_profile.profiled_seq),
                    # Stop index of features (real stop position, not Pythonic stop index for slicing).
                    output_seq_length - trna_profile.num_extra_threeprime - 1,
                    trna_profile.num_conserved,
                    trna_profile.num_unconserved,
                    trna_profile.num_paired,
                    trna_profile.num_unpaired,
                    trna_profile.num_in_extrapolated_fiveprime_feature,
                    trna_profile.num_extra_fiveprime,
                    trna_profile.num_extra_threeprime,
                    alpha_start,
                    alpha_stop,
                    beta_start,
                    beta_stop))

                trnaseq_features_table_entries.append(
                    (output_id, )
                    # When tRNA features were not found at the 5' end of the read,
                    # their start and stop indices also were not found.
                    + tuple(['?' * 2 for _ in range((len(TRNA_FEATURE_NAMES) - len(trna_profile.features)))]) * 2
                    + tuple(itertools.chain(*zip(
                        [str(feature.start_index) if hasattr(feature, 'start_index')
                            else ','.join(map(str, feature.start_indices))
                            for feature in trna_profile.features],
                        # Convert Pythonic stop index for slicing to real stop position of feature.
                        [str(feature.stop_index - 1) if hasattr(feature, 'stop_index')
                            else ','.join(map(str, [stop_index - 1 for stop_index in feature.stop_indices]))
                            for feature in trna_profile.features]))))

                for unconserved_tuple in unconserved_info:
                    trnaseq_unconserved_table_entries.append((output_id, ) + unconserved_tuple)

                for unpaired_tuple in unpaired_info:
                    trnaseq_unpaired_table_entries.append((output_id, ) + unpaired_tuple)

            if len(trnaseq_sequences_table_entries) > 0:
                trnaseq_db.db._exec_many(
                    '''INSERT INTO %s VALUES (%s)'''
                    % ('sequences', ','.join('?' * len(tables.trnaseq_sequences_table_structure))),
                    trnaseq_sequences_table_entries)
                trnaseq_db.db._exec_many(
                    '''INSERT INTO %s VALUES (%s)'''
                    % ('basic_info', ','.join('?' * len(tables.trnaseq_info_table_structure))),
                    trnaseq_info_table_entries)
                trnaseq_db.db._exec_many(
                    '''INSERT INTO %s VALUES (%s)'''
                    % ('features', ','.join('?' * len(tables.trnaseq_features_table_structure))),
                    trnaseq_features_table_entries)
                trnaseq_db.db._exec_many(
                    '''INSERT INTO %s VALUES (%s)'''
                    % ('unconserved_nucleotides', ','.join('?' * len(tables.trnaseq_unconserved_table_structure))),
                    trnaseq_unconserved_table_entries)
                trnaseq_db.db._exec_many(
                    '''INSERT INTO %s VALUES (%s)'''
                    % ('unpaired_nucleotides', ','.join('?' * len(tables.trnaseq_unpaired_table_structure))),
                    trnaseq_unpaired_table_entries)


        def set_meta_values():
            trnaseq_db.db.set_meta_value('num_input_reads_processed', self.input_seq_count)
            trnaseq_db.db.set_meta_value('num_trna_reads', self.trna_count)
            trnaseq_db.db.set_meta_value('num_unique_trna_seqs', self.unique_trna_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_containing_anticodon', self.trna_containing_anticodon_count)
            trnaseq_db.db.set_meta_value('num_mature_trna_reads', self.mature_trna_count)
            trnaseq_db.db.set_meta_value('num_trna_with_one_to_five_extra_fiveprime_bases', self.trna_with_one_to_five_extra_fiveprime_bases_count)
            trnaseq_db.db.set_meta_value('num_trna_with_more_than_five_extra_fiveprime_bases', self.trna_with_more_than_five_extra_fiveprime_bases_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_with_extrapolated_fiveprime_feature', self.trna_with_extrapolated_fiveprime_feature_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_cca', self.trna_with_threeprime_cca_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_cc', self.trna_with_threeprime_cc_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_c', self.trna_with_threeprime_c_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_nca_cna_ccn', self.trna_with_threeprime_nca_cna_ccn_count)
            trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_ccan_ccann', self.trna_with_threeprime_ccan_ccann_count)


        def output_run_info():
            print() # without this, run info will write on the same line as progress
            self.run.info("Reads processed", self.input_seq_count)
            self.run.info("Reads identified as tRNA", self.trna_count)
            self.run.info("Unique tRNA sequences identified", self.unique_trna_count)
            self.run.info("tRNA reads containing anticodon", self.trna_containing_anticodon_count)
            self.run.info("tRNA reads spanning acceptor stem (\"mature\")", self.mature_trna_count)
            self.run.info("tRNA reads with 1-5 extra 5' bases", self.trna_with_one_to_five_extra_fiveprime_bases_count)
            self.run.info("tRNA reads with more than 5 extra 5' bases", self.trna_with_more_than_five_extra_fiveprime_bases_count)
            self.run.info("tRNA reads containing extrapolated 5' feature", self.trna_with_extrapolated_fiveprime_feature_count)
            self.run.info("tRNA reads ending in 3'-CCA", self.trna_with_threeprime_cca_count)
            self.run.info("tRNA reads ending in 3'-CC", self.trna_with_threeprime_cc_count)
            self.run.info("tRNA reads ending in 3'-C", self.trna_with_threeprime_c_count)
            self.run.info("tRNA reads ending in 3'-NCA/CNA/CCN", self.trna_with_threeprime_nca_cna_ccn_count)
            self.run.info("tRNA reads ending in 3'-CCAN/CCANN", self.trna_with_threeprime_ccan_ccann_count)


        ########################
        # profile_trna main body
        ########################

        self.progress.new("Profiling input sequences for tRNA features")

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)

        unique_representative_ids, unique_member_ids, unique_input_seqs = get_unique_input_seqs()
        unique_input_seq_count = len(unique_representative_ids)

        # Reads are identified as tRNA if they have the minimum set of features from the 3' acceptor through the T loop.
        num_features_through_t_loop = TRNA_FEATURE_NAMES[::-1].index('t_loop') + 1
        num_features_through_t_stem = num_features_through_t_loop + 1

        # Determine the points when to write to the database.
        # Write points are the cumulative number of sequences processed.
        num_chunks = unique_input_seq_count // self.write_buffer_size
        if num_chunks == 0:
            num_chunks += 1
            write_points = [unique_input_seq_count]
        else:
            write_points = [self.write_buffer_size * (i + 1) for i in range(num_chunks)]
            if unique_input_seq_count % num_chunks > 0:
                num_chunks += 1
                write_points.append(unique_input_seq_count)
        write_point_index = 0

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        processes = [multiprocessing.Process(target=trnaidentifier.profile_wrapper, args=(input_queue, output_queue))
                     for _ in range(self.num_threads)]
        for process in processes:
            process.start()


        replicate_dict = {}
        nontrna_report = open(self.nontrna_report_path, 'w')
        nontrna_report.write("\t".join(self.NONTRNA_REPORT_HEADER) + "\n")
        for input_id, input_seq, replicate_ids in zip(unique_representative_ids, unique_input_seqs, unique_member_ids):

            # Profile each unique input sequence.
            input_queue.put((input_id, input_seq))
            self.profiled_unique_seq_count += 1
            replicate_dict[input_id] = replicate_ids

            if self.profiled_unique_seq_count != write_points[write_point_index]:
                # Keep adding sequences to the profiling queue until the write point is hit.
                continue

            write_results()
            replicate_dict = {}

            self.progress.update("%d of %d unique sequences have been profiled"
                                 % (self.profiled_unique_seq_count, len(unique_input_seqs)))

            write_point_index += 1
            if write_point_index == len(write_points):
                break

        for process in processes:
            process.terminate()
        nontrna_report.close()

        set_meta_values()

        output_run_info()

        trnaseq_db.disconnect()

        self.progress.end()


    def normalize_ends(self):
        self.progress.new("Normalizing the 3' and 5' ends of sequences")
        self.progress.update("...")

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)

        # Split each sequence into the 3' acceptor sequence and the rest of the sequence.
        normalized_seqs, num_extra_fiveprime_list, fiveprime_seqs, threeprime_seqs = zip(*
            [(seq[num_extra_fiveprime: acceptor_start], num_extra_fiveprime, seq[: num_extra_fiveprime], seq[acceptor_start: ])
             for seq, num_extra_fiveprime, acceptor_start
             in zip(trnaseq_db.db.get_single_column_from_table('sequences', 'sequence'),
                    trnaseq_db.db.get_single_column_from_table('basic_info', 'num_extra_fiveprime'),
                    trnaseq_db.db.get_single_column_from_table('features', 'acceptor_start'))])

        extra_dicts = []
        for num_extra_fiveprime, fiveprime_seq, threeprime_seq, input_seq_count in zip(
            num_extra_fiveprime_list, fiveprime_seqs, threeprime_seqs, trnaseq_db.db.get_single_column_from_table('sequences', 'replicate_count')):
            extra_dict = {}
            extra_dict['num_extra_fiveprime'] = num_extra_fiveprime
            extra_dict['fiveprime_seq'] = fiveprime_seq
            extra_dict['threeprime_seq'] = threeprime_seq
            extra_dict['input_seq_count'] = input_seq_count
            extra_dicts.append(extra_dict)

        clusters = SequenceDereplicator(
            trnaseq_db.db.get_single_column_from_table('sequences', 'name'),
            normalized_seqs,
            extra_list=extra_dicts).full_length_dereplicate()

        new_clusters = []
        for cluster in clusters:
            new_member_seq_dicts = []
            cluster_input_seq_count = 0
            for member_id, member_seq_dict in zip(cluster[1], cluster[2]):
                member_seq_dict['member_id'] = member_id
                new_member_seq_dicts.append(member_seq_dict)
                cluster_input_seq_count += member_seq_dict['input_seq_count']
            new_clusters.append({'seq': cluster[0],
                                 'input_seq_count': cluster_input_seq_count,
                                 'member_seq_dicts': new_member_seq_dicts})
        clusters = new_clusters
        clusters.sort(key=lambda d: -d['input_seq_count'])

        # Sort member IDs in descending order of number of extra 5' bases and then descending order of input sequence count.
        for cluster in clusters:
            cluster['member_seq_dicts'].sort(key=lambda d: (-d['num_extra_fiveprime'], -d['input_seq_count']))
        # The representative member ID of each cluster is chosen as follows, based on sequences present in the cluster:
        # 1. Most abundant full-length tRNA (no extra 5' bases), ignoring acceptor sequence
        # 2. Most abundant longer-than-full-length tRNA
        # 3. Most abundant fragmentary tRNA
        representative_ids = []
        for cluster in clusters:
            member_seq_dicts = cluster['member_seq_dicts']
            if member_seq_dicts[0]['num_extra_fiveprime'] > 0: # check for extra 5' bases
                if member_seq_dicts[-1]['num_extra_fiveprime'] == 0:
                    # If there is also a sequence without extra 5' bases in the cluster, it must be a full-length sequence.
                    representative_id = member_seq_dicts[-1]['member_id']
                else:
                    representative_id = member_seq_dicts[0]['member_id']
            else:
                # All sequences in the cluster are full-length or a fragment.
                representative_id = member_seq_dicts[-1]['member_id']
            cluster['representative_id'] = representative_id

        if self.write_normalization_reports:
            fiveprime_normalization_report = open(self.fiveprime_normalization_report_path, 'w')
            fiveprime_normalization_report.write("\t".join(self.FIVEPRIME_NORMALIZATION_REPORT_HEADER) + "\n")
            threeprime_normalization_report = open(self.threeprime_normalization_report_path, 'w')
            threeprime_normalization_report.write("\t".join(self.THREEPRIME_NORMALIZATION_REPORT_HEADER) + "\n")
        normalization_table_entries = []
        for cluster in clusters:
            representative_id = cluster['representative_id']
            member_seq_dicts = cluster['member_seq_dicts']
            fiveprime_unique_count = 0
            fiveprime_input_count = 0
            threeprime_variant_count_dict = OrderedDict(
                [(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
            for member_seq_dict in member_seq_dicts:
                member_input_seq_count = member_seq_dict['input_seq_count']
                if member_seq_dict['fiveprime_seq']:
                    # Normalization trimmed the 5' end of the sequence.
                    fiveprime_unique_count += 1
                    fiveprime_input_count += member_input_seq_count
                    if self.write_normalization_reports:
                        fiveprime_normalization_report.write(
                            representative_id + "\t"
                            + member_seq_dict['member_id'] + "\t"
                            + member_seq_dict['fiveprime_seq'] + "\t"
                            + member_seq_dict['threeprime_seq'] + "\t"
                            + str(member_seq_dict['input_seq_count']) + "\n")
                elif self.write_normalization_reports:
                    threeprime_normalization_report.write(
                        representative_id + "\t"
                        + member_seq_dict['member_id'] + "\t"
                        + member_seq_dict['threeprime_seq'] + "\t"
                        + str(member_seq_dict['input_seq_count']) + "\n")
                threeprime_variant_count_dict[member_seq_dict['threeprime_seq']] += member_input_seq_count

            normalization_table_entries.append(
                (representative_id, len(member_seq_dicts), cluster['input_seq_count'], cluster['seq'], fiveprime_unique_count, fiveprime_input_count)
                + tuple([v for v in threeprime_variant_count_dict.values()]))

        trnaseq_db.db._exec_many(
            '''INSERT INTO %s VALUES (%s)'''
            % ('normalization', ','.join('?' * len(tables.trnaseq_normalization_table_structure))),
            normalization_table_entries)

        normalized_seq_count = len(clusters)
        trnaseq_db.db.set_meta_value('num_unique_normalized_trna_seqs', normalized_seq_count)
        print()
        self.run.info("Unique tRNA sequences, normalizing the ends", normalized_seq_count)

        trnaseq_db.disconnect()

        self.progress.end()


    def dereplicate_subsequences(self):

        def dereplicate_from_threeprime_end():
            # Track certain information associated with each normalized sequence.
            extra_dicts = []
            for t in trnaseq_db.db.get_some_columns_from_table(
                'normalization',
                'input_seq_count,fiveprime_input_seq_count,'
                + ','.join([threeprime_variant + '_input_seq_count' for threeprime_variant in THREEPRIME_VARIANTS])):
                extra_dict = {}
                extra_dict['end_input_seq_count'] = t[0]
                extra_dict['fiveprime_input_seq_count'] = t[1]
                for i, threeprime_variant in enumerate(THREEPRIME_VARIANTS, 2):
                    extra_dict[threeprime_variant + '_input_seq_count'] = t[i]
                extra_dicts.append(extra_dict)

            clusters, memberships = SequenceDereplicator(
                trnaseq_db.db.get_single_column_from_table('normalization', 'name'),
                [seq[::-1] for seq in trnaseq_db.db.get_single_column_from_table('normalization', 'sequence')], # change seq orientation to 3'-5' to dereplicate from 3' end
                extra_list=extra_dicts).prefix_dereplicate()

            # Here is the cluster format as it is returned by prefix_dereplicate:
            # clusters = [(seed seq A ID, seed sequence A, [(member seq X ID, member seq X length, member seq X extra dict), ...]), ...]
            # Change the cluster format to this:
            # clusters = [{'representative_id': seed seq A ID, 'seq': seed sequence A, 'member_seq_dicts': [{'member_id': member seq X ID, 'start_position': ...}, {...}, ...]}, ...]
            new_clusters = []
            for cluster in clusters:
                new_cluster = {}
                new_cluster['representative_id'] = cluster[0]
                new_cluster['seq'] = seq = cluster[1][::-1] # revert seq orientation to 5'-3'
                reference_length = len(seq)
                new_member_seq_dicts = []
                for member_seq_info in cluster[2]:
                    new_member_seq_dict = {}
                    new_member_seq_dict['member_id'] = member_seq_info[0]
                    new_member_seq_dict['start_position'] = reference_length - member_seq_info[1]
                    new_member_seq_dict['end_position'] = reference_length
                    new_member_seq_dict['interior_input_seq_count'] = 0
                    new_member_seq_dict.update(member_seq_info[2])
                    new_member_seq_dicts.append(new_member_seq_dict)
                new_cluster['member_seq_dicts'] = new_member_seq_dicts
                new_clusters.append(new_cluster)
            clusters = new_clusters

            # Here is the membership format as it is returned by prefix_dereplicate:
            # memberships = [(member sequence X, [seed seq A ID, seed seq B ID, ...], member seq X extra dict), ...]
            # Change the membership format to this:
            # memberships = [{'member_id': member sequence X, 'representative_ids': [seed seq A ID, seed seq B ID, ...], 'end_input_seq_count': member seq X end input seq count, ...}, ...]
            new_memberships = []
            for membership in memberships:
                new_membership_dict = {}
                new_membership_dict['member_id'] = membership[0]
                new_membership_dict['representative_ids'] = membership[1]
                new_membership_dict.update(membership[2])
                new_memberships.append(new_membership_dict)
            memberships = new_memberships

            return clusters, memberships


        ####################################
        # dereplicate_subsequences main body
        ####################################

        self.progress.new("Dereplicating normalized tRNA sequences")

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)

        self.progress.update("Dereplicating normalized tRNA sequences from the 3' end")

        clusters, memberships = dereplicate_from_threeprime_end()

        # "Multiplicity" of a clustered normalized tRNA sequence
        # = number of clusters that the sequence is in
        #   * number of input sequences represented by the normalized sequence
        member_multiplicity_dict = {}
        for membership in memberships:
            member_multiplicity_dict[membership['member_id']] = len(membership['representative_ids']) * membership['end_input_seq_count']

        if self.write_normalization_reports:
            subsequence_report_dict = OrderedDict()
            for cluster in clusters:
                end_subsequence_ids = []
                subsequence_report_dict[cluster['representative_id']] = {
                    'end_subsequence_ids': end_subsequence_ids,
                    'interior_subsequence_ids': []}
                for member_seq_dict in cluster['member_seq_dicts']:
                    end_subsequence_ids.append(member_seq_dict['member_id'])

        print()
        self.run.info("Number of normalized sequences after dereplicating from the 3' end", len(clusters))


        # The following dereplication steps rely on mapping by Bowtie2.
        os.mkdir(self.derep_mapping_temp_dir)

        self.progress.update("Dereplicating short unprofiled tRNA fragments from the 3' end")

        # Use normalized tRNA sequences as alignment targets.
        # Note: targets could also include sequences with extra 5' bases.
        # This would require relating the these sequences to their normalized sequence cluster (representative ID).
        # This would also require storing the extra 5' sequence of the aligned query
        # and adding member information, including the remainder of the query sequence, to the cluster.
        os.mkdir(self.derep_interior_temp_dir)
        target_fasta_path = os.path.join(self.derep_interior_temp_dir, "target.fasta")
        target_fasta = fastalib.FastaOutput(target_fasta_path)
        # Only consider 3' variants that are subsequence-unique.
        target_threeprime_additions = ['CCC', 'CCG', 'CCT',
                                       'CAA', 'CGA', 'CTA',
                                       'ACA', 'GCA', 'TCA',
                                       'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT',
                                       'CCACA', 'CCACC', 'CCACG', 'CCACT',
                                       'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT',
                                       'CCATA', 'CCATC', 'CCATG', 'CCATT']
        for cluster in clusters:
            for threeprime_variant in target_threeprime_additions:
                target_fasta.write_id(cluster['representative_id'] + "_" + threeprime_variant)
                target_fasta.write_seq(cluster['seq'] + threeprime_variant)
        target_fasta.close()

        # Write a FASTA file of unique sequences that were not successfully profiled as tRNA.
        # Store the input sequence count of each unique sequence.
        nontrna_fasta_path = os.path.join(self.derep_interior_temp_dir, "nontrna.fasta")
        nontrna_fasta = fastalib.FastaOutput(nontrna_fasta_path)
        interior_input_seq_count_dict = {}
        nontrna_report = open(self.nontrna_report_path) # a tab-delimited file with four columns
        nontrna_report.readline() # ignore the header
        for line in nontrna_report:
            representative_id, _, interior_input_seq_count, seq = line.rstrip().split("\t")
            interior_input_seq_count_dict[representative_id] = interior_input_seq_count
            nontrna_fasta.write_id(representative_id)
            nontrna_fasta.write_seq(seq)
        nontrna_report.close()
        nontrna_fasta.close()

        # Build a bowtie2 index from the target FASTA file.
        target_index_path = os.path.join(self.derep_interior_temp_dir, "normalized_trna")
        cmd_line = "bowtie2-build -f %s %s" % (target_fasta_path, target_index_path)
        utils.run_command(cmd_line, self.log_path, remove_log_file_if_exists=False)

        # Map "non-tRNA" queries to normalized sequence targets.
        # Note: --score-min C,0,0 ensures that end-to-end alignments have no mismatches or gaps.
        # An end-to-end alignment with mismatches or gaps has a negative score.
        sam_path = os.path.join(self.derep_interior_temp_dir, "nontrna.sam")
        cmd_line = "bowtie2 -x %s -f %s -S %s --score-min C,0,0 -k %d -p %d" % (
            target_index_path, nontrna_fasta_path, sam_path, len(clusters), self.num_threads)
        # REMOVE
        # cmd_line = "bowtie2 -x %s -f %s -S %s -k 1 -p %d" % (
        #     target_index_path, nontrna_fasta_path, sam_path, self.num_threads)
        utils.run_command(cmd_line, self.log_path, remove_log_file_if_exists=False)

        # Convert SAM to BAM.
        raw_bam_path = os.path.join(self.derep_interior_temp_dir, "nontrna-RAW.bam")
        pysam.view('-bS', sam_path, '-F', '4', '-o', raw_bam_path, catch_stdout=False)

        # Sort and index BAM.
        bam_path = os.path.join(self.derep_interior_temp_dir, "nontrna.bam")
        cmd_line = "anvi-init-bam %s -o %s -T %d" % (raw_bam_path, bam_path, self.num_threads)
        utils.run_command(cmd_line, self.log_path, remove_log_file_if_exists=False)

        cluster_dict = {cluster['representative_id']: cluster for cluster in clusters}
        # Add a cluster entry for the number of unique interior sequences.
        for cluster in clusters:
            cluster['interior_unique_count'] = 0

        # Parse the alignments, which are grouped by reference sequence.
        bam_file = pysam.AlignmentFile(bam_path, 'rb')
        interior_reference_count_dict = {}
        unique_interior_subsequence_count = 0
        target_reference_id = ''
        tiny_threeprime_trna_ids = []
        for aligned_segment in bam_file.fetch():
            member_id = aligned_segment.query_name
            target_reference_start = aligned_segment.reference_start
            target_reference_end = aligned_segment.reference_end # pythonic stop position

            if aligned_segment.reference_name == target_reference_id:
                if target_reference_start > true_seq_length or target_reference_end > true_seq_length:
                    tiny_threeprime_trna_ids.append(member_id)
                    continue
            else:
                target_reference_id = aligned_segment.reference_name
                added_threeprime_seq = target_reference_id.split('_')[-1]
                true_reference_id = target_reference_id[: -len(added_threeprime_seq) - 1]
                cluster = cluster_dict[true_reference_id]
                true_seq = cluster['seq']
                true_seq_length = len(true_seq)
                member_seq_dicts = cluster['member_seq_dicts']

                if target_reference_start > true_seq_length or target_reference_end > true_seq_length:
                    # The query matches the 3' end of the reference.
                    # This implies that the query is a short tRNA fragment that could not be successfully profiled.
                    # Ignore these sequences, but write their IDs to a text file.
                    tiny_threeprime_trna_ids.append(member_id)
                    continue

                try:
                    interior_reference_count_dict[true_reference_id] += 1
                except KeyError:
                    interior_reference_count_dict[true_reference_id] = 0

            print(true_reference_id)
            cluster['interior_unique_count'] += 1 # count the unique subsequence
            member_seq_dict = {} # add information on the interior subsequence to the cluster
            member_seq_dict['member_id'] = member_id
            if self.write_normalization_reports:
                subsequence_report_dict[true_reference_id]['interior_subsequence_ids'].append(member_id)
            member_seq_dict['end_position'] = target_reference_end
            member_seq_dict['end_input_seq_count'] = 0
            member_seq_dict['fiveprime_input_seq_count'] = 0
            for threeprime_variant in THREEPRIME_VARIANTS:
                member_seq_dict[threeprime_variant + '_input_seq_count'] = 0
            interior_input_seq_count = interior_input_seq_count_dict[member_id]
            member_seq_dict['interior_input_seq_count'] = interior_input_seq_count
            try:
                member_multiplicity_dict[member_id] += interior_input_seq_count
            except KeyError:
                member_multiplicity_dict[member_id] = interior_input_seq_count
                unique_interior_subsequence_count += 1

        with open(self.tiny_threeprime_trna_ids_path, 'w') as tiny_threeprime_trna_ids_file:
            tiny_threeprime_trna_ids_file.write("unique_seq_id\n")
            for tiny_threeprime_trna_id in set(tiny_threeprime_trna_ids):
                tiny_threeprime_trna_ids_file.write(tiny_threeprime_trna_id + "\n")

        import sys
        sys.exit(0)

        subseq_table_entries = []
        threeprime_variant_keys = [threeprime_variant + 'input_seq_count' for threeprime_variant in THREEPRIME_VARIANTS]
        for cluster in clusters:
            member_seq_dicts = cluster[2]
            end_input_seq_count = 0
            interior_input_seq_count = 0
            fiveprime_input_seq_count = 0
            member_multiplicity_sum = 0
            threeprime_variant_input_seq_counts = [0 for _ in threeprime_variant_keys]
            for member_seq_dict in member_seq_dicts:
                end_input_seq_count += member_seq_dict['end_input_seq_count']
                interior_input_seq_count += member_seq_dict['interior_input_seq_count']
                fiveprime_input_seq_count += member_seq_dict['fiveprime_input_seq_count']
                member_multiplicity_sum += member_multiplicity_dict[member_seq_dict['member_id']]
                for i, threeprime_variant_key in enumerate(threeprime_variant_keys):
                    threeprime_variant_input_seq_counts[i] += member_seq_dict[threeprime_variant_key]
            input_seq_count = end_input_seq_count + interior_input_seq_count
            average_multiplicity = round(member_multiplicity_sum / input_seq_count, 1)
            subseq_table_entries.append((
                cluster[0], # representative ID
                len(member_seq_dicts), # number of unique dereplicated sequences
                input_seq_count,
                average_multiplicity,
                cluster[-1], # count of unique interior sequences
                interior_input_seq_count,
                fiveprime_input_seq_count)
                + tuple(threeprime_variant_input_seq_counts))

        trnaseq_db.disconnect()

        if self.write_normalization_reports:
            with open(self.subsequence_report_path, 'w') as subsequence_report:
                subsequence_report.write("\t".join(self.SUBSEQUENCE_REPORT_HEADER) + "\n")
                for representative_id, subsequence_id_lists in subsequence_report_dict.items():
                    subsequence_report.write(representative_id + "\t"
                                             + ",".join(subsequence_id_lists[0]) + "\t"
                                             + ",".join(subsequence_id_lists[1]) + "\n")

        self.run.info("Number of 3' dereplicated sequences with interior subsequences", sum(interior_reference_count_dict))
        self.run.info("Number of unique interior subsequences found", unique_interior_subsequence_count)

        self.progress.end()


    def process(self):
        self.sanity_check()

        self.generate_trnaseq_db()

        # Profile tRNA features in sequences.
        # tRNA results are stored in the database, while non-tRNA is written (optionally) to a file.
        self.profile_trna()

        self.normalize_ends()
        self.dereplicate_subsequences()