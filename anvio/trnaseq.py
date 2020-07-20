# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes for tRNA-seq dataset operations. anvi-trnaseq is the default client using this module. """

import itertools
import multiprocessing
import os
import pysam
import shutil

from collections import OrderedDict

import anvio
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.trnaidentifier as trnaidentifier
import anvio.utils as utils

from anvio.agglomeration import Agglomerator
from anvio.constants_package.trnaseq import THREEPRIME_VARIANTS, TRNA_FEATURE_NAMES
from anvio.dbops_package.trnaseq import TRNASeqDatabase
from anvio.drivers.bowtie2 import Bowtie2
from anvio.errors import ConfigError, FilesNPathsError
from anvio.sequence import SequenceDereplicator

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"



class UniqueSeq:
    def __init__(self, seq, input_ids, identification_method=None, acceptor_length=None, extra_fiveprime_length=None, skip_init=False):
        self.string = seq
        self.input_ids = input_ids
        self.identification_method = identification_method # 'profiled' or 'mapped' if tRNA
        self.acceptor_length = acceptor_length
        self.extra_fiveprime_length = extra_fiveprime_length

        if skip_init:
            self.input_count = None
            self.representative_id = None
        else:
            self.init()


    def init(self):
        self.input_count = len(self.input_ids)
        self.representative_id = self.input_ids[0]


class TrimmedSeq:
    def __init__(self, seq, unique_seqs, skip_init=False):
        self.string = seq
        self.unique_seqs = unique_seqs # list of UniqueSeq objects

        if skip_init:
            self.input_count = None
            self.unique_with_extra_fiveprime_count = None
            self.input_with_extra_fiveprime_count = None
            self.representative_id = None
            self.input_acceptor_variant_count_dict = None
            self.identification_method = None
        else:
            self.init()


    def init(self):
        self.input_count = sum([len(unique_seq.input_ids) for unique_seq in self.unique_seqs])

        self.unique_with_extra_fiveprime_count = sum([1 if unique_seq.extra_fiveprime_length else 0
                                                      for unique_seq in self.unique_seqs])
        self.input_with_extra_fiveprime_count = sum([len(unique_seq.input_ids) if unique_seq.extra_fiveprime_length else 0
                                                     for unique_seq in self.unique_seqs])


        # The representative ID is chosen as follows:
        # 1. Most abundant full-length tRNA (no extra 5' bases), ignoring acceptor sequence
        # 2. Most abundant longer-than-full-length tRNA
        # 3. Most abundant fragmentary tRNA
        unique_seqs = sorted(self.unique_seqs, key=lambda unique_seq: (-unique_seq.extra_fiveprime_length, -len(unique_seq.input_ids)))

        if unique_seqs[0].extra_fiveprime_length > 0:
            # If there is also a unique sequence that was ultimately trimmed down
            # to the same sequence as the sequence with extra 5' bases, it must be a full-length sequence.
            if unique_seqs[-1].extra_fiveprime_length == 0:
                representative_id = unique_seqs[-1].representative_id
            else:
                representative_id = unique_seqs[0].representative_id
        else:
            # ALL unique sequences are EITHER full-length OR a fragment.
            representative_id = unique_seqs[0].representative_id

        self.representative_id = representative_id

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for unique_seq in self.unique_seqs:
            if unique_seq.acceptor_length: # unique_seq need not have an acceptor
                acceptor_string = unique_seq.string[-unique_seq.acceptor_length: ]
                input_acceptor_variant_count_dict[acceptor_string] += unique_seq.input_count
        self.input_acceptor_variant_count_dict = input_acceptor_variant_count_dict

        identification_methods = set(unique_seq.identification_method for unique_seq in self.unique_seqs)
        if len(identification_methods) == 1:
            self.identification_method = identification_methods.pop()
        else:
            raise ConfigError("A TrimmedSeq should not be made from UniqueSeq objects "
                              "with different identification methods. "
                              "Trimmed tRNA sequences will EITHER be formed from "
                              "\"profiled\" tRNA sequences or \"mapped\" tRNA sequences, "
                              "because they are of different lengths and are fragments from different parts of the tRNA.")


class NormalizedSeq:
    def __init__(self, trimmed_seqs, start_positions=None, end_positions=None, skip_init=False):
        self.trimmed_seqs = trimmed_seqs # list of TrimmedSeq objects
        self.representative_id = trimmed_seqs[0].representative_id
        self.string = trimmed_seqs[0].string
        if start_positions and end_positions:
            self.start_positions = start_positions
            self.end_positions = end_positions
        elif (not start_positions) and (not end_positions):
            # Trimmed seqs were dereplicated from the 3' end of the normalized sequence.
            normalized_seq_length = len(self.string)
            self.start_positions = [normalized_seq_length - len(trimmed_seq.string) for trimmed_seq in self.trimmed_seqs]
            self.end_positions = [normalized_seq_length] * len(trimmed_seqs)
        else:
            self.start_positions = None
            self.end_positions = None

        if skip_init:
            self.input_count = None
            self.input_with_extra_fiveprime_count = None
            self.input_acceptor_variant_count_dict = None
            self.count_of_trimmed_seqs_mapped_to_threeprime_end = None
            self.count_of_input_seqs_mapped_to_threeprime_end = None
            self.count_of_trimmed_seqs_mapped_to_interior = None
            self.count_of_input_seqs_mapped_to_interior = None
            self.count_of_trimmed_seqs_mapped_to_fiveprime_end = None
            self.count_of_input_seqs_mapped_to_fiveprime_end = None
        else:
            self.init()


    def init(self):
        self.input_count = sum([trimmed_seq.input_count for trimmed_seq in self.trimmed_seqs])

        self.input_with_extra_fiveprime_count = sum([trimmed_seq.input_with_extra_fiveprime_count for trimmed_seq in self.trimmed_seqs])

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        input_acceptor_variant_counts = []
        for trimmed_seq in self.trimmed_seqs:
            for acceptor_string, input_count in trimmed_seq.input_acceptor_variant_count_dict.items():
                if input_count > 0:
                    input_acceptor_variant_count_dict[acceptor_string] += input_count
        self.input_acceptor_variant_count_dict = input_acceptor_variant_count_dict

        count_of_trimmed_seqs_mapped_to_threeprime_end = 0
        count_of_input_seqs_mapped_to_threeprime_end = 0
        count_of_trimmed_seqs_mapped_to_interior = 0
        count_of_input_seqs_mapped_to_interior = 0
        count_of_trimmed_seqs_mapped_to_fiveprime_end = 0
        count_of_input_seqs_mapped_to_fiveprime_end = 0
        for trimmed_seq, start_position, end_position in zip(self.trimmed_seqs, self.start_positions, self.end_positions):
            if trimmed_seq.identification_method == 'mapped':
                # TrimmedSeqs are comprised of EITHER "profiled" OR "mapped" UniqueSeqs.
                # Mapped short 3' sequences (<30 nt) can form TrimmedSeqs from multiple UniqueSeqs
                # that only differ in the acceptor sequence variant.
                if end_position == len(self.string):
                    count_of_trimmed_seqs_mapped_to_threeprime_end += 1
                    count_of_input_seqs_mapped_to_threeprime_end += trimmed_seq.input_count
                elif trimmed_seq.unique_with_extra_fiveprime_count > 0:
                    count_of_trimmed_seqs_mapped_to_fiveprime_end += 1
                    count_of_input_seqs_mapped_to_fiveprime_end += trimmed_seq.input_count
                else:
                    count_of_trimmed_seqs_mapped_to_interior += 1
                    count_of_input_seqs_mapped_to_interior += trimmed_seq.input_count
        self.count_of_trimmed_seqs_mapped_to_threeprime_end = count_of_trimmed_seqs_mapped_to_threeprime_end
        self.count_of_input_seqs_mapped_to_threeprime_end = count_of_input_seqs_mapped_to_threeprime_end
        self.count_of_trimmed_seqs_mapped_to_interior = count_of_trimmed_seqs_mapped_to_interior
        self.count_of_input_seqs_mapped_to_interior = count_of_input_seqs_mapped_to_interior
        self.count_of_trimmed_seqs_mapped_to_fiveprime_end = count_of_trimmed_seqs_mapped_to_fiveprime_end
        self.count_of_input_seqs_mapped_to_fiveprime_end = count_of_input_seqs_mapped_to_fiveprime_end


class AgglomeratedSeq:
    def __init__(self, normalized_seqs):
        self.normalized_seqs = normalized_seqs


class TRNASeqDataset:
    UNIQUED_NONTRNA_HEADER = ["representative_id", "input_ids", "input_count", "sequence"]
    UNIQUED_TRNA_HEADER = ["representative_id", "input_ids"]
    TRIMMED_ENDS_HEADER = ["representative_id", "unique_id", "fiveprime_sequence", "threeprime_sequence", "input_seq_count"]
    NORMALIZED_FRAGMENTS_HEADER = ["representative_id", "trimmed_id", "start", "end"]
    AGGLOMERATED_SEQS_HEADER = ["representative_id", "normalized_ids"]

    SHORT_FRAGMENT_LENGTH = 30 # number of nucleotides from 3'-CCANN through T arm
    TARGET_ACCEPTOR_VARIANTS = ['CCC', 'CCG', 'CCT', # used in mapping fragments to tRNA interior/5' end
                                'CAA', 'CGA', 'CTA',
                                'ACA', 'GCA', 'TCA',
                                'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT',
                                'CCACA', 'CCACC', 'CCACG', 'CCACT',
                                'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT',
                                'CCATA', 'CCATC', 'CCATG', 'CCATT']

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_fasta_path = A('fasta_file')
        self.project_name = A('project_name')
        self.output_dir = os.path.abspath(A('output_dir')) if A('output_dir') else None
        self.description_file_path = os.path.abspath(A('description')) if A('description') else None
        self.skip_fasta_check = A('skip_fasta_check')
        self.num_threads = A('num_threads')
        self.write_buffer_size = A('write_buffer_size')
        self.max_possible_alignments = A('max_possible_alignments')

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using --fasta-file or -f.")
        if not self.project_name:
            raise ConfigError("Please set a project name using --project-name or -n.")
        if not self.output_dir:
            raise ConfigError("Please provide a non-existent output directory using --output-dir or -o.")

        self.input_seq_count = 0

        self.trnaseq_db_path = os.path.join(self.output_dir, self.project_name + "-TRNASEQ.db")

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
        # The following file paths are used for all Bowtie2 query and target files.
        self.query_fasta_path = os.path.join(self.derep_mapping_temp_dir, "query.fasta")
        self.target_fasta_path = os.path.join(self.derep_mapping_temp_dir, "target.fasta")

        self.uniqued_nontrna_path = os.path.join(self.output_dir, self.project_name + "-UNIQUED_NONTRNA.txt")
        self.uniqued_trna_path = os.path.join(self.output_dir, self.project_name + "-UNIQUED_TRNA.txt")
        self.trimmed_ends_path = os.path.join(self.output_dir, self.project_name + "-TRIMMED_ENDS.txt")
        self.normalized_fragments_path = os.path.join(self.output_dir, self.project_name + "-NORMALIZED_FRAGMENTS.txt")
        self.agglomerated_seqs_path = os.path.join(self.output_dir, self.project_name + "-AGGLOMERATED_SEQS.txt")

        self.log_path = os.path.join(self.output_dir, "log.txt")

        self.unique_nontrna_seqs = []
        self.unique_trna_seqs = []
        self.trimmed_trna_seqs = []
        self.normalized_trna_seqs = []
        self.agglomerated_trna_seqs = []

        self.counts_of_normalized_seqs_containing_trimmed_seqs = [] # same length as self.trimmed_trna_seqs
        # "Multiplicity" of a trimmed tRNA sequence
        # = number of normalized sequences that the sequence is in * number of input sequences represented by the trimmed sequence
        self.multiplicities_of_trimmed_seqs_among_normalized_seqs = [] # same length as self.trimmed_trna_seqs
        self.average_multiplicities_of_normalized_seqs = [] # same length as self.normalized_trna_seqs


    def check_programs(self):
        utils.is_program_exists('bowtie2')


    def check_output_dir(self):
        if os.path.exists(self.output_dir):
            raise ConfigError("The directory that was specified by --output-dir or -o, %s, already exists. "
                              "Please try again with a non-existent directory." % self.output_dir)
        os.mkdir(self.output_dir)
        filesnpaths.is_output_dir_writable(self.output_dir)


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
            self.description = self.description_file_path.read()
        else:
            self.description = None
        self.run.info("Description", self.description_file_path if self.description_file_path else "No description given")

        self.run.log_path = self.log_path

        if self.max_possible_alignments < 2:
            raise ConfigError("The maximum number of possible alignments in Bowtie2 clustering "
                              "of normalized sequences to themselves cannot be less than 2, "
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


    def get_unique_input_seqs(self):
        fasta = fastalib.SequenceSource(self.input_fasta_path)
        ids = []
        seqs = []
        input_seq_count = 0
        while next(fasta):
            ids.append(fasta.id)
            seqs.append(fasta.seq)
            input_seq_count += 1
        fasta.close()
        self.input_seq_count = input_seq_count

        clusters = SequenceDereplicator(ids, seqs).full_length_dereplicate()

        unique_input_seqs = [UniqueSeq(t[0], t[1]) for t in clusters]

        return unique_input_seqs


    def generate_trnaseq_db(self):
        meta_values = {'project_name': self.project_name,
                       'description': self.description if self.description else '_No description is provided_',
                       'max_possible_alignments': self.max_possible_alignments}
        TRNASeqDatabase(self.trnaseq_db_path, quiet=False).create(meta_values)


    def get_write_points(self, item_count):
        # Determine the points when to write processed items.
        # Write points are cumulative numbers of items processed.
        num_chunks = item_count // self.write_buffer_size
        if num_chunks == 0:
            num_chunks += 1
            write_points = [item_count]
        else:
            write_points = [self.write_buffer_size * (i + 1) for i in range(num_chunks)]
            if item_count % num_chunks > 0:
                num_chunks += 1
                write_points.append(item_count)
        return write_points


    def write_results(self, output_queue, chunk_dict, trnaseq_db):

        # List of entries for each table
        trnaseq_sequences_table_entries = []
        trnaseq_info_table_entries = []
        trnaseq_features_table_entries = []
        trnaseq_unconserved_table_entries = []
        trnaseq_unpaired_table_entries = []

        retrieved_profile_count = 0
        unique_trna_seqs = []
        unique_nontrna_seqs = []

        while retrieved_profile_count < self.profiled_unique_seq_count:
            # Retrieve profiles from the output queue and write tRNA profiles to the database.
            trna_profile = output_queue.get()
            retrieved_profile_count += 1

            output_id = trna_profile.name
            output_seq = trna_profile.input_seq

            if not trna_profile.is_predicted_trna:
                unique_nontrna_seqs.append(chunk_dict[output_id])
                continue

            unique_seq = chunk_dict[output_id]
            unique_seq.identification_method = 'profiled'
            unique_seq.acceptor_length = len(trna_profile.acceptor_variant_string)
            unique_seq.extra_fiveprime_length = trna_profile.num_extra_fiveprime
            unique_trna_seqs.append(unique_seq)

            self.unique_trna_count += 1
            output_seq_length = len(output_seq)

            num_replicates = len(unique_seq.input_ids)
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
            elif trna_profile.acceptor_variant_string == 'CCA':
                self.trna_with_threeprime_cca_count += num_replicates
            elif trna_profile.acceptor_variant_string == 'CC':
                self.trna_with_threeprime_cc_count += num_replicates
            elif trna_profile.acceptor_variant_string == 'C':
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

            trnaseq_info_table_entries.append(
                (output_id,
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

        return unique_trna_seqs, unique_nontrna_seqs


    def set_meta_values(self, trnaseq_db):
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


    def output_run_info(self):
        self.run.info("Reads processed", self.input_seq_count)
        self.run.info("Reads profiled as tRNA", self.trna_count)
        self.run.info("Unique profiled tRNA sequences", self.unique_trna_count)
        self.run.info("Profiled reads with anticodon", self.trna_containing_anticodon_count)
        self.run.info("Profiled reads spanning acceptor stem", self.mature_trna_count)
        self.run.info("Profiled reads with 1-5 extra 5' bases", self.trna_with_one_to_five_extra_fiveprime_bases_count)
        self.run.info("Profiled reads with >5 extra 5' bases", self.trna_with_more_than_five_extra_fiveprime_bases_count)
        self.run.info("Profiled reads with extrapolated 5' feature", self.trna_with_extrapolated_fiveprime_feature_count)
        self.run.info("Profiled reads ending in 3'-CCA", self.trna_with_threeprime_cca_count)
        self.run.info("Profiled reads ending in 3'-CC", self.trna_with_threeprime_cc_count)
        self.run.info("Profiled reads ending in 3'-C", self.trna_with_threeprime_c_count)
        self.run.info("Profiled reads ending in 3'-NCA/CNA/CCN", self.trna_with_threeprime_nca_cna_ccn_count)
        self.run.info("Profiled reads ending in 3'-CCAN/CCANN", self.trna_with_threeprime_ccan_ccann_count)


    def profile_trna(self, unique_input_seqs):
        self.progress.new("Profiling input sequences for tRNA features")

        write_points = self.get_write_points(len(unique_input_seqs))
        write_point_index = 0

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        processes = [multiprocessing.Process(target=trnaidentifier.profile_wrapper, args=(input_queue, output_queue))
                     for _ in range(self.num_threads)]
        for process in processes:
            process.start()


        chunk_dict = {}
        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        for unique_seq in unique_input_seqs:
            # Profile each unique input sequence.
            input_queue.put((unique_seq.representative_id, unique_seq.string))
            self.profiled_unique_seq_count += 1
            chunk_dict[unique_seq.representative_id] = unique_seq

            if self.profiled_unique_seq_count != write_points[write_point_index]:
                # Keep adding sequences to the profiling queue until the write point is hit.
                continue

            chunk_unique_trna_seqs, chunk_unique_nontrna_seqs = self.write_results(output_queue, chunk_dict, trnaseq_db)
            self.unique_trna_seqs.extend(chunk_unique_trna_seqs)
            self.unique_nontrna_seqs.extend(chunk_unique_nontrna_seqs)
            chunk_dict = {}

            self.progress.update("%d of %d unique sequences have been profiled"
                                 % (self.profiled_unique_seq_count, len(unique_input_seqs)))

            write_point_index += 1
            if write_point_index == len(write_points):
                break

        for process in processes:
            process.terminate()

        self.set_meta_values(trnaseq_db)

        trnaseq_db.disconnect()

        self.progress.end()

        self.output_run_info()


    def trim_ends(self, unique_trna_seqs):
        self.progress.new("Normalizing the 3' and 5' ends of sequences")

        self.trimmed_trna_seqs.extend([TrimmedSeq(seq, clustered_unique_seqs) for seq, _, clustered_unique_seqs in
                                       SequenceDereplicator(
                                           [unique_seq.representative_id for unique_seq in unique_trna_seqs],
                                           [unique_seq.string[unique_seq.extra_fiveprime_length: -unique_seq.acceptor_length]
                                            for unique_seq in unique_trna_seqs],
                                           extra_list=unique_trna_seqs).full_length_dereplicate()])
        self.trimmed_trna_seqs.sort(key=lambda trimmed_seq: -trimmed_seq.input_count)

        self.progress.end()


    def write_unprofiled_threeprime_target_fasta(self):
        # The alignment targets for short 3' tRNA fragments
        # are trimmed profiled tRNA sequences with different 3'-acceptor variants added.
        target_fasta = fastalib.FastaOutput(self.target_fasta_path)
        for trimmed_seq in self.trimmed_trna_seqs:
            representative_id = trimmed_seq.representative_id
            for acceptor_variant in self.TARGET_ACCEPTOR_VARIANTS:
                target = trimmed_seq.string + acceptor_variant
                target_fasta.write_id(representative_id + "_" + acceptor_variant + "_" + str(len(target)))
                target_fasta.write_seq(target)
        target_fasta.close()


    def write_unprofiled_threeprime_query_fasta(self):
        # Only consider sequences shorter than the span of the acceptor through the T arm.
        # (For sequences ending 3'-CCA, the minimum feature set needed to profile the sequence as tRNA
        # goes from the acceptor through the T loop.
        # For sequences ending in variants of 3'-CCA, such as 3'-C or 3'-CCAAA,
        # the minimum feature set goes from the acceptor through the full T arm.)
        # A query must end in CCA or a recognized 3' variant.
        query_fasta = fastalib.FastaOutput(self.query_fasta_path)
        for i, unique_seq in enumerate(self.unique_nontrna_seqs):
            if len(unique_seq.string) < self.SHORT_FRAGMENT_LENGTH:
                query_fasta.write_id(unique_seq.representative_id + "_" + str(i))
                query_fasta.write_seq(unique_seq.string)
        query_fasta.close()


    def convert_sam_to_bam(self, sam_path, bam_path):
        self.progress.update("Converting SAM to BAM file")
        pysam.view('-bS', sam_path, '-F', '4', '-o', bam_path, catch_stdout=False)
        self.run.info("Raw BAM file", bam_path, quiet=(not anvio.DEBUG))


    def sort_bam(self, raw_bam_path, sorted_bam_path):
        self.progress.update("Sorting BAM file")
        # -@ is number of _additional_ threads (i.e the default is 0), so subtract 1
        pysam.sort('-o', sorted_bam_path, '-@', str(self.num_threads - 1), raw_bam_path)
        self.run.info("Sorted BAM file", sorted_bam_path, quiet=(not anvio.DEBUG))


    def index_bam(self, sorted_bam_path):
        self.progress.update("Indexing BAM file")
        pysam.index(sorted_bam_path)
        self.run.info("BAM file index", sorted_bam_path + '.bai', quiet=(not anvio.DEBUG))


    def do_alignment_steps(self,
                           query_fasta_path=None,
                           score_min='C,0,0',
                           max_alignments_per_query=1,
                           seed_length=None,
                           seed_mismatches_allowed=None,
                           rdg=None,
                           rfg=None):
        if query_fasta_path is None:
            query_fasta_path = self.query_fasta_path
        bowtie2 = Bowtie2(query_fasta=query_fasta_path, target_fasta=self.target_fasta_path, num_threads=self.num_threads)
        bowtie2.build_index()
        bowtie2.score_min = score_min # only allow end-to-end alignments with no mismatches or gaps
        # Query sequences can theoretically map to this maximum number of targets.
        bowtie2.k = max_alignments_per_query
        bowtie2.L = seed_length
        bowtie2.N = seed_mismatches_allowed
        bowtie2.rdg = rdg
        bowtie2.rfg = rfg
        bowtie2.align()
        raw_bam_path = os.path.splitext(bowtie2.sam)[0] + "-RAW.bam"
        self.convert_sam_to_bam(bowtie2.sam, raw_bam_path)
        sorted_bam_path = os.path.splitext(bowtie2.sam)[0] + "-SORTED.bam"
        self.sort_bam(raw_bam_path, sorted_bam_path)
        os.remove(raw_bam_path)
        self.index_bam(sorted_bam_path)
        return pysam.AlignmentFile(sorted_bam_path, 'rb')


    def map_threeprime(self):
        """
        When 3' tRNA fragments are not long enough to be profiled successfully
        (not spanning 3'-CCA through T loop or 3'-acceptor variant through T arm) they are left over as non-tRNA.
        We try to salvage these sequences to increase the accuracy of our inference of the relative abundances of tRNA species
        and to remove unprofiled tRNA sequences from the file of other potentially interesting RNA sequences.

        EXAMPLE:
        target tRNA  : TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG(CCA)
        mapped tRNA 1:                                                    GTTCAATTCCCCGTCGCGGAG CC
        mapped tRNA 2:                                                 GGGGTTCAATTCCCCGTCGCGGAG CCA
        """
        self.progress.new("Mapping short unprofiled tRNA fragments to profiled tRNA")

        os.mkdir(self.derep_mapping_temp_dir)

        self.write_unprofiled_threeprime_target_fasta()

        self.write_unprofiled_threeprime_query_fasta()

        bam_file = self.do_alignment_steps()

        unique_trna_seqs = []
        nontrna_indices = []
        mapped_count = 0
        for aligned_segment in bam_file.fetch():
            alignment_start = aligned_segment.reference_start
            alignment_end = aligned_segment.reference_end # pythonic stop position

            reference_id = aligned_segment.reference_name
            reference_length = int(reference_id.split('_')[-1])
            reference_id = reference_id[: -len(str(reference_length)) - 1]
            reference_acceptor_string = reference_id.split('_')[-1]
            reference_id = reference_id[: -len(reference_acceptor_string) - 1]
            acceptorless_reference_length = reference_length - len(reference_acceptor_string)

            if (alignment_end <= acceptorless_reference_length) or (alignment_start >= acceptorless_reference_length):
                # The query does not align with the reference sequence's acceptor tail.
                continue

            query_id = aligned_segment.query_name
            nontrna_index = int(query_id.split('_')[-1])
            query_id = query_id[: -len(str(nontrna_index)) - 1]
            query_string = aligned_segment.query_sequence
            query_acceptor_string = reference_acceptor_string[: alignment_end - reference_length]

            if query_acceptor_string not in THREEPRIME_VARIANTS:
                # The acceptor variant sequence is not recognized, e.g., "A" or "TC".
                continue

            unique_seq = self.unique_nontrna_seqs[nontrna_index]
            unique_seq.identification_method = 'mapped'
            unique_seq.acceptor_length = len(query_acceptor_string)
            unique_seq.extra_fiveprime_length = 0
            unique_trna_seqs.append(unique_seq)

            nontrna_indices.append(nontrna_index)
            mapped_count += unique_seq.input_count

        for nontrna_index in sorted(nontrna_indices, reverse=True):
            del self.unique_nontrna_seqs[nontrna_index]

        self.unique_trna_seqs.extend(unique_trna_seqs)

        self.progress.end() # the call below to trim_ends calls self.progress.new, so end the current progress

        self.trim_ends(unique_trna_seqs)

        shutil.rmtree(self.derep_mapping_temp_dir)

        self.run.info("Mapped short reads to 3' end of tRNA", mapped_count)


    def dereplicate_threeprime(self):
        """
        Dereplicate trimmed sequences from the 3' end of longer trimmed sequences.

        EXAMPLE:
        normalized tRNA (trimmed tRNA 1): TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        trimmed tRNA 2                  :                       AATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        trimmed tRNA 3                  :                                     GCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        """

        self.progress.new("Dereplicating trimmed tRNA sequences from the 3' end")

        # SequenceDereplicator.prefix_dereplicate returns two lists of clusters and cluster memberships.
        # We are interested in the clusters list.
        # The format of this list is as follows:
        # [(seed ID A, seed seq string A, [(trimmed seq ID A, trimmed seq string A, TrimmedSeq A),
        #                                  (trimmed seq ID B, trimmed seq string B, TrimmedSeq B), ...]),
        #  (seed ID X, ...), ...]

        self.normalized_trna_seqs = [NormalizedSeq([trimmed_seqs_info[2] for trimmed_seqs_info in normalized_seq_info[2]], skip_init=True)
                                     for normalized_seq_info
                                     in SequenceDereplicator(
                                        [trimmed_seq.representative_id for trimmed_seq in self.trimmed_trna_seqs],
                                        [trimmed_seq.string[::-1] for trimmed_seq in self.trimmed_trna_seqs], # change sequence string orientation to 3'-5' to dereplicate from 3' end
                                        extra_list=self.trimmed_trna_seqs).prefix_dereplicate()[0]]

        self.progress.end()


    def write_unprofiled_nonthreeprime_target_fasta(self):
        # Map leftover unique non-tRNA sequences to trimmed tRNA sequences
        # with extra 5' bases added when present in underlying unique tRNA sequences.
        target_fasta = fastalib.FastaOutput(self.target_fasta_path)
        for normalized_index, normalized_seq in enumerate(self.normalized_trna_seqs):
            normalized_id = normalized_seq.representative_id
            normalized_string = normalized_seq.string
            longest_trimmed_seq = normalized_seq.trimmed_seqs[0]
            if longest_trimmed_seq.unique_with_extra_fiveprime_count > 0:
                fiveprime_string_additions = []
                for unique_seq in longest_trimmed_seq.unique_seqs:
                    if unique_seq.extra_fiveprime_length > 0:
                        fiveprime_string_additions.append(unique_seq.string[: unique_seq.extra_fiveprime_length])
                fiveprime_string_additions = list(set(fiveprime_string_additions))
                for fiveprime_index, fiveprime_string in enumerate(fiveprime_string_additions):
                    # It is possible to have multiple sequences with different 5' extensions of the same length,
                    # so distinguish them by an index.
                    target_fasta.write_id(normalized_id + "_"
                                          + str(len(fiveprime_string)) + "_"
                                          + str(fiveprime_index) + "_"
                                          + str(normalized_index))
                    target_fasta.write_seq(fiveprime_string + normalized_string)
            else:
                target_fasta.write_id(normalized_id + "_0_0_" + str(normalized_index)) # no extra 5' bases
                target_fasta.write_seq(normalized_string)
        target_fasta.close()


    def write_unprofiled_nonthreeprime_query_fasta(self):
        query_fasta = fastalib.FastaOutput(self.query_fasta_path)
        for i, unique_seq in enumerate(self.unique_nontrna_seqs):
            query_fasta.write_id(unique_seq.representative_id + "_" + str(i))
            query_fasta.write_seq(unique_seq.string)
        query_fasta.close()


    def map_nonthreeprime(self):
        """
        tRNA-seq can produce transcripts that do not start at the 3' end of the tRNA molecule,
        which are left over after profiling and lumped in with non-tRNA.
        We try to salvage these sequences to increase the accuracy of our inference of the relative abundances of tRNA species
        and to remove unprofiled tRNA sequences from the file of other potentially interesting RNA sequences.

        EXAMPLE:
        normalized tRNA:                 (GT)TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        mapped tRNA 1 (extra 5' bases) :   T TCCGTGATAGTTTAATGGTCAGAATGG
        mapped tRNA 2 (interior)  :                 TAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGG
        """

        self.progress.new("Mapping unprofiled interior tRNA fragments to profiled tRNA")

        os.mkdir(self.derep_mapping_temp_dir)

        self.write_unprofiled_nonthreeprime_target_fasta()

        self.write_unprofiled_nonthreeprime_query_fasta()

        # The queries can map to multiple targets.
        bam_file = self.do_alignment_steps(max_alignments_per_query=len(self.normalized_trna_seqs))

        nontrna_indices = []
        mapping_dict = {}
        interior_mapped_count = 0
        fiveprime_mapped_count = 0
        for aligned_segment in bam_file.fetch():
            alignment_start = aligned_segment.reference_start
            alignment_end = aligned_segment.reference_end # pythonic stop position

            reference_id = aligned_segment.reference_name
            normalized_index = int(reference_id.split('_')[-1])
            reference_id = reference_id[: -len(str(normalized_index)) - 1]
            reference_fiveprime_index = reference_id.split('_')[-1]
            reference_id = reference_id[: -len(reference_fiveprime_index) - 1]
            reference_fiveprime_length = int(reference_id.split('_')[-1])
            reference_id = reference_id[: -len(str(reference_fiveprime_length)) - 1]

            normalized_end_position = alignment_end - reference_fiveprime_length
            if normalized_end_position < 0:
                # Ignore queries that align entirely to extra 5' bases.
                continue

            query_id = aligned_segment.query_name
            nontrna_index = int(query_id.split('_')[-1])
            query_id = query_id[: -len(str(nontrna_index)) - 1]

            unique_mapped_seq = self.unique_nontrna_seqs[nontrna_index]
            unique_mapped_seq.identification_method = 'mapped'
            unique_mapped_seq.acceptor_length = 0
            if reference_fiveprime_length - alignment_start > 0:
                # Assume that the extra 5' sequences are the same regardless of the reference sequence
                # the query sequence is mapped to.
                # This could be false if tRNA profiling erroneously identifies the end of the acceptor stem,
                # which seems very unlikely,
                # or if the query sequence maps to different places at the end of the acceptor stem in different tRNA,
                # which also seems very unlikely.
                unique_mapped_seq.extra_fiveprime_length = reference_fiveprime_length - alignment_start
                normalized_start_position = 0
            else:
                unique_mapped_seq.extra_fiveprime_length = 0
                normalized_start_position = alignment_start - reference_fiveprime_length

            trimmed_seq = TrimmedSeq(unique_mapped_seq.string[unique_mapped_seq.extra_fiveprime_length: ], [unique_mapped_seq])

            if normalized_index in mapping_dict:
                normalized_dict = mapping_dict[normalized_index]
                # There were multiple targets created from normalized sequences
                # for each 5' extension found among the sequences trimmed down to get the normalized sequence.
                # So a query sequence can map to the same normalized sequence
                # that formed target sequences with 5' extensions that are subsequences of each other.
                # Therefore, we need to make sure that the query maps to the normalized sequence only once.
                if nontrna_index not in normalized_dict:
                    normalized_dict[nontrna_index] = (trimmed_seq, normalized_start_position, normalized_end_position)
            else:
                mapping_dict[normalized_index] = {
                    nontrna_index: (trimmed_seq, normalized_start_position, normalized_end_position)}

            nontrna_indices.append(nontrna_index)

        for normalized_index, normalized_dict in mapping_dict.items():
            normalized_seq = self.normalized_trna_seqs[normalized_index]
            for trimmed_seq, start_position, end_position in normalized_dict.values():
                normalized_seq.trimmed_seqs.append(trimmed_seq)
                normalized_seq.start_positions.append(start_position)
                normalized_seq.end_positions.append(end_position)

            if trimmed_seq.unique_seqs[0].extra_fiveprime_length > 0:
                fiveprime_mapped_count += trimmed_seq.input_count
            else:
                interior_mapped_count += trimmed_seq.input_count

        for normalized_seq in self.normalized_trna_seqs:
            normalized_seq.init()

        for nontrna_index in sorted(set(nontrna_indices), reverse=True):
            unique_mapped_seq = self.unique_nontrna_seqs.pop(nontrna_index)
            self.trimmed_trna_seqs.append(
                TrimmedSeq(unique_mapped_seq.string[unique_mapped_seq.extra_fiveprime_length: ], [unique_mapped_seq]))

        self.progress.end()

        shutil.rmtree(self.derep_mapping_temp_dir)

        self.run.info("Mapped reads to interior of tRNA", interior_mapped_count)
        self.run.info("Mapped reads to tRNA with extra 5' bases", fiveprime_mapped_count)


    # def write_agglomerate_query_and_target(self):
    #     fasta = fastalib.FastaOutput(self.target_fasta_path)
    #     for normalized_seq in self.normalized_trna_seqs:
    #         fasta.write_id(normalized_seq.representative_id)
    #         fasta.write_seq(normalized_seq.string)
    #     fasta.close()


    # def agglomerate(self):
    #     self.progress.new("A")

    #     os.mkdir(self.derep_mapping_temp_dir)

    #     self.write_agglomerate_query_and_target()
    #     bam_file = self.do_alignment_steps(query_fasta_path=self.target_fasta_path,
    #                                        score_min=None,
    #                                        max_alignments_per_query=self.max_possible_alignments,
    #                                        seed_length=32,
    #                                        seed_mismatches_allowed=1,
    #                                        rdg='100,3',
    #                                        rfg='100,3')

    #     os.path.join(self.derep_mapping_temp_dir, "target-SORTED.bam")
    #     Agglomerator(input_fasta_path=self.target_fasta_path,
    #                  input_bam_path=os.path.join(self.derep_mapping_temp_dir, "target-SORTED.bam"),
    #                  output_fasta_path=os.path.join(self.derep_mapping_temp_dir, "agglomerated.fasta"),
    #                  output_bam_path=os.path.join(self.derep_mapping_temp_dir, "target-AGGLOMERATED.bam"),
    #                  max_possible_alignments=self.max_possible_alignments).agglomerate()

    #     self.progress.end()
    #     sys.exit()

    #     shutil.rmtree(self.derep_mapping_temp_dir)


    def calc_normalization_stats(self):
        self.progress.new("Calculating normalized tRNA statistics")

        self.progress.update("For each trimmed sequence, counting the normalized sequences containing it")
        normalized_count_dict = OrderedDict([(trimmed_seq.representative_id, 0) for trimmed_seq in self.trimmed_trna_seqs])
        for normalized_seq in self.normalized_trna_seqs:
            for trimmed_seq in self.trimmed_trna_seqs:
                normalized_count_dict[trimmed_seq.representative_id] += 1
        self.counts_of_normalized_seqs_containing_trimmed_seqs = [normalized_count for normalized_count in normalized_count_dict.values()]

        self.progress.update("Finding the \"multiplicity\" of each trimmed sequence among normalized sequences")
        multiplicity_dict = OrderedDict()
        for trimmed_seq, normalized_count_item in zip(self.trimmed_trna_seqs, normalized_count_dict.items()):
            representative_id, normalized_count = normalized_count_item
            multiplicity_dict[representative_id] = trimmed_seq.input_count * normalized_count
        self.multiplicities_of_trimmed_seqs_among_normalized_seqs = [multiplicity for multiplicity in multiplicity_dict.values()]

        self.progress.update("Finding the \"average multiplicity\" of each normalized sequence")
        for normalized_seq in self.normalized_trna_seqs:
            multiplicity_sum = 0
            for trimmed_seq in normalized_seq.trimmed_seqs:
                multiplicity_sum += multiplicity_dict[trimmed_seq.representative_id]
            self.average_multiplicities_of_normalized_seqs.append(round(multiplicity_sum / len(normalized_seq.trimmed_seqs), 1))

        self.progress.end()


    def write_trimmed_table(self):
        self.progress.new("Writing tRNA-seq database table of trimmed tRNA sequences")
        trimmed_table_entries = []
        for trimmed_seq, normalized_seq_count in zip(self.trimmed_trna_seqs, self.counts_of_normalized_seqs_containing_trimmed_seqs):
            trimmed_table_entries.append(
                (trimmed_seq.representative_id,
                 len(trimmed_seq.unique_seqs),
                 trimmed_seq.input_count,
                 trimmed_seq.string,
                 normalized_seq_count,
                 trimmed_seq.unique_with_extra_fiveprime_count,
                 trimmed_seq.input_with_extra_fiveprime_count)
                + tuple([v for v in trimmed_seq.input_acceptor_variant_count_dict.values()]))

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db._exec_many(
            '''INSERT INTO %s VALUES (%s)'''
            % ('trimmed', ','.join('?' * len(tables.trnaseq_trimmed_table_structure))),
            trimmed_table_entries)

        trimmed_seq_count = len(self.trimmed_trna_seqs)
        trnaseq_db.db.set_meta_value('num_trimmed_trna_seqs', trimmed_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()
        self.run.info("Trimmed tRNA, removing 5'/3' ends", trimmed_seq_count)


    def write_normalized_table(self):
        self.progress.new("Writing tRNA-seq database table of normalized tRNA sequences")

        normalized_table_entries = []
        threeprime_variant_keys = [threeprime_variant + 'input_seq_count' for threeprime_variant in THREEPRIME_VARIANTS]
        for normalized_seq, average_multiplicity in zip(self.normalized_trna_seqs, self.average_multiplicities_of_normalized_seqs):
            normalized_table_entries.append(
                (normalized_seq.representative_id,
                 len(normalized_seq.trimmed_seqs),
                 normalized_seq.input_count,
                 average_multiplicity,
                 normalized_seq.count_of_trimmed_seqs_mapped_to_threeprime_end,
                 normalized_seq.count_of_input_seqs_mapped_to_threeprime_end,
                 normalized_seq.count_of_trimmed_seqs_mapped_to_interior,
                 normalized_seq.count_of_input_seqs_mapped_to_interior,
                 normalized_seq.count_of_trimmed_seqs_mapped_to_fiveprime_end,
                 normalized_seq.count_of_input_seqs_mapped_to_fiveprime_end)
                + tuple(normalized_seq.input_acceptor_variant_count_dict))

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db._exec_many(
            '''INSERT INTO %s VALUES (%s)'''
            % ('normalized', ','.join('?' * len(tables.trnaseq_normalized_table_structure))),
            normalized_table_entries)

        normalized_seq_count = len(self.normalized_trna_seqs)
        trnaseq_db.db.set_meta_value('num_normalized_trna_seqs', normalized_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()
        self.run.info("Normalized tRNA, consolidating tRNA fragments", normalized_seq_count)


    # def write_agglomerated_table(self):
    #     self.progress.new("Writing tRNA-seq database table of agglomerated tRNA sequences")

    #     agglomerated_table_entries = []

    #     trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
    #     trnaseq_db.db._exec_many(
    #         '''INSERT INTO %s VALUES (%s)'''
    #         % ('agglomerated', ','.join('?' * len(tables.trnaseq_agglomerated_table_structure))),
    #         agglomerated_table_entries)

    #     agglomerated_seq_count = len(self.agglomerated_trna_seqs)
    #     trnaseq_db.db.set_meta_value('num_agglomerated_trna_seqs', agglomerated_seq_count)
    #     trnaseq_db.disconnect()

    #     self.progress.end()
    #     self.run.info("Agglomerated tRNA, from normalized tRNA", agglomerated_seq_count)


    def write_uniqued_nontrna_supplement(self):
        self.progress.new("Writing a file of sequences not identified as tRNA.")
        self.run.info("Output non-tRNA file", self.uniqued_nontrna_path)

        with open(self.uniqued_nontrna_path, 'w') as nontrna_file:
            nontrna_file.write("\t".join(self.UNIQUED_NONTRNA_HEADER) + "\n")
            for unique_seq in self.unique_nontrna_seqs:
                nontrna_file.write(unique_seq.representative_id + "\t"
                                   + ",".join(unique_seq.input_ids) + "\t"
                                   + str(unique_seq.input_count) + "\t"
                                   + unique_seq.string + "\n")

        self.progress.end()


    def write_uniqued_trna_supplement(self):
        self.progress.new("Writing a file showing which input sequences formed unique tRNA sequences")
        self.run.info("Output uniqued tRNA file", self.uniqued_trna_path)

        with open(self.uniqued_trna_path, 'w') as uniqued_file:
            uniqued_file.write("\t".join(self.UNIQUED_TRNA_HEADER) + "\n")
            for unique_seq in self.unique_trna_seqs:
                uniqued_file.write(unique_seq.representative_id + "\t"
                                   + ",".join(unique_seq.input_ids) + "\n")

        self.progress.end()


    def write_trimmed_supplement(self):
        self.progress.new("Writing a file showing how trimmed tRNA sequences were formed from unique sequences")
        self.run.info("Output trimmed tRNA file", self.trimmed_ends_path)

        with open(self.trimmed_ends_path, 'w') as trimmed_file:
            trimmed_file.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for trimmed_seq in sorted(self.trimmed_trna_seqs,
                                      key=lambda trimmed_seq: -trimmed_seq.input_count):
                representative_id = trimmed_seq.representative_id
                for unique_seq in sorted(trimmed_seq.unique_seqs,
                                         key=lambda unique_seq: (-unique_seq.extra_fiveprime_length,
                                                                 -unique_seq.acceptor_length)):
                    trimmed_file.write(representative_id + "\t"
                                       + unique_seq.representative_id + "\t"
                                       + unique_seq.string[: unique_seq.extra_fiveprime_length] + "\t"
                                       + unique_seq.string[-unique_seq.acceptor_length: ] + "\t"
                                       + str(unique_seq.input_count) + "\n")

        self.progress.end()


    def write_normalized_supplement(self):
        self.progress.new("Writing a file showing how normalized tRNA sequences were formed from trimmed sequences")
        self.run.info("Output normalized tRNA file", self.normalized_fragments_path)

        with open(self.normalized_fragments_path, 'w') as normalized_file:
            normalized_file.write("\t".join(self.NORMALIZED_FRAGMENTS_HEADER) + "\n")
            for normalized_seq in sorted(self.normalized_trna_seqs,
                                         key=lambda normalized_seq: -normalized_seq.input_count):
                representative_id = normalized_seq.representative_id
                for trimmed_seq, start_position, end_position in sorted(
                    zip(normalized_seq.trimmed_seqs, normalized_seq.start_positions, normalized_seq.end_positions),
                    key=lambda t: (t[1], -t[2])):
                    normalized_file.write(representative_id + "\t"
                                          + trimmed_seq.representative_id + "\t"
                                          + str(start_position) + "\t"
                                          + str(end_position) + "\n")

        self.progress.end()


    # def write_agglomerated_supplement(self):
    #     self.progress.new("Writing a file showing which normalized sequences agglomerated")
    #     self.run.info("Output agglomerated tRNA file", self.agglomerated_seqs_path)

    #     with open(self.agglomerated_seqs_path, 'w') as agglomerated_file:
    #         agglomerated_file.write("\t".join(self.AGGLOMERATED_SEQS_HEADER) + "\n")
    #         for agglomerated_seq in sorted(self.agglomerated_trna_seqs,
    #                                        key=lambda agglomerated_seq: -agglomerated_seq.input_count):
    #             agglomerated_file.write(agglomerated_seq.representative_id + "\t"
    #                                     + ",".join(normalized_seq.representative_id
    #                                                for normalized_seq in agglomerated_seq.normalized_seqs) + "\n")

    #     self.progress.end()


    def process(self):
        self.sanity_check()

        self.generate_trnaseq_db()
        self.profile_trna(self.get_unique_input_seqs())

        self.trim_ends(self.unique_trna_seqs)

        self.map_threeprime()

        self.dereplicate_threeprime()

        self.map_nonthreeprime()

        self.calc_normalization_stats()

        self.write_trimmed_table()
        self.write_normalized_table()

        # self.agglomerate()

        # self.write_agglomerated_table()

        self.write_uniqued_nontrna_supplement()
        self.write_uniqued_trna_supplement()
        self.write_trimmed_supplement()
        self.write_normalized_supplement()
        # self.write_agglomerated_supplement()
