# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes for tRNA-seq dataset operations. anvi-trnaseq is the default client using this module. """

import os
import itertools
import multiprocessing

from itertools import combinations
from collections import OrderedDict

import anvio
import anvio.utils as utils
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths
import anvio.trnaidentifier as trnaidentifier

from anvio.errors import ConfigError
from anvio.agglomeration import Agglomerator
from anvio.sequence import Aligner, Dereplicator
from anvio.dbops_package.trnaseq import TRNASeqDatabase
from anvio.constants_package.trnaseq import THREEPRIME_VARIANTS, TRNA_FEATURE_NAMES, DISCRIMINATOR_THROUGH_FIVEPRIME_T_STEM_STRAND_LENGTH


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


class UniqueSeq:
    __slots__ = ['seq_string',
                 'input_names',
                 'identification_method',
                 'acceptor_length',
                 'extra_fiveprime_length',
                 'input_count',
                 'representative_name']

    def __init__(self, seq_string, input_names, identification_method=None, acceptor_length=None, extra_fiveprime_length=None, skip_init=False):
        """A dereplicated tRNA-seq read, with information from tRNA feature profiling"""

        self.seq_string = seq_string
        self.input_names = input_names
        self.identification_method = identification_method # 'profiled' or 'mapped' if tRNA
        self.acceptor_length = acceptor_length
        self.extra_fiveprime_length = extra_fiveprime_length

        if skip_init:
            self.input_count = None
            self.representative_name = None
        else:
            self.init()


    def init(self):
        """Set attributes representative of a final set of dereplicated reads"""

        self.input_count = len(self.input_names)
        self.representative_name = self.input_names[0]


class TrimmedSeq:
    __slots__ = ['seq_string',
                 'unique_seqs',
                 'input_count',
                 'unique_with_extra_fiveprime_count',
                 'input_with_extra_fiveprime_count',
                 'representative_name',
                 'input_acceptor_variant_count_dict',
                 'identification_method']

    def __init__(self, seq_string, unique_seqs, skip_init=False):
        """A tRNA sequence with bases trimmed 5' of the acceptor stem and 3' of the discriminator"""

        self.seq_string = seq_string
        self.unique_seqs = unique_seqs # list of UniqueSeq objects

        if skip_init:
            self.input_count = None
            self.unique_with_extra_fiveprime_count = None
            self.input_with_extra_fiveprime_count = None
            self.representative_name = None
            self.input_acceptor_variant_count_dict = None
            self.identification_method = None
        else:
            self.init()


    def init(self):
        """Set attributes representative of a final set of input `UniqueSeq` objects"""

        self.input_count = sum([len(unique_seq.input_names) for unique_seq in self.unique_seqs])

        self.unique_with_extra_fiveprime_count = sum([1 if unique_seq.extra_fiveprime_length else 0
                                                      for unique_seq in self.unique_seqs])
        self.input_with_extra_fiveprime_count = sum([len(unique_seq.input_names) if unique_seq.extra_fiveprime_length else 0
                                                     for unique_seq in self.unique_seqs])

        # The representative name is chosen as follows:
        # 1. Most abundant full-length tRNA (no extra 5' bases), ignoring acceptor sequence
        # 2. Most abundant longer-than-full-length tRNA
        # 3. Most abundant fragmentary tRNA
        # Sort such that the first sequence is the most abundant longest and the last is the least abundant shortest.
        unique_seqs = sorted(self.unique_seqs, key=lambda unique_seq: (-unique_seq.extra_fiveprime_length, -unique_seq.input_count))

        if unique_seqs[0].extra_fiveprime_length > 0:
            # If there is also a unique sequence that was ultimately trimmed down
            # to the same sequence as the sequence with extra 5' bases, it must be a full-length sequence.
            if unique_seqs[-1].extra_fiveprime_length == 0:
                # Sort such that the last sequence is the most abundant shortest.
                representative_name = sorted(unique_seqs,
                                             key=lambda unique_seq: (-unique_seq.extra_fiveprime_length,
                                                                     unique_seq.input_count))[-1].representative_name
            else:
                representative_name = unique_seqs[0].representative_name
        else:
            # ALL unique sequences are EITHER full-length OR a fragment.
            representative_name = unique_seqs[0].representative_name

        self.representative_name = representative_name

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for unique_seq in self.unique_seqs:
            if unique_seq.acceptor_length: # unique_seq need not have an acceptor
                acceptor_seq_string = unique_seq.seq_string[-unique_seq.acceptor_length: ]
                input_acceptor_variant_count_dict[acceptor_seq_string] += unique_seq.input_count
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
    __slots__ = ['trimmed_seqs',
                 'representative_name',
                 'seq_string',
                 'start_positions',
                 'end_positions',
                 'input_count',
                 'input_with_extra_fiveprime_count',
                 'input_acceptor_variant_count_dict',
                 'count_of_trimmed_seqs_mapped_to_threeprime_end',
                 'count_of_input_seqs_mapped_to_threeprime_end',
                 'count_of_trimmed_seqs_mapped_to_interior',
                 'count_of_input_seqs_mapped_to_interior',
                 'count_of_trimmed_seqs_mapped_to_fiveprime_end',
                 'count_of_input_seqs_mapped_to_fiveprime_end']

    def __init__(self, trimmed_seqs, start_positions=None, end_positions=None, skip_init=False):
        """A longer tRNA sequence consolidated from shorter tRNA fragments"""

        self.trimmed_seqs = trimmed_seqs # list of TrimmedSeq objects
        self.representative_name = trimmed_seqs[0].representative_name
        self.seq_string = trimmed_seqs[0].seq_string
        if start_positions and end_positions:
            self.start_positions = start_positions
            self.end_positions = end_positions
        elif (not start_positions) and (not end_positions):
            # Trimmed seqs were dereplicated from the 3' end of the normalized sequence.
            normalized_seq_length = len(self.seq_string)
            self.start_positions = [normalized_seq_length - len(trimmed_seq.seq_string) for trimmed_seq in self.trimmed_seqs]
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
        """Set attributes representative of a final set of `TrimmedSeq` objects"""

        self.input_count = sum([trimmed_seq.input_count for trimmed_seq in self.trimmed_seqs])

        self.input_with_extra_fiveprime_count = sum([trimmed_seq.input_with_extra_fiveprime_count for trimmed_seq in self.trimmed_seqs])

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for trimmed_seq in self.trimmed_seqs:
            for acceptor_seq_string, input_count in trimmed_seq.input_acceptor_variant_count_dict.items():
                if input_count > 0:
                    input_acceptor_variant_count_dict[acceptor_seq_string] += input_count
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
                if end_position == len(self.seq_string):
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


class ModifiedSeq:
    __slots__ = ['normalized_seqs',
                 'representative_name',
                 'modification_indices',
                 'seq_string',
                 'input_count',
                 'input_with_extra_fiveprime_count',
                 'input_acceptor_variant_count_dict',
                 'count_of_input_seqs_mapped_to_threeprime_end',
                 'count_of_input_seqs_mapped_to_interior',
                 'count_of_input_seqs_mapped_to_fiveprime_end']

    def __init__(self, normalized_seqs, modification_dict=None, skip_init=False):
        """A tRNA sequence with sites of predicted modification-induced mutations"""

        self.normalized_seqs = normalized_seqs # list of NormalizedSeq objects
        self.representative_name = normalized_seqs[0].representative_name
        if modification_dict:
            self.modification_indices = tuple(modification_dict.keys())
            modified_seq_string = normalized_seqs[0].seq_string
            for modification_index, abundant_nucleotide in modification_dict.items():
                modified_seq_string = (modified_seq_string[: modification_index]
                                   + abundant_nucleotide
                                   + modified_seq_string[modification_index + 1: ])
            self.seq_string = modified_seq_string
        else:
            self.modification_indices = tuple()
            self.seq_string = normalized_seqs[0].seq_string

        if skip_init:
            self.input_count = None
            self.input_with_extra_fiveprime_count = None
            self.input_acceptor_variant_count_dict = None
            self.count_of_input_seqs_mapped_to_threeprime_end = None
            self.count_of_input_seqs_mapped_to_interior = None
            self.count_of_input_seqs_mapped_to_fiveprime_end = None
        else:
            self.init()


    def init(self):
        """Set attributes representative of a final set of `NormalizedSeq` objects"""

        self.input_count = sum([normalized_seq.input_count for normalized_seq in self.normalized_seqs])

        self.input_with_extra_fiveprime_count = sum([normalized_seq.input_with_extra_fiveprime_count for normalized_seq in self.normalized_seqs])

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for normalized_seq in self.normalized_seqs:
            for acceptor_seq_string, input_count in normalized_seq.input_acceptor_variant_count_dict.items():
                if input_count > 0:
                    input_acceptor_variant_count_dict[acceptor_seq_string] += input_count
        self.input_acceptor_variant_count_dict = input_acceptor_variant_count_dict

        count_of_input_seqs_mapped_to_threeprime_end = 0
        count_of_input_seqs_mapped_to_interior = 0
        count_of_input_seqs_mapped_to_fiveprime_end = 0
        for normalized_seq in self.normalized_seqs:
            count_of_input_seqs_mapped_to_threeprime_end += normalized_seq.count_of_input_seqs_mapped_to_threeprime_end
            count_of_input_seqs_mapped_to_interior += normalized_seq.count_of_input_seqs_mapped_to_interior
            count_of_input_seqs_mapped_to_fiveprime_end += normalized_seq.count_of_input_seqs_mapped_to_fiveprime_end
        self.count_of_input_seqs_mapped_to_threeprime_end = count_of_input_seqs_mapped_to_threeprime_end
        self.count_of_input_seqs_mapped_to_interior = count_of_input_seqs_mapped_to_interior
        self.count_of_input_seqs_mapped_to_fiveprime_end = count_of_input_seqs_mapped_to_fiveprime_end


class TRNASeqDataset:
    UNIQUED_NONTRNA_HEADER = ["representative_name", "input_names", "input_count", "sequence"]
    UNIQUED_TRNA_HEADER = ["representative_name", "input_names"]
    TRIMMED_ENDS_HEADER = ["representative_name", "unique_name", "fiveprime_sequence", "threeprime_sequence", "input_seq_count"]
    NORMALIZED_FRAGMENTS_HEADER = ["representative_name", "trimmed_name", "start", "end"]

    SHORT_FRAGMENT_LENGTH = 30 # number of nucleotides from 3'-CCANN through T arm
    TARGET_ACCEPTOR_VARIANTS = ['CCC', 'CCG', 'CCT', # used in mapping fragments to tRNA interior/5' end
                                'CAA', 'CGA', 'CTA',
                                'ACA', 'GCA', 'TCA',
                                'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT',
                                'CCACA', 'CCACC', 'CCACG', 'CCACT',
                                'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT',
                                'CCATA', 'CCATC', 'CCATG', 'CCATT']

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        """Class for processing tRNA-seq dataset"""

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

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using --fasta-file or -f.")
        if not self.project_name:
            raise ConfigError("Please set a project name using --project-name or -n.")
        if not self.output_dir:
            raise ConfigError("Please provide a non-existent output directory using --output-dir or -o.")

        self.input_seq_count = 0

        self.trnaseq_db_path = os.path.join(self.output_dir, self.project_name + "-TRNASEQ.db")

        self.profiled_unique_seq_count = 0 # tracks reads placed in profile input queue
        self.retrieved_profile_count = 0 # tracks reads retrieved from profile output queue

        self.trna_count = 0
        self.unique_trna_count = 0
        self.trna_containing_anticodon_count = 0
        self.full_length_trna_count = 0
        self.trna_with_one_to_three_extra_fiveprime_bases_count = 0
        self.trna_with_more_than_three_extra_fiveprime_bases_count = 0
        self.trna_with_extrapolated_fiveprime_feature_count = 0
        self.trna_with_threeprime_cca_count = 0
        self.trna_with_threeprime_cc_count = 0
        self.trna_with_threeprime_c_count = 0
        self.trna_with_threeprime_nca_cna_ccn_count = 0
        self.trna_with_threeprime_ccan_ccann_count = 0

        self.uniqued_nontrna_path = os.path.join(self.output_dir, self.project_name + "-UNIQUED_NONTRNA.txt")
        self.uniqued_trna_path = os.path.join(self.output_dir, self.project_name + "-UNIQUED_TRNA.txt")
        self.trimmed_ends_path = os.path.join(self.output_dir, self.project_name + "-TRIMMED_ENDS.txt")
        self.normalized_fragments_path = os.path.join(self.output_dir, self.project_name + "-NORMALIZED_FRAGMENTS.txt")

        self.log_path = os.path.join(self.output_dir, "log.txt")

        self.unique_nontrna_seqs = []
        self.unique_trna_seqs = []
        self.trimmed_trna_seqs = []
        self.normalized_trna_seqs = []
        self.modified_trna_seqs = []

        self.counts_of_normalized_seqs_containing_trimmed_seqs = [] # same length as self.trimmed_trna_seqs
        # "Multiplicity" of a trimmed tRNA sequence
        # = number of normalized sequences that the sequence is in * number of input sequences represented by the trimmed sequence
        self.multiplicities_of_trimmed_seqs_among_normalized_seqs = [] # same length as self.trimmed_trna_seqs
        self.average_multiplicities_of_normalized_seqs = [] # same length as self.normalized_trna_seqs
        self.average_multiplicities_of_modified_seqs = []


    def sanity_check(self):
        """Check inputs before proceeding."""

        if os.path.exists(self.output_dir):
            raise ConfigError("The directory that was specified by --output-dir or -o, %s, already exists. "
                              "Please try again with a non-existent directory." % self.output_dir)
        os.mkdir(self.output_dir)
        filesnpaths.is_output_dir_writable(self.output_dir)

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = self.description_file_path.read()
        else:
            self.description = None
        self.run.info("Description", self.description_file_path if self.description_file_path else "No description given")

        self.run.log_path = self.log_path

        if not 1 < self.num_threads < multiprocessing.cpu_count():
            ConfigError("The number of threads to use must be a positive integer "
                        "less than or equal to %d. Try again!" % multiprocessing.cpu_count())

        if self.write_buffer_size < 1:
            ConfigError("The write buffer size must be a positive integer. Try again!")

        self.run.info("Input FASTA file", self.input_fasta_path)

        utils.check_fasta_id_uniqueness(self.input_fasta_path)

        if not self.skip_fasta_check:
            self.progress.new("Checking input FASTA defline format")
            self.progress.update("...")

            utils.check_fasta_id_formatting(self.input_fasta_path)

            self.progress.end()

            self.run.info_single("FASTA deflines were found to be anvi'o-compliant", mc='green')


    def create_trnaseq_db(self):
        """Create an unfilled tRNA-seq database."""

        meta_values = {'project_name': self.project_name,
                       'description': self.description if self.description else '_No description is provided_'}
        TRNASeqDatabase(self.trnaseq_db_path, quiet=False).create(meta_values)


    def get_unique_input_seqs(self):
        """Dereplicate input reads"""

        self.progress.new("Finding unique input sequences")
        self.progress.update("...")

        fasta = fastalib.SequenceSource(self.input_fasta_path)
        names = []
        seqs = []
        input_seq_count = 0
        while next(fasta):
            names.append(fasta.id)
            seqs.append(fasta.seq)
            input_seq_count += 1
        fasta.close()
        self.input_seq_count = input_seq_count

        clusters = Dereplicator(names, seqs, progress=self.progress).full_length_dereplicate()

        unique_input_seqs = [UniqueSeq(cluster.representative_seq_string, cluster.member_names) for cluster in clusters]

        self.progress.end()

        return unique_input_seqs


    def profile_trna(self, unique_input_seqs):
        """Profile tRNA features in input sequences.

        Parameters
        ==========
        unique_input_seqs : list
            List of UniqueSeq objects
        """

        self.progress.new("Profiling input sequences for tRNA features")
        self.progress.update("...")

        write_points = self.get_write_points(len(unique_input_seqs))
        write_point_index = 0

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        processes = [multiprocessing.Process(target=trnaidentifier.profile_wrapper, args=(input_queue, output_queue))
                     for _ in range(self.num_threads)]
        for process in processes:
            process.start()


        chunk_dict = {} # map unique seq names to UniqueSeq objects
        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        for unique_seq in unique_input_seqs:
            # Profile each unique input sequence.
            input_queue.put((unique_seq.representative_name, unique_seq.seq_string))
            self.profiled_unique_seq_count += 1
            chunk_dict[unique_seq.representative_name] = unique_seq

            if self.profiled_unique_seq_count != write_points[write_point_index]:
                # Keep adding sequences to the profiling queue until the write point is hit.
                continue

            unique_trna_seqs_chunk, unique_nontrna_seqs_chunk = self.write_profile_results(output_queue, chunk_dict, trnaseq_db)
            self.unique_trna_seqs.extend(unique_trna_seqs_chunk)
            self.unique_nontrna_seqs.extend(unique_nontrna_seqs_chunk)
            chunk_dict = {}

            self.progress.update("%d of %d unique sequences have been profiled"
                                 % (self.profiled_unique_seq_count, len(unique_input_seqs)))

            write_point_index += 1
            if write_point_index == len(write_points):
                break

        for process in processes:
            process.terminate()

        trnaseq_db.db.set_meta_value('num_input_reads_processed', self.input_seq_count)
        trnaseq_db.db.set_meta_value('num_trna_reads', self.trna_count)
        trnaseq_db.db.set_meta_value('num_unique_trna_seqs', self.unique_trna_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_containing_anticodon', self.trna_containing_anticodon_count)
        trnaseq_db.db.set_meta_value('num_full_length_trna_reads', self.full_length_trna_count)
        trnaseq_db.db.set_meta_value('num_trna_with_one_to_three_extra_fiveprime_bases', self.trna_with_one_to_three_extra_fiveprime_bases_count)
        trnaseq_db.db.set_meta_value('num_trna_with_more_than_three_extra_fiveprime_bases', self.trna_with_more_than_three_extra_fiveprime_bases_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_with_extrapolated_fiveprime_feature', self.trna_with_extrapolated_fiveprime_feature_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_cca', self.trna_with_threeprime_cca_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_cc', self.trna_with_threeprime_cc_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_c', self.trna_with_threeprime_c_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_nca_cna_ccn', self.trna_with_threeprime_nca_cna_ccn_count)
        trnaseq_db.db.set_meta_value('num_trna_reads_with_threeprime_ccan_ccann', self.trna_with_threeprime_ccan_ccann_count)
        trnaseq_db.disconnect()

        # Profiled seqs were added to the output queue as they were processed, so sort by name.
        self.unique_trna_seqs.sort(key=lambda unique_seq: unique_seq.representative_name)
        self.unique_nontrna_seqs.sort(key=lambda unique_seq: unique_seq.representative_name)

        self.progress.end()

        self.run.info("Reads processed", self.input_seq_count)
        self.run.info("Reads profiled as tRNA", self.trna_count)
        self.run.info("Unique profiled tRNA sequences", self.unique_trna_count)
        self.run.info("Profiled reads with anticodon", self.trna_containing_anticodon_count)
        self.run.info("Profiled reads spanning acceptor stem", self.full_length_trna_count)
        self.run.info("Profiled reads with 1-3 extra 5' bases", self.trna_with_one_to_three_extra_fiveprime_bases_count)
        self.run.info("Profiled reads with >3 extra 5' bases", self.trna_with_more_than_three_extra_fiveprime_bases_count)
        self.run.info("Profiled reads with extrapolated 5' feature", self.trna_with_extrapolated_fiveprime_feature_count)
        self.run.info("Profiled reads ending in 3'-CCA", self.trna_with_threeprime_cca_count)
        self.run.info("Profiled reads ending in 3'-CC", self.trna_with_threeprime_cc_count)
        self.run.info("Profiled reads ending in 3'-C", self.trna_with_threeprime_c_count)
        self.run.info("Profiled reads ending in 3'-NCA/CNA/CCN", self.trna_with_threeprime_nca_cna_ccn_count)
        self.run.info("Profiled reads ending in 3'-CCAN/CCANN", self.trna_with_threeprime_ccan_ccann_count)


    def get_write_points(self, item_count):
        """Helper function to determine the points when to write multiprocessed outputs.
        Write points are cumulative numbers of items processed.
        """

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


    def write_profile_results(self, output_queue, chunk_dict, trnaseq_db):
        """Helper function for `get_unique_input_seqs` to write profiling results to database."""

        # List of entries for each table
        trnaseq_sequences_table_entries = []
        trnaseq_info_table_entries = []
        trnaseq_features_table_entries = []
        trnaseq_unconserved_table_entries = []
        trnaseq_unpaired_table_entries = []

        unique_trna_seqs = []
        unique_nontrna_seqs = []

        while self.retrieved_profile_count < self.profiled_unique_seq_count:
            # Retrieve profiles from the output queue and write tRNA profiles to the database.
            trna_profile = output_queue.get()
            self.retrieved_profile_count += 1

            output_name = trna_profile.name
            output_seq = trna_profile.input_seq

            if not trna_profile.is_predicted_trna:
                unique_nontrna_seqs.append(chunk_dict[output_name])
                continue

            unique_seq = chunk_dict[output_name]
            unique_seq.identification_method = 'profiled'
            unique_seq.acceptor_length = len(trna_profile.acceptor_variant_string)
            unique_seq.extra_fiveprime_length = trna_profile.num_extra_fiveprime
            unique_trna_seqs.append(unique_seq)

            self.unique_trna_count += 1
            output_seq_length = len(output_seq)

            num_replicates = len(unique_seq.input_names)
            self.trna_count += num_replicates

            # Recover nucleotides that did not fit expectation,
            # either by not being the expected nucleotide or type of nucleotide
            # or by not base pairing in a stem.
            unconserved_info = trna_profile.get_unconserved_positions()
            unpaired_info = trna_profile.get_unpaired_positions()

            if trna_profile.anticodon_seq:
                self.trna_containing_anticodon_count += num_replicates
            if trna_profile.has_complete_feature_set:
                self.full_length_trna_count += num_replicates
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

            if trna_profile.num_extra_fiveprime > 3:
                self.trna_with_more_than_three_extra_fiveprime_bases_count += num_replicates
            elif trna_profile.num_extra_fiveprime > 0:
                self.trna_with_one_to_three_extra_fiveprime_bases_count += num_replicates

            trnaseq_sequences_table_entries.append((output_name, num_replicates, output_seq))

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
                (output_name,
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
                (output_name, )
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
                trnaseq_unconserved_table_entries.append((output_name, ) + unconserved_tuple)

            for unpaired_tuple in unpaired_info:
                trnaseq_unpaired_table_entries.append((output_name, ) + unpaired_tuple)

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


    def trim_ends(self, unique_trna_seqs, report_progress=True):
        """Trim any nucleotides 5' of the acceptor stem and 3' of the discriminator.

        Appends TrimmedSeq objects formed from input UniqueSeq objects to `self.trimmed_trna_seqs`

        Parameters
        ==========
        unique_trna_seqs : list
            List of UniqueSeq objects

        report_progress : bool, True
            Whether to start a new Progress report
        """

        if report_progress:
            self.progress.new("Trimming the 3' and 5' ends of sequences")
            self.progress.update("...")

        representative_names = [unique_seq.representative_name for unique_seq in unique_trna_seqs]
        trimmed_seq_strings = [unique_seq.seq_string[unique_seq.extra_fiveprime_length: -unique_seq.acceptor_length]
                               for unique_seq in unique_trna_seqs]

        clusters = Dereplicator(representative_names,
                                trimmed_seq_strings,
                                extras=unique_trna_seqs,
                                progress=self.progress).full_length_dereplicate()

        trimmed_seqs = [TrimmedSeq(cluster.representative_seq_string, cluster.member_extras) for cluster in clusters]

        self.trimmed_trna_seqs.extend(trimmed_seqs)

        if report_progress:
            self.progress.end()


    def map_threeprime(self):
        """Map any unprofiled short 3' fragments to longer profiled tRNA sequences.

        Input sequences not successfully profiled as tRNA may include 3' tRNA fragments not long enough for profiling,
        i.e., fragments not spanning the 3' end through the T arm.
        Try to salvage some of these sequences, excluding exceedingly short fragments,
        i.e., fragments not spanning the 3' acceptor (variant) through the 5' strand of the T loop.

        EXAMPLE:
        target profiled trimmed tRNA    : TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG(CCA)
        mapped unprofiled tRNA 1:                                                            GTTCAATTCCCCGTCGCGGAG CC
        mapped unprofiled tRNA 2:                                                          GGGTTCAATTCCCCGTCGCGGAG CCA
        """

        self.progress.new("Mapping short leftovers to 3' ends of profiled tRNA")

        self.progress.update("Finding query sequences")
        query_ids = []
        query_seq_strings = []
        for nontrna_index, unique_seq in enumerate(self.unique_nontrna_seqs):
            unique_seq_string = unique_seq.seq_string

            # Unique sequences may end in an ambiguous 3' acceptor variant.
            # An ending of CCC may be CCC, CC, or C.
            # An ending of CCACC may be CCACC, CC, or C.
            # Therefore, make a separate query for each possibility by trimming the different possible variants.
            # These possibilities are resolved when one and only one of them matches a trimmed profiled tRNA sequence.
            for threeprime_variant_seq_string in THREEPRIME_VARIANTS:
                if unique_seq_string[-len(threeprime_variant_seq_string): ] == threeprime_variant_seq_string:
                    query_ids.append((nontrna_index, threeprime_variant_seq_string))
                    # Reverse the sequence string for prefix dereplication.
                    query_seq_strings.append(unique_seq_string[: -len(threeprime_variant_seq_string)][::-1])

        self.progress.update("Finding target sequences")
        target_names = []
        target_seq_strings = []
        for trimmed_seq in self.trimmed_trna_seqs:
            target_names.append(trimmed_seq.representative_name)
            # Reverse the sequence string for prefix dereplication.
            target_seq_strings.append(trimmed_seq.seq_string[::-1])

        results = Aligner(query_ids,
                          query_seq_strings,
                          target_names,
                          target_seq_strings,
                          num_threads=self.num_threads,
                          progress=self.progress).match_prefixes(max_matches_per_query=1)

        self.progress.update("Evaluating mapping results")
        # Since there can be multiple queries with different 3' acceptor variants derived from the same unique sequence,
        # a chunk of results for the unique sequence is processed then evaluated after moving to the next unique sequence.

        num_mapped_nontrna_reads = 0
        nontrna_indices_to_remove = []
        new_unique_trna_seqs = []

        if query_ids: # just in case
            if results[0]:
                num_matching_nontrna = 1 # we only care that the query matched ANY target or not
                min_nontrna_seq_length = len(query_seq_strings[0])
            else: # an empty list since the query did not match any targets
                num_matching_nontrna = 0
                min_nontrna_seq_length = float('inf')
            prev_nontrna_index = query_ids[0][0]
            prev_threeprime_variant_seq_string = query_ids[0][1]

        for query_id, query_seq_string, matched_target_names in zip(query_ids[1:], query_seq_strings[1:], results[1:]):
            nontrna_index, threeprime_variant_seq_string = query_id

            if nontrna_index == prev_nontrna_index:
                if matched_target_names:
                    num_matching_nontrna += 1
                    min_nontrna_seq_length = min(min_nontrna_seq_length, len(query_seq_string))
                continue

            if num_matching_nontrna == 1:
                if min_nontrna_seq_length >= DISCRIMINATOR_THROUGH_FIVEPRIME_T_STEM_STRAND_LENGTH:
                    unique_seq = self.unique_nontrna_seqs[prev_nontrna_index]
                    unique_seq.identification_method = 'mapped'
                    unique_seq.acceptor_length = len(prev_threeprime_variant_seq_string)
                    unique_seq.extra_fiveprime_length = 0
                    new_unique_trna_seqs.append(unique_seq)
                    nontrna_indices_to_remove.append(prev_nontrna_index)
                    num_mapped_nontrna_reads += unique_seq.input_count

            if matched_target_names:
                num_matching_nontrna = 1
                min_nontrna_seq_length = len(query_seq_string)
            else:
                num_matching_nontrna = 0
                min_nontrna_seq_length = float('inf')
            prev_nontrna_index = nontrna_index
            prev_threeprime_variant_seq_string = threeprime_variant_seq_string

        # Handle the last query.
        if query_ids: # just in case
            if num_matching_nontrna == 1:
                if min_nontrna_seq_length >= DISCRIMINATOR_THROUGH_FIVEPRIME_T_STEM_STRAND_LENGTH:
                    unique_seq = self.unique_nontrna_seqs[prev_nontrna_index]
                    unique_seq.identification_method = 'mapped'
                    unique_seq.acceptor_length = len(prev_threeprime_variant_seq_string)
                    unique_seq.extra_fiveprime_length = 0
                    new_unique_trna_seqs.append(unique_seq)
                    nontrna_indices_to_remove.append(prev_nontrna_index)
                    num_mapped_nontrna_reads += unique_seq.input_count

        for nontrna_index in sorted(nontrna_indices_to_remove, reverse=True):
            del self.unique_nontrna_seqs[nontrna_index]

        self.unique_trna_seqs.extend(new_unique_trna_seqs)

        # Remake TrimmedSeq objects including new UniqueSeq objects.
        self.progress.update("Trimming sequences")
        self.trim_ends(new_unique_trna_seqs, report_progress=False)

        self.progress.end()

        self.run.info("Mapped short reads to 3' end of tRNA", num_mapped_nontrna_reads)


    def dereplicate_threeprime(self):
        """Dereplicate trimmed tRNA sequences from the 3' end of longer trimmed sequences.

        EXAMPLE:
        normalized tRNA (trimmed tRNA 1): TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        trimmed tRNA 2                  :                       AATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        trimmed tRNA 3                  :                                     GCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        """

        self.progress.new("Dereplicating trimmed tRNA sequences from the 3' end")
        self.progress.update("...")

        representative_names = [trimmed_seq.representative_name for trimmed_seq in self.trimmed_trna_seqs]
        # Reverse sequence orientation to dereplicate from the 3' end.
        reversed_seq_strings = [trimmed_seq.seq_string[::-1] for trimmed_seq in self.trimmed_trna_seqs]
        clusters = Dereplicator(representative_names,
                                reversed_seq_strings,
                                extras=self.trimmed_trna_seqs,
                                progress=self.progress).prefix_dereplicate()

        self.normalized_trna_seqs = [NormalizedSeq(cluster.member_extras, skip_init=True) for cluster in clusters]

        self.progress.end()


    def map_nonthreeprime(self):
        """Map unprofiled tRNA fragments not from the 3' end to longer profiled tRNA sequences.

        tRNA fragments may not start at the 3' end of the tRNA molecule and thereby go unprofiled.
        Try to salvage fragments from the interior and 5' end of tRNA by mapping to profiled tRNA.

        EXAMPLE:
        normalized tRNA:                 (GT)TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        mapped tRNA 1 (extra 5' bases) :   T TCCGTGATAGTTTAATGGTCAGAATGG
        mapped tRNA 2 (interior)  :                 TAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGG
        """

        self.progress.new("Mapping leftovers to the interior/5' end of profiled tRNA")

        self.progress.update("Finding query sequences")
        nonthreeprime_query_names = []
        nonthreeprime_query_seqs = []
        for nontrna_index, unique_seq in enumerate(self.unique_nontrna_seqs):
            # Avoid "very short" sequences, using the same length threshold as in `map_threeprime`.
            if len(unique_seq.seq_string) >= DISCRIMINATOR_THROUGH_FIVEPRIME_T_STEM_STRAND_LENGTH:
                nonthreeprime_query_names.append((unique_seq.representative_name, nontrna_index))
                nonthreeprime_query_seqs.append(unique_seq.seq_string)

        self.progress.update("Finding target sequences")
        # Leftover non-tRNA sequences are mapped to normalized tRNA sequences
        # with extra 5' bases added when present in underlying unique tRNA sequences.
        # Multiple targets for each normalized sequence are therefore produced for different 5' sequences.
        nonthreeprime_target_names = []
        nonthreeprime_target_seqs = []
        for normalized_index, normalized_seq in enumerate(self.normalized_trna_seqs):
            normalized_name = normalized_seq.representative_name
            normalized_seq_string = normalized_seq.seq_string
            # The longest trimmed sequence (the first in the list) is by design
            # the only one of the trimmed sequences forming the normalized sequence that may have extra 5' bases.
            longest_trimmed_seq = normalized_seq.trimmed_seqs[0]
            if longest_trimmed_seq.unique_with_extra_fiveprime_count > 0:
                fiveprime_seq_string_additions = []
                for unique_seq in longest_trimmed_seq.unique_seqs:
                    if unique_seq.extra_fiveprime_length > 0:
                        fiveprime_seq_string_additions.append(unique_seq.seq_string[: unique_seq.extra_fiveprime_length])
                fiveprime_seq_string_additions = sorted(list(set(fiveprime_seq_string_additions)), key=lambda s: len(s))
                for fiveprime_index, fiveprime_seq_string in enumerate(fiveprime_seq_string_additions):
                    # It is possible to have multiple sequences with different 5' extensions of the same length,
                    # so distinguish them by an index.
                    nonthreeprime_target_names.append((normalized_name, len(fiveprime_seq_string), fiveprime_index, normalized_index))
                    nonthreeprime_target_seqs.append(fiveprime_seq_string + normalized_seq_string)
            else:
                nonthreeprime_target_names.append((normalized_name, 0, 0, normalized_index)) # no extra 5' bases
                nonthreeprime_target_seqs.append(normalized_seq_string)

        aligned_query_dict, aligned_target_dict = Aligner(nonthreeprime_query_names,
                                                          nonthreeprime_query_seqs,
                                                          nonthreeprime_target_names,
                                                          nonthreeprime_target_seqs,
                                                          num_threads=self.num_threads,
                                                          progress=self.progress).align(max_mismatch_freq=0)

        self.progress.update("Evaluating mapping results")
        nontrna_indices = []
        mapping_dict = {}
        for query_name, aligned_query in aligned_query_dict.items():
            for alignment in aligned_query.alignments:
                reference_alignment_start = alignment.target_start
                reference_alignment_end = alignment.target_start + alignment.alignment_length

                (reference_name,
                 reference_fiveprime_length,
                 reference_fiveprime_index,
                 normalized_index) = alignment.aligned_target.name

                normalized_end_position = reference_alignment_end - reference_fiveprime_length
                if normalized_end_position < 0:
                    # Ignore queries that align entirely to extra 5' bases.
                    continue

                nontrna_index = query_name[1]

                unique_mapped_seq = self.unique_nontrna_seqs[nontrna_index]
                unique_mapped_seq.identification_method = 'mapped'
                unique_mapped_seq.acceptor_length = 0
                if reference_fiveprime_length - reference_alignment_start > 0:
                    # Assume that the extra 5' sequences are the same
                    # regardless of the reference this query is mapped to.
                    # This could be false if tRNA profiling erroneously identifies the end of the acceptor stem,
                    # which seems very unlikely,
                    # or if the query maps to different places at the end of the acceptor stem in different tRNA,
                    # which also seems very unlikely.
                    unique_mapped_seq.extra_fiveprime_length = reference_fiveprime_length - reference_alignment_start
                    normalized_start_position = 0
                else:
                    unique_mapped_seq.extra_fiveprime_length = 0
                    normalized_start_position = reference_alignment_start - reference_fiveprime_length

                trimmed_seq = TrimmedSeq(unique_mapped_seq.seq_string[unique_mapped_seq.extra_fiveprime_length: ],
                                         [unique_mapped_seq])

                if normalized_index in mapping_dict:
                    normalized_dict = mapping_dict[normalized_index]
                    # There were multiple targets created from normalized sequences
                    # for each 5' extension found among the sequences trimmed down to get the normalized sequence.
                    # So a query sequence can map to the same normalized sequence
                    # via target sequences with 5' extensions that are subsequences of each other.
                    # Therefore, we need to make sure that the query maps to the normalized sequence only once.
                    if nontrna_index not in normalized_dict:
                        normalized_dict[nontrna_index] = (trimmed_seq, normalized_start_position, normalized_end_position)
                else:
                    mapping_dict[normalized_index] = {
                        nontrna_index: (trimmed_seq, normalized_start_position, normalized_end_position)}

                nontrna_indices.append(nontrna_index)

        for normalized_index, normalized_dict in mapping_dict.items():
            normalized_seq = self.normalized_trna_seqs[normalized_index]
            for trimmed_seq, normalized_start_position, normalized_end_position in normalized_dict.values():
                normalized_seq.trimmed_seqs.append(trimmed_seq)
                normalized_seq.start_positions.append(normalized_start_position)
                normalized_seq.end_positions.append(normalized_end_position)

        for normalized_seq in self.normalized_trna_seqs:
            normalized_seq.init()

        interior_mapped_count = 0
        fiveprime_mapped_count = 0
        for nontrna_index in sorted(set(nontrna_indices), reverse=True):
            unique_mapped_seq = self.unique_nontrna_seqs.pop(nontrna_index)
            self.trimmed_trna_seqs.append(
                TrimmedSeq(unique_mapped_seq.seq_string[unique_mapped_seq.extra_fiveprime_length: ], [unique_mapped_seq]))

            if unique_mapped_seq.extra_fiveprime_length > 0:
                fiveprime_mapped_count += unique_mapped_seq.input_count
            else:
                interior_mapped_count += unique_mapped_seq.input_count

        self.progress.end()

        self.run.info("Mapped reads to interior of tRNA", interior_mapped_count)
        self.run.info("Mapped reads to tRNA with extra 5' bases", fiveprime_mapped_count)


    def find_modifications(self):
        self.progress.new("Finding modifications")

        normalized_names, normalized_seq_strings = zip(*[(normalized_seq.representative_name, normalized_seq.seq_string)
                                                  for normalized_seq in self.normalized_trna_seqs])

        agglomerator = Agglomerator(normalized_names,
                                    normalized_seq_strings,
                                    max_mismatch_freq=2/75,
                                    num_threads=self.num_threads,
                                    progress=self.progress)
        _, agglomerated_aligned_reference_dict = agglomerator.agglomerate()

        self.progress.update("Separating modification-induced mutations from \"inter-strain\" variants")

        normalized_seq_dict = {normalized_seq.representative_name: normalized_seq
                               for normalized_seq in self.normalized_trna_seqs}
        normalized_names_in_multiple_modified_seqs = []
        unique_seq_input_counts = []
        for reference_name, aligned_reference in agglomerated_aligned_reference_dict.items():
            reference_seq_string = aligned_reference.seq_string
            reference_length = len(reference_seq_string)

            query_normalized_seqs = []
            query_input_counts = []
            padded_query_seq_strings = []
            query_seq_starts = []
            for alignment in aligned_reference.alignments:
                if reference_length - alignment.target_start - alignment.alignment_length == 0:
                    # Normalized sequences should only align at the 3' end.
                    # Avoid false positive alignments in the interior of a sequence.
                    query_normalized_seq = normalized_seq_dict[alignment.aligned_query.name]
                    query_normalized_seqs.append(query_normalized_seq)
                    query_input_counts.append(query_normalized_seq.input_count)
                    padded_query_seq_strings.append(' ' * alignment.target_start + alignment.aligned_query.seq_string)
                    query_seq_starts.append(alignment.target_start)

            reference_normalized_seq = normalized_seq_dict[reference_name]

            if not padded_query_seq_strings:
                continue

            nucleotide_counts = []
            for reference_nucleotide in reference_seq_string:
                nucleotide_counts.append({reference_nucleotide: reference_normalized_seq.input_count})
            for padded_query_seq_string, query_input_count in zip(padded_query_seq_strings, query_input_counts):
                for query_nucleotide, nucleotide_count_dict in zip(padded_query_seq_string, nucleotide_counts):
                    if query_nucleotide == ' ':
                        continue
                    if query_nucleotide in nucleotide_count_dict:
                        nucleotide_count_dict[query_nucleotide] += query_input_count
                    else:
                        nucleotide_count_dict[query_nucleotide] = query_input_count

            candidate_modification_indices = []
            for reference_index, nucleotide_count_dict in enumerate(nucleotide_counts):
                # modification test
                if len(nucleotide_count_dict) > 2:
                    input_counts = sorted(nucleotide_count_dict.values(), reverse=True)
                    total_input_count = sum(input_counts)
                    if total_input_count >= 50:
                        if sum(input_counts[1: ]) / total_input_count >= 0.1:
                            candidate_modification_indices.append(reference_index)

            modification_configurations = []
            for num_modifications in range(len(candidate_modification_indices), 0, -1):
                for modification_configuration in combinations(candidate_modification_indices, num_modifications):
                    modification_configurations.append(modification_configuration)

            sorted_normalized_seq_items = sorted(zip([reference_normalized_seq] + query_normalized_seqs,
                                                     [reference_seq_string] + padded_query_seq_strings,
                                                     [0] + query_seq_starts),
                                                 key=lambda normalized_seq_item: (-len(normalized_seq_item[0].seq_string),
                                                                                  -normalized_seq_item[0].input_count))
            for modification_configuration in modification_configurations:
                modified_candidate_superdict = OrderedDict()
                normalized_seq_index = 0
                for normalized_seq, padded_seq_string, seq_start in sorted_normalized_seq_items:
                    wildcard_seq_string = padded_seq_string
                    for modification_index in modification_configuration:
                        wildcard_seq_string = (wildcard_seq_string[: modification_index]
                                               + '*'
                                               + wildcard_seq_string[modification_index + 1: ])

                    is_seq_clustered = False
                    for wildcard_seed_seq_string, modified_candidate_dicts in modified_candidate_superdict.items():
                        if wildcard_seq_string[seq_start: ] == wildcard_seed_seq_string[seq_start: ]:
                            modified_candidate_dicts.append({'normalized_seq': normalized_seq,
                                                             'padded_seq_string': padded_seq_string,
                                                             'seq_start': seq_start,
                                                             'normalized_seq_index': normalized_seq_index})
                            is_seq_clustered = True

                    if not is_seq_clustered:
                        modified_candidate_superdict[wildcard_seq_string] = [{'normalized_seq': normalized_seq,
                                                                              'padded_seq_string': padded_seq_string,
                                                                              'seq_start': seq_start,
                                                                              'normalized_seq_index': normalized_seq_index}]
                    normalized_seq_index += 1

                normalized_seq_indices_assigned_to_modified_seqs = []
                for modified_candidate_dicts in modified_candidate_superdict.values():
                    abundant_nucleotides = []
                    for modification_index in modification_configuration:
                        nucleotide_count_dict = {}
                        for modified_candidate_dict in modified_candidate_dicts:
                            nucleotide = modified_candidate_dict['padded_seq_string'][modification_index]
                            if nucleotide == ' ':
                                continue
                            if nucleotide in nucleotide_count_dict:
                                nucleotide_count_dict[nucleotide] += modified_candidate_dict['normalized_seq'].input_count
                            else:
                                nucleotide_count_dict[nucleotide] = modified_candidate_dict['normalized_seq'].input_count
                        if not nucleotide_count_dict:
                            break
                        (sorted_nucleotides,
                         sorted_input_counts) = zip(*sorted(nucleotide_count_dict.items(),
                                                            key=lambda nucleotide_count_item: -nucleotide_count_item[1]))
                        total_input_count = sum(sorted_input_counts)
                        if total_input_count < 50:
                            break
                        if len(sorted_input_counts) <= 2:
                            break
                        if sum(sorted_input_counts[1: ]) / total_input_count < 0.1:
                            break
                        abundant_nucleotides.append(sorted_nucleotides[0])
                    else:
                        normalized_seqs = [modified_candidate_dict['normalized_seq']
                                           for modified_candidate_dict in modified_candidate_dicts]
                        if normalized_seqs[0].representative_name in normalized_names_in_multiple_modified_seqs:
                            break
                        seed_seq_start = modified_candidate_dicts[0]['seq_start']
                        modification_dict = OrderedDict([(modification_index - seed_seq_start, abundant_nucleotide)
                                                         for modification_index, abundant_nucleotide
                                                         in zip(modification_configuration, abundant_nucleotides)])
                        self.modified_trna_seqs.append(ModifiedSeq(normalized_seqs, modification_dict=modification_dict))
                        for normalized_seq in normalized_seqs:
                            if normalized_seq.representative_name not in normalized_names_in_multiple_modified_seqs:
                                for trimmed_seq in normalized_seq.trimmed_seqs:
                                    for unique_seq in trimmed_seq.unique_seqs:
                                        unique_seq_input_counts.append((unique_seq.representative_name, unique_seq.input_count))
                                normalized_names_in_multiple_modified_seqs.append(normalized_seq.representative_name)
                        for normalized_seq_index in [modified_candidate_dict['normalized_seq_index']
                                                     for modified_candidate_dict in modified_candidate_dicts]:
                            normalized_seq_indices_assigned_to_modified_seqs.append(normalized_seq_index)

                for normalized_seq_index in sorted(set(normalized_seq_indices_assigned_to_modified_seqs), reverse=True):
                    sorted_normalized_seq_items.pop(normalized_seq_index)
                if len(sorted_normalized_seq_items) <= 2:
                    break

        modified_read_count = sum(t[1] for t in list(set(unique_seq_input_counts)))

        self.progress.end()

        self.run.info("Reads from modified tRNA", modified_read_count)


    def calc_normalization_stats(self):
        self.progress.new("Calculating normalized tRNA stats")

        # self.progress.update("Counting normalized sequences containing each trimmed sequence")
        normalized_count_dict = OrderedDict([(trimmed_seq.representative_name, 0) for trimmed_seq in self.trimmed_trna_seqs])
        for normalized_seq in self.normalized_trna_seqs:
            for trimmed_seq in normalized_seq.trimmed_seqs:
                normalized_count_dict[trimmed_seq.representative_name] += 1
        self.counts_of_normalized_seqs_containing_trimmed_seqs = [normalized_count for normalized_count in normalized_count_dict.values()]

        # self.progress.update("Finding the \"multiplicity\" of trimmed sequences among normalized sequences")
        multiplicity_dict = OrderedDict()
        for trimmed_seq, normalized_count_item in zip(self.trimmed_trna_seqs, normalized_count_dict.items()):
            trimmed_representative_name, normalized_count = normalized_count_item
            multiplicity_dict[trimmed_representative_name] = trimmed_seq.input_count * normalized_count
        self.multiplicities_of_trimmed_seqs_among_normalized_seqs = [multiplicity for multiplicity in multiplicity_dict.values()]

        # self.progress.update("Finding the \"average multiplicity\" of normalized sequences")
        for normalized_seq in self.normalized_trna_seqs:
            multiplicity_sum = 0
            for trimmed_seq in normalized_seq.trimmed_seqs:
                multiplicity_sum += multiplicity_dict[trimmed_seq.representative_name]
            self.average_multiplicities_of_normalized_seqs.append(round(multiplicity_sum / normalized_seq.input_count, 1))

        for modified_seq in self.modified_trna_seqs:
            multiplicity_sum = 0
            for normalized_seq in modified_seq.normalized_seqs:
                for trimmed_seq in normalized_seq.trimmed_seqs:
                    multiplicity_sum += multiplicity_dict[trimmed_seq.representative_name]
            self.average_multiplicities_of_modified_seqs.append(round(multiplicity_sum / modified_seq.input_count, 1))

        self.progress.end()


    def write_trimmed_table(self):
        self.progress.new("Writing tRNA-seq database table of trimmed tRNA sequences")
        trimmed_table_entries = []
        for trimmed_seq, normalized_seq_count in zip(self.trimmed_trna_seqs, self.counts_of_normalized_seqs_containing_trimmed_seqs):
            trimmed_table_entries.append(
                (trimmed_seq.representative_name,
                 len(trimmed_seq.unique_seqs),
                 trimmed_seq.input_count,
                 trimmed_seq.seq_string,
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
        for normalized_seq, average_multiplicity in zip(self.normalized_trna_seqs,
                                                        self.average_multiplicities_of_normalized_seqs):
            normalized_table_entries.append(
                (normalized_seq.representative_name,
                 len(normalized_seq.trimmed_seqs),
                 normalized_seq.input_count,
                 average_multiplicity,
                 normalized_seq.count_of_trimmed_seqs_mapped_to_threeprime_end,
                 normalized_seq.count_of_input_seqs_mapped_to_threeprime_end,
                 normalized_seq.count_of_trimmed_seqs_mapped_to_interior,
                 normalized_seq.count_of_input_seqs_mapped_to_interior,
                 normalized_seq.count_of_trimmed_seqs_mapped_to_fiveprime_end,
                 normalized_seq.count_of_input_seqs_mapped_to_fiveprime_end)
                + tuple(normalized_seq.input_acceptor_variant_count_dict.values()))

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


    def write_modified_table(self):
        self.progress.new("Writing tRNA-seq database table of modified tRNA sequences")

        modified_table_entries = []
        for modified_seq, average_multiplicity in zip(self.modified_trna_seqs,
                                                      self.average_multiplicities_of_modified_seqs):
            modified_table_entries.append(
                (modified_seq.representative_name,
                 ','.join([str(modification_index) for modification_index in modified_seq.modification_indices]),
                 modified_seq.seq_string,
                 ','.join([normalized_seq.representative_name for normalized_seq in modified_seq.normalized_seqs]),
                 len(modified_seq.normalized_seqs),
                 modified_seq.input_count,
                 average_multiplicity,
                 modified_seq.count_of_input_seqs_mapped_to_threeprime_end,
                 modified_seq.count_of_input_seqs_mapped_to_interior,
                 modified_seq.count_of_input_seqs_mapped_to_fiveprime_end)
                + tuple(modified_seq.input_acceptor_variant_count_dict.values()))

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db._exec_many(
            '''INSERT INTO %s VALUES (%s)'''
            % ('modified', ','.join('?' * len(tables.trnaseq_modified_table_structure))),
            modified_table_entries)

        modified_seq_count = len(self.modified_trna_seqs)
        trnaseq_db.db.set_meta_value('num_modified_trna_seqs', modified_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()
        self.run.info("Modified tRNA", modified_seq_count)


    def write_uniqued_nontrna_supplement(self):
        self.progress.new("Writing a file of sequences not identified as tRNA.")
        self.run.info("Output non-tRNA file", self.uniqued_nontrna_path)

        with open(self.uniqued_nontrna_path, 'w') as nontrna_file:
            nontrna_file.write("\t".join(self.UNIQUED_NONTRNA_HEADER) + "\n")
            for unique_seq in self.unique_nontrna_seqs:
                nontrna_file.write(unique_seq.representative_name + "\t"
                                   + ",".join(unique_seq.input_names) + "\t"
                                   + str(unique_seq.input_count) + "\t"
                                   + unique_seq.seq_string + "\n")

        self.progress.end()


    def write_uniqued_trna_supplement(self):
        self.progress.new("Writing a file showing which input sequences formed unique tRNA sequences")
        self.run.info("Output uniqued tRNA file", self.uniqued_trna_path)

        with open(self.uniqued_trna_path, 'w') as uniqued_file:
            uniqued_file.write("\t".join(self.UNIQUED_TRNA_HEADER) + "\n")
            for unique_seq in self.unique_trna_seqs:
                uniqued_file.write(unique_seq.representative_name + "\t"
                                   + ",".join(unique_seq.input_names) + "\n")

        self.progress.end()


    def write_trimmed_supplement(self):
        self.progress.new("Writing a file showing how trimmed tRNA sequences were formed from unique sequences")
        self.run.info("Output trimmed tRNA file", self.trimmed_ends_path)

        with open(self.trimmed_ends_path, 'w') as trimmed_file:
            trimmed_file.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for trimmed_seq in sorted(self.trimmed_trna_seqs,
                                      key=lambda trimmed_seq: -trimmed_seq.input_count):
                representative_name = trimmed_seq.representative_name
                for unique_seq in sorted(trimmed_seq.unique_seqs,
                                         key=lambda unique_seq: (-unique_seq.extra_fiveprime_length,
                                                                 -unique_seq.acceptor_length)):
                    trimmed_file.write(representative_name + "\t"
                                       + unique_seq.representative_name + "\t"
                                       + unique_seq.seq_string[: unique_seq.extra_fiveprime_length] + "\t"
                                       + unique_seq.seq_string[-unique_seq.acceptor_length: ] + "\t"
                                       + str(unique_seq.input_count) + "\n")

        self.progress.end()


    def write_normalized_supplement(self):
        self.progress.new("Writing a file showing how normalized tRNA sequences were formed from trimmed sequences")
        self.run.info("Output normalized tRNA file", self.normalized_fragments_path)

        with open(self.normalized_fragments_path, 'w') as normalized_file:
            normalized_file.write("\t".join(self.NORMALIZED_FRAGMENTS_HEADER) + "\n")
            for normalized_seq in sorted(self.normalized_trna_seqs,
                                         key=lambda normalized_seq: -normalized_seq.input_count):
                representative_name = normalized_seq.representative_name
                for trimmed_seq, start_position, end_position in sorted(
                    zip(normalized_seq.trimmed_seqs, normalized_seq.start_positions, normalized_seq.end_positions),
                    key=lambda t: (t[1], -t[2])):
                    normalized_file.write(representative_name + "\t"
                                          + trimmed_seq.representative_name + "\t"
                                          + str(start_position) + "\t"
                                          + str(end_position) + "\n")

        self.progress.end()


    def process(self):
        self.sanity_check()

        self.create_trnaseq_db()

        # Profile each input sequence for tRNA features.
        self.profile_trna(self.get_unique_input_seqs())

        # Trim 5' and 3' ends of profiled tRNA.
        self.trim_ends(self.unique_trna_seqs)

        # Map short 3' fragments and then trim their 3' ends.
        self.map_threeprime()

        # Sort the combined list of profiled and mapped trimmed sequences.
        self.trimmed_trna_seqs.sort(key=lambda trimmed_seq: trimmed_seq.representative_name)

        # Consolidate 3' fragments of longer tRNA sequences.
        self.dereplicate_threeprime()

        # Map fragments derived from the interior and 5' end of tRNA.
        self.map_nonthreeprime()

        self.find_modifications()

        self.calc_normalization_stats()
        self.write_trimmed_table()
        self.write_normalized_table()
        self.write_modified_table()

        self.write_uniqued_nontrna_supplement()
        self.write_uniqued_trna_supplement()
        self.write_trimmed_supplement()
        self.write_normalized_supplement()
