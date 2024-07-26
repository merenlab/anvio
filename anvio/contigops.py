# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes and functions for handling, storing, and retrieving atomic data
   from contigs and splits. Also includes classes to deal with external
   contig data such as GenbankToAnvio."""

import os
import re
import io
import gzip
import numpy as np
import string
import argparse

from Bio import SeqIO
from collections import OrderedDict

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.variability import VariablityTestFactory, ProcessNucleotideCounts, ProcessCodonCounts, ProcessIndelCounts


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = ["Mike Lee", "Faruk Uzun"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


OK_CHARS_FOR_ORGANISM_NAME = string.ascii_letters + string.digits + '_'
OK_CHARS_FOR_ACCESSION = OK_CHARS_FOR_ORGANISM_NAME

# These filter SNVs and INDELs, respectively, based on a coverage-dependent departure from reference value
variability_test_class_default = VariablityTestFactory(params={'b': 2, 'm': 1.45, 'c': 0.05})
indel_test_class_default = VariablityTestFactory(params={'b': 2, 'm': 1.45, 'c': 0.05})

# These are null filters, which do not filter SNVs and INDELs based on a coverage-dependent departure from reference value
variability_test_class_null = VariablityTestFactory(params=None)
indel_test_class_null = VariablityTestFactory(params=None)


def gen_split_name(parent_name, order):
    return '_'.join([parent_name, 'split', '%05d' % (order + 1)])


def get_atomic_data(sample_id, contigs, atomic_data_field):
    """Takes a list of contigops.Contig objects, and returns views for an atomic_data_field"""

    atomic_data_contigs = []
    atomic_data_splits = []

    # this loop will get atomic_data information from Contig instanes and store them into the db
    # at once. this was broken down into about 10 functions, but this structure seems to be the most efficient
    # although it looks crappy:
    for contig in contigs:
        contig_atomic_data = contig.get_atomic_data_dict(atomic_data_field)

        for split in contig.splits:
            atomic_data_contigs.append((split.name, sample_id, contig_atomic_data), )

        # contig is done, deal with splits in it:
        for split in contig.splits:
            split_atomic_data = split.get_atomic_data_dict(atomic_data_field)
            atomic_data_splits.append((split.name, sample_id, split_atomic_data), )

    return atomic_data_splits, atomic_data_contigs


class Contig:
    def __init__(self, name):
        self.name = name
        self.sequence = None
        self.parent = None
        self.splits = []
        self.length = 0
        self.abundance = 0.0
        self.coverage = anvio.bamops.Coverage()

        self.skip_SNV_profiling = False


    def get_atomic_data_dict(self, atomic_data_field):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'mean_coverage_Q2Q3': self.coverage.mean_Q2Q3,
             'detection': self.coverage.detection,
             'abundance': self.abundance,
             'variability': sum(s.auxiliary.variation_density for s in self.splits) if not self.skip_SNV_profiling else None}

        return d[atomic_data_field]


    def analyze_coverage(self, bam, min_percent_identity):

        if min_percent_identity is None:
            self.coverage.run(
                bam,
                self,
                read_iterator='fetch',
                max_coverage=anvio.auxiliarydataops.DEFAULT_COVERAGE_MAX_VALUE,
            )
        else:
            self.coverage.run(
                bam,
                self,
                read_iterator='fetch_filter_and_trim',
                max_coverage=anvio.auxiliarydataops.DEFAULT_COVERAGE_MAX_VALUE,
                percent_id_cutoff=min_percent_identity,
            )

        if len(self.splits) == 1:
            # Coverage.process_c is a potentially expensive operation (taking up ~90% of the time of
            # analyze_coverage when coverage is low (but when coverage is >500X, the time taken is
            # ~1%)). Regardless, this clause exists to catch the somewhat common occurence when a
            # contig only has one split. In this case, the contig _is_ the split, and so all of the
            # split.coverage attributes can simply be referenced directly from contig.coverage
            split = self.splits[0]
            split.coverage = anvio.bamops.Coverage()
            split.coverage.c = self.coverage.c
            split.coverage.min = self.coverage.min
            split.coverage.max = self.coverage.max
            split.coverage.median = self.coverage.median
            split.coverage.mean = self.coverage.mean
            split.coverage.std = self.coverage.std
            split.coverage.detection = self.coverage.detection
            split.coverage.is_outlier = self.coverage.is_outlier
            split.coverage.is_outlier_in_parent = self.coverage.is_outlier
            split.coverage.mean_Q2Q3 = self.coverage.mean_Q2Q3
            return

        for split in self.splits:
            split.coverage = anvio.bamops.Coverage()
            split.coverage.c = self.coverage.c[split.start:split.end]
            split.coverage.is_outlier_in_parent = self.coverage.is_outlier[split.start:split.end]
            split.coverage.process_c(split.coverage.c)


class Split:
    def __init__(self, name, sequence, parent, order, start=0, end=0):
        self.name = name
        self.sequence = sequence
        self.parent = parent
        self.end = end
        self.order = order
        self.start = start
        self.length = end - start
        self.explicit_length = 0
        self.abundance = 0.0
        self.auxiliary = None
        self.num_SNV_entries = 0
        self.num_INDEL_entries = 0
        self.num_SCV_entries = {}
        self.SNV_profiles = {}
        self.SCV_profiles = {}
        self.INDEL_profiles = {}
        self.per_position_info = {} # stores per nt info that is not coverage


    def get_atomic_data_dict(self, atomic_data_field):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'mean_coverage_Q2Q3': self.coverage.mean_Q2Q3,
             'detection': self.coverage.detection,
             'abundance': self.abundance,
             'variability': self.auxiliary.variation_density if self.auxiliary else None}

        return d[atomic_data_field]


class Auxiliary:
    def __init__(self, split, min_coverage_for_variability=10, report_variability_full=False,
                 profile_SCVs=False, skip_INDEL_profiling=False, skip_SNV_profiling=False,
                 min_percent_identity=None, skip_edges=0):

        if anvio.DEBUG:
            self.run = terminal.Run()

        self.split = split
        self.variation_density = 0.0
        self.min_coverage_for_variability = min_coverage_for_variability
        self.min_percent_identity = min_percent_identity
        self.skip_SNV_profiling = skip_SNV_profiling
        self.profile_SCVs = profile_SCVs
        self.skip_INDEL_profiling = skip_INDEL_profiling
        self.report_variability_full = report_variability_full
        self.skip_edges = skip_edges

        # used during array processing
        self.nt_to_array_index = {nt: i for i, nt in enumerate(constants.nucleotides)}
        self.cdn_to_array_index = {cdn: i for i, cdn in enumerate(constants.codons)}

        if self.profile_SCVs:
            if not all([necessary in self.split.per_position_info for necessary in ['forward', 'gene_start', 'gene_stop']]):
                raise ConfigError("Auxiliary :: self.split.per_position_info does not contain the info required for SCV profiling")


    def process(self, bam):
        self.run_SNVs_and_indels(bam)

        if self.profile_SCVs:
            self.run_SCVs(bam)


    def run_SCVs(self, bam):
        """Profile SCVs

        Parameters
        ==========
        bam : bamops.BAMFileObject

        Notes
        =====
        - Loop through reads and then finding the genes it overlaps with has proven to be the
          fastest way to do this. Looping through genes with fetch(self.split.parent, gene_start,
          gene_stop) is approximately twice as slow because trimming is expensive, and I/O
          operations on the BAM file suffer
        """

        reference_codon_sequences = {}
        gene_allele_counts = {}
        gene_calls = {}

        genes_with_SNVs = set(self.split.SNV_profiles['corresponding_gene_call'])

        # Decide how we want to iterate through reads
        if self.min_percent_identity:
            read_iterator = bam.fetch_filter_and_trim
            kwargs = {'percent_id_cutoff': self.min_percent_identity}
        else:
            read_iterator = bam.fetch_and_trim
            kwargs = {}

        read_count = 0
        for read in read_iterator(self.split.parent, self.split.start, self.split.end, **kwargs):
            # This loop will be making extensive use of the vectorized form of Read, so it is well
            # worth the time to vectorize it right off the bat
            read.vectorize()

            # Although rare, the read can overlap multiple genes. To generalize, we loop through the
            # genes the read overlaps with. If the read overlaps just one, shortcuts are taken
            # _within_ the loop.

            gene_id_per_nt_in_read = self.split.per_position_info['corresponding_gene_call'][
                (read.reference_start - self.split.start):(read.reference_end - self.split.start)
            ]

            genes_in_read = set(gene_id_per_nt_in_read)

            for gene_id in genes_in_read:

                if gene_id == -1:
                    continue

                if gene_id not in genes_with_SNVs:
                    # By design, we include SCVs only if they contain a SNV. In this case, there were no
                    # SNVs in the gene so there is nothing to do.
                    continue

                if len(genes_in_read) == 1:
                    # The read maps entirely to 1 gene. Easy peasy.
                    gene_overlap_start = read.reference_start
                    segment_that_overlaps_gene = read.v
                else:
                    # Okay, we need to do some work to get the segment that overlaps

                    # FIXME There is something extremely rare that can happen that leads to an
                    # inaccuracy. Here's the situation: A read completely covers a very small gene
                    # that is _fully_ inside another gene (such a situation can happen when running
                    # tRNA HMMs). Since the genes overlap, the nt positions in the small gene are
                    # given a value of -1, so gene_id_per_nt_in_read looks like [42, 42, 42, 42, -1,
                    # ..., -1, 42, 42]. The portion of the read after the consecutive -1's will be
                    # trimmed and not included. The solution to this problem is to better manage
                    # which genes we include for SCV analysis, but the problem is that the
                    # vectorized data we use for gene information (nt_position_info) has these
                    # overlapping gene segments pre-baked in.
                    gene_overlap_start, gene_overlap_stop = utils.get_constant_value_blocks(gene_id_per_nt_in_read, gene_id)[0]
                    gene_overlap_start += read.reference_start
                    gene_overlap_stop += read.reference_start - 1
                    start_index = utils.find_value_index(read[:, 0], gene_overlap_start)
                    stop_index = utils.find_value_index(read[:, 0], gene_overlap_stop)
                    segment_that_overlaps_gene = read[start_index:stop_index+1]

                if gene_id not in gene_calls:
                    # We make an on-the-fly gene call dict. See the NOTE in
                    # profiler.BAMProfiler.populate_gene_info_for_splits if you are confused by why
                    # we do not pass this information to split beforehand.
                    #
                    # We need to access gene-wide attributes from per-nt arrays, so any index in the
                    # array will suffice, so long as it corresponds to the gene_id. We arbitrarily
                    # pick the index corresponding to the genes starting position
                    accessor = gene_overlap_start - self.split.start
                    gene_calls[gene_id] = {
                        'contig': self.split.parent,
                        'start': self.split.per_position_info['gene_start'][accessor],
                        'stop': self.split.per_position_info['gene_stop'][accessor],
                        'direction': 'f' if self.split.per_position_info['forward'][accessor] else 'r',
                        'call_type': constants.gene_call_types['CODING'] if self.split.per_position_info['in_coding_gene_call'][accessor] else 2,
                    }

                gene_call = gene_calls[gene_id]

                if gene_call['call_type'] != constants.gene_call_types['CODING']:
                    # We cannot handle non-coding genes because they have no frame
                    continue

                for gapless_segment in read.iterate_blocks_by_mapping_type(mapping_type=0, array=segment_that_overlaps_gene):
                    block_start, block_end = gapless_segment[0, 0], gapless_segment[-1, 0] + 1

                    if block_end - block_start < 3:
                        # This block does not contain a full codon
                        continue

                    block_start_split = block_start - self.split.start
                    block_end_split = block_end - self.split.start

                    # this read must not contribute to codons it does not fully cover. Hence, we
                    # must determine by how many nts on each side we must trim
                    base_positions = self.split.per_position_info['base_pos_in_codon'][block_start_split:block_end_split]

                    first_pos = utils.find_value_index(base_positions, (1 if gene_call['direction'] == 'f' else 3))
                    last_pos = utils.find_value_index(base_positions, (3 if gene_call['direction'] == 'f' else 1), reverse_search=True)

                    if last_pos - first_pos < 3:
                        # the required trimming creates a sequence that is less than a codon long.
                        # We cannot use this read.
                        continue

                    # At this point, we are 100% sure this segment of the read will contribute to SCVs
                    gapless_segment = gapless_segment[first_pos:(last_pos+1), :]

                    # Update these for posterity
                    block_start_split += first_pos
                    block_end_split -= block_end - block_start - last_pos - 1

                    if gene_id not in gene_allele_counts:
                        # This is the first time a read has contributed to this gene_id, so we log
                        # its reference codon sequence and initialize an allele counts array
                        reference_codon_sequences[gene_id] = self.get_codon_sequence_for_gene(gene_call)
                        gene_allele_counts[gene_id] = self.init_allele_counts_array(gene_call)

                    codon_sequence_as_index = (
                        utils.nt_seq_to_codon_num_array(gapless_segment[:, 1], is_ord=True)
                        if gene_call['direction'] == 'f'
                        else utils.nt_seq_to_RC_codon_num_array(gapless_segment[:, 1], is_ord=True)
                    )

                    start, stop = self.split.per_position_info['codon_order_in_gene'][[block_start_split, block_end_split - 1]]
                    if gene_call['direction'] == 'r': start, stop = stop, start
                    codon_orders = np.arange(start, stop + 1)

                    # Codons with ambiguous characters have index values of 64. Remove them here
                    codon_orders = codon_orders[codon_sequence_as_index <= 63]
                    codon_sequence_as_index = codon_sequence_as_index[codon_sequence_as_index <= 63]

                    gene_allele_counts[gene_id] = utils.add_to_2D_numeric_array(
                        codon_sequence_as_index,
                        codon_orders,
                        gene_allele_counts[gene_id]
                    )

            read_count += 1

        if anvio.DEBUG: self.run.info_single('Done SCVs for %s (%d reads processed)' % (self.split.name, read_count), nl_before=0, nl_after=0)

        for gene_id in gene_allele_counts:
            allele_counts = gene_allele_counts[gene_id]

            if not np.sum(allele_counts):
                # There were no viable reads that landed on the gene
                continue

            cdn_profile = ProcessCodonCounts(
                allele_counts=allele_counts,
                allele_to_array_index=self.cdn_to_array_index,
                sequence=reference_codon_sequences[gene_id],
                min_coverage_for_variability=1,
            )

            # By design, we include SCVs only if they contain a SNV--filter out those that do not
            codon_orders_that_contain_SNVs = self.get_codon_orders_that_contain_SNVs(gene_id)
            if len(codon_orders_that_contain_SNVs):
                # Only keep those that contain SNVs
                cdn_profile.filter(codon_orders_that_contain_SNVs)
            else:
                # No codon orders contain SNVs
                continue

            if cdn_profile.process(skip_competing_items=True):
                # Processing allele counts yielded a non-zero number of codons being kept
                self.split.SCV_profiles[gene_id] = cdn_profile.d
                self.split.num_SCV_entries[gene_id] = len(cdn_profile.d['coverage'])
            else:
                # There were no codon positions worth keeping. We do not add this gene to
                # self.split.SCV_profiles
                pass

        if anvio.DEBUG: self.run.info_single('%d SCVs to report' % (sum(self.split.num_SCV_entries.values())), nl_before=0, nl_after=0, level=2)


    def get_codon_orders_that_contain_SNVs(self, gene_id):
        """Helper for run_SCVs"""

        matches_gene_boolean = self.split.SNV_profiles['corresponding_gene_call'] == gene_id

        return np.unique(self.split.SNV_profiles['codon_order_in_gene'][matches_gene_boolean])


    def get_codon_sequence_for_gene(self, gene_call):
        """Helper for run_SCVs"""

        seq_dict = {self.split.parent: {'sequence': self.split.sequence}}

        return utils.get_list_of_codons_for_gene_call(gene_call, seq_dict, subtract_by=self.split.start)


    def init_allele_counts_array(self, gene_call):
        """Create the array that will house the codon allele counts for a gene"""

        allele_counts_array_shape = (len(constants.codons), (gene_call['stop'] - gene_call['start']) // 3)

        return np.zeros(allele_counts_array_shape)


    def run_SNVs_and_indels(self, bam):
        """Profile SNVs (and indels if not self.skip_INDEL_profiling)

        Parameters
        ==========
        bam : bamops.BAMFileObject

        Notes
        =====
        See `variability-profile` artifact under anvio/docs/artifacts for details.
        """

        additional_per_position_data = self.split.per_position_info
        additional_per_position_data.update({
            'cov_outlier_in_split': self.split.coverage.is_outlier.astype(int),
            'cov_outlier_in_contig': self.split.coverage.is_outlier_in_parent.astype(int),
        })

        if not self.skip_INDEL_profiling:
            indels = {}
            get_indel_entry = lambda indel_type, seq, pos, length: OrderedDict([
                ('split_name', self.split.name),
                ('pos', pos),
                ('pos_in_contig', pos + self.split.start),
                ('corresponding_gene_call', additional_per_position_data['corresponding_gene_call'][pos]),
                ('in_noncoding_gene_call', additional_per_position_data['in_noncoding_gene_call'][pos]),
                ('in_coding_gene_call', additional_per_position_data['in_coding_gene_call'][pos]),
                ('base_pos_in_codon', additional_per_position_data['base_pos_in_codon'][pos]),
                ('codon_order_in_gene', additional_per_position_data['codon_order_in_gene'][pos]),
                ('cov_outlier_in_split', additional_per_position_data['cov_outlier_in_split'][pos]),
                ('cov_outlier_in_contig', additional_per_position_data['cov_outlier_in_contig'][pos]),
                ('reference', self.split.sequence[pos]),
                ('type', indel_type),
                ('sequence', seq),
                ('length', length),
                ('count', 1),
            ])

        # make an array with as many rows as there are nucleotides in the split, and as many rows as
        # there are nucleotide types. Each nucleotide (A, C, T, G, N) gets its own row which is
        # defined by the self.nt_to_array_index dictionary
        allele_counts_array_shape = (len(constants.nucleotides), self.split.length)
        allele_counts_array = np.zeros(allele_counts_array_shape)

        # Decide how we want to iterate through reads
        if self.min_percent_identity:
            read_iterator = bam.fetch_filter_and_trim
            kwargs = {'percent_id_cutoff': self.min_percent_identity}
        else:
            read_iterator = bam.fetch_and_trim
            kwargs = {}

        read_count = 0
        for read in read_iterator(self.split.parent, self.split.start, self.split.end, **kwargs):
            aligned_sequence_as_ord, reference_positions = read.get_aligned_sequence_and_reference_positions()

            # if the user is asking some nucleotides to be excluded from the calculation of
            # single-nucleotide variants due to DNA damage or other reasons, don't take them
            # into consideration:
            if self.skip_edges > 0:
                aligned_sequence_as_ord = aligned_sequence_as_ord[self.skip_edges:-self.skip_edges]
                reference_positions = reference_positions[self.skip_edges:-self.skip_edges]

            aligned_sequence_as_index = utils.nt_seq_to_nt_num_array(aligned_sequence_as_ord, is_ord=True)
            reference_positions_in_split = reference_positions - self.split.start

            allele_counts_array = utils.add_to_2D_numeric_array(aligned_sequence_as_index, reference_positions_in_split, allele_counts_array)

            if not self.skip_INDEL_profiling:
                read.vectorize()

                for ins_segment in read.iterate_blocks_by_mapping_type(mapping_type=1):
                    # Get the position and sequence of the insertion, create hash as a key for storage
                    ins_seq = ''.join([chr(x) for x in ins_segment[:, 1]])
                    ins_pos = ins_segment[0, 0] - self.split.start
                    indel_hash = hash((ins_pos, ins_seq))

                    if indel_hash in indels:
                        indels[indel_hash]['count'] += 1
                    else:
                        indels[indel_hash] = get_indel_entry(
                            indel_type='INS',
                            seq=ins_seq,
                            pos=ins_pos,
                            length=len(ins_seq),
                        )

                for del_segment in read.iterate_blocks_by_mapping_type(mapping_type=2):
                    # Get the position and length of the deletion, create hash as a key for storage
                    del_len = del_segment.shape[0]
                    del_pos = del_segment[0, 0] - self.split.start
                    indel_hash = hash((del_pos, del_len))

                    if indel_hash in indels:
                        indels[indel_hash]['count'] += 1
                    else:
                        indels[indel_hash] = get_indel_entry(
                            indel_type='DEL',
                            seq='',
                            pos=del_pos,
                            length=del_len,
                        )

            read_count += 1

        if anvio.DEBUG:
            self.run.info_single('Done SNVs for %s (%d reads processed)' % (self.split.name, read_count), nl_before=0, nl_after=0)

        split_as_index = utils.nt_seq_to_nt_num_array(self.split.sequence)
        nt_profile = ProcessNucleotideCounts(
            allele_counts=allele_counts_array,
            allele_to_array_index=self.nt_to_array_index,
            sequence=self.split.sequence,
            sequence_as_index=split_as_index,
            min_coverage_for_variability=self.min_coverage_for_variability,
            test_class=variability_test_class_null if self.report_variability_full else variability_test_class_default,
            additional_per_position_data=additional_per_position_data,
        )
        nt_profile.process()
        self.split.SNV_profiles = nt_profile.d

        if not self.skip_INDEL_profiling:
            indel_profile = ProcessIndelCounts(
                indels=indels,
                coverage=allele_counts_array.sum(axis=0),
                test_class=variability_test_class_null if self.report_variability_full else variability_test_class_default,
                min_coverage_for_variability=self.min_coverage_for_variability if not self.report_variability_full else 1,
            )
            indel_profile.process()
            self.split.INDEL_profiles = indel_profile.indels

        self.split.num_SNV_entries = len(nt_profile.d['coverage'])
        self.split.num_INDEL_entries = len(self.split.INDEL_profiles)
        self.variation_density = self.split.num_SNV_entries * 1000.0 / self.split.length

        if anvio.DEBUG:
            self.run.info_single('%d SNVs to report' % (self.split.num_SNV_entries), nl_before=0, nl_after=0, level=2)
            self.run.info_single('%d INDELs to report' % (self.split.num_INDEL_entries), nl_before=0, nl_after=0, level=2)


class GenbankToAnvioWrapper:
    """An ad hoc class that takes a metadata file generated by `ncbi-genome-download`"""

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args

        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.metadata_file_path = A('metadata')
        self.output_directory_path = os.path.abspath(A('output_dir') or os.path.curdir)
        self.output_fasta_descriptor = A('output_fasta_txt') or os.path.join(self.output_directory_path, 'fasta-input.txt')
        self.exclude_gene_calls_from_fasta_txt = A('exclude_gene_calls_from_fasta_txt')


    def sanity_check(self):
        filesnpaths.is_file_tab_delimited(self.metadata_file_path)

        if os.path.exists(self.output_directory_path):
            filesnpaths.is_output_dir_writable(self.output_directory_path)
        else:
            filesnpaths.gen_output_directory(self.output_directory_path)

        filesnpaths.is_output_file_writable(self.output_fasta_descriptor)


    def process(self):
        self.sanity_check()

        self.run.info('Input metadata file', self.metadata_file_path)
        self.run.info('Output directory', self.output_directory_path)

        columns = utils.get_columns_of_TAB_delim_file(self.metadata_file_path)
        if 'organism_name' not in columns or 'local_filename' not in columns:
            raise ConfigError("The metadata file you provided does not look like a metadata "
                              "file output from the program `ncbi-genome-download` :/ Why? "
                              "Because anvi'o expects that file to have at least the following "
                              "two columns in it: 'organism_name' and 'local_filename'.")

        metadata = utils.get_TAB_delimited_file_as_dictionary(self.metadata_file_path)

        for entry in metadata:
            if not os.path.exists(metadata[entry]['local_filename']):
                raise ConfigError("At least one of the files in your metadata input does not seem to be "
                                  "where they think they are :/ Please make sure the entry %s and others "
                                  "point to proper local file paths..." % entry)

        self.run.info('Num entries in metadata', len(metadata))

        output_fasta_dict = {}
        self.progress.new("GenBank to anvi'o", progress_total_items=len(metadata))
        for entry in metadata:
            self.progress.increment()
            self.progress.update('Processing %s ...' % entry)

            # set the organism name and accession id and clean them from weird
            # characters.
            organism_name = metadata[entry]['organism_name']
            for char in [c for c in organism_name if c not in OK_CHARS_FOR_ORGANISM_NAME]:
                organism_name = organism_name.replace(char, '_')

            accession_id = entry
            for char in [c for c in accession_id if c not in OK_CHARS_FOR_ACCESSION]:
                accession_id = accession_id.replace(char, '_')

            final_name = '_'.join([organism_name, accession_id])

            args = argparse.Namespace(input_genbank=metadata[entry]['local_filename'],
                                      output_file_prefix=os.path.join(self.output_directory_path, final_name))
            g = GenbankToAnvio(args, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))

            if final_name in output_fasta_dict:
                raise ConfigError("The final name '%s' for your genome has alrady been used by "
                                  "another one :/ This should never happen unless your metadata "
                                  "contains entries with identical accession numbers...")
            output_fasta_dict[final_name] = g.process()

        self.progress.end()

        headers = ['name', 'path']
        if not self.exclude_gene_calls_from_fasta_txt:
            headers.extend(['external_gene_calls', 'gene_functional_annotation'])

        utils.store_dict_as_TAB_delimited_file(output_fasta_dict, self.output_fasta_descriptor, headers=headers)

        self.run.info('Output FASTA descriptor', self.output_fasta_descriptor)


class GenbankToAnvio:
    """A class to deal with GenBank files. The initial code was implemented by Mike Lee"""

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args

        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_genbank_path = A('input_genbank')
        self.output_file_prefix = A('output_file_prefix')
        self.output_fasta_path = A('output_fasta')
        self.output_functions_path = A('output_functions')
        self.output_gene_calls_path = A('output_gene_calls')
        self.source = A('annotation_source') or 'NCBI_PGAP'
        self.version = A('annotation_version') or 'v4.6'
        self.omit_aa_sequences_column = A('omit_aa_sequences_column') or False
        self.include_locus_tags_as_functions = A('include_locus_tags_as_functions')

        # gene callers id start from 0. you can change your instance
        # prior to processing the genbank file to start from another
        # value
        self.gene_callers_id = 0

        # dumping gene if noted as these in the "note" section of the call
        self.note_terms_to_exclude = ["frameshifted", "internal stop", "incomplete"]

        # dumping gene if "location" section contains any of these: "join" means the
        # gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig
        # FIXME join does not necessarily mean the gene call spans multiple columns. See the first
        # gene here to see an instance where this is not true: https://www.ncbi.nlm.nih.gov/nuccore/MN908947
        self.location_terms_to_exclude = ["join", "<", ">"]


    def sanity_check(self):
        if self.output_file_prefix and (self.output_fasta_path or self.output_functions_path or self.output_gene_calls_path):
            raise ConfigError("Your arguments contain an output file prefix, and other output file paths. You can either "
                              "define a prefix, and the output files would be named accordingly (such as 'PREFIX-extenral-gene-calls', "
                              "'PREFIX-external-functions.txt', and 'PREFIX-contigs.fa'), ORRR you can set output file names "
                              "or paths for each of these files independently. You can also leave it as is for default file names to "
                              "be used. But you can't mix everything together and confuse us here.")

        self.output_fasta_path = self.output_fasta_path or 'contigs.fa'
        self.output_functions_path = self.output_functions_path or 'external-functions.txt'
        self.output_gene_calls_path = self.output_gene_calls_path or 'external-gene-calls.txt'

        if self.output_file_prefix:
            J = lambda x: '-'.join([self.output_file_prefix, x])
            self.output_fasta_path = J(self.output_fasta_path)
            self.output_functions_path = J(self.output_functions_path)
            self.output_gene_calls_path = J(self.output_gene_calls_path)

        filesnpaths.is_output_file_writable(self.output_fasta_path)
        filesnpaths.is_output_file_writable(self.output_functions_path)
        filesnpaths.is_output_file_writable(self.output_gene_calls_path)
        filesnpaths.is_file_exists(self.input_genbank_path)

        files_already_exist = [f for f in [self.output_fasta_path, self.output_functions_path, self.output_gene_calls_path] if os.path.exists(f)]
        if len(files_already_exist):
            raise ConfigError("Some of the output files already exist :/ Anvi'o feels uneasy about simply overwriting "
                              "them and would like to outsource that risk to you. Please either use different output "
                              "file names, or delete these files and come back: '%s'" % (', '.join(files_already_exist)))


    def get_genbank_file_object(self):
        try:
            if self.input_genbank_path.endswith('.gz'):
                genbank_file_object = SeqIO.parse(io.TextIOWrapper(gzip.open(self.input_genbank_path, 'r')), "genbank")
            else:
                genbank_file_object = SeqIO.parse(open(self.input_genbank_path, "r"), "genbank")
        except Exception as e:
            raise ConfigError("Someone didn't like your unput 'genbank' file :/ Here's what they said "
                              "about it: '%s'." % e)

        return genbank_file_object


    def process(self):
        self.sanity_check()

        output_fasta = {}
        output_gene_calls = {}
        output_functions = []
        num_genbank_records_processed = 0
        num_genes_found = 0
        num_genes_reported = 0
        num_genes_with_AA_sequences = 0
        num_partial_genes = 0
        num_genes_with_functions = 0

        num_genes_excluded = 0
        genes_excluded = set([])

        # A very quick look at the genbank file to see if translated sequences are present in it
        aa_sequences_present = False
        for genbank_record in self.get_genbank_file_object():
            genes = [gene for gene in genbank_record.features if gene.type =="CDS"]

            # clearly, not every genebank record has to have genes.
            # just like not every bank has to have money I guess.
            # IRONIES OF LIFE BUT IN CODING SPACE BECAUSE BIOINFORMATICS.
            if not len(genes):
                continue

            # do we have AA sequences in this?
            if 'translation' in genes[0].qualifiers:
                aa_sequences_present = True

            if aa_sequences_present and self.omit_aa_sequences_column:
                aa_sequences_present = False
                self.run.info_single("Amino acid sequences seem to be present in this GFF file, but you wanted anvi'o "
                                     "to not report them. FINE. They shall be ignored.", nl_after=1)
                break

        # The main loop to go through all records forreals.
        for genbank_record in self.get_genbank_file_object():
            if genbank_record.name in output_fasta:
                raise ConfigError("Anvi'o is not able to convert this GenBank file because it contains sequences with identical "
                                   "locus names :/. An example is locus '%s'." % genbank_record.name)
            else:
                num_genbank_records_processed += 1
                output_fasta[genbank_record.name] = str(genbank_record.seq)

            genes = [gene for gene in genbank_record.features if gene.type =="CDS"] # focusing on features annotated as "CDS" by NCBI's PGAP

            for gene in genes:
                num_genes_found += 1
                location = str(gene.location)
                # "join" in `location` means that the gene call spans multiple contigs
                # so we exclude it here:
                if 'join' in location:
                    num_genes_excluded += 1
                    if 'protein_id' in gene.qualifiers:
                        genes_excluded.add(gene.qualifiers['protein_id'][0])

                    continue

                if "note" in gene.qualifiers:
                    note = str(gene.qualifiers["note"][0])

                    # dumping gene if noted as any of these in the "note" section set above
                    if any(exclusion_term in note for exclusion_term in self.note_terms_to_exclude):
                        continue

                # dumping if overlapping translation frame
                if "transl_except" in gene.qualifiers:
                    continue

                # dumping if gene declared a pseudogene
                if "pseudo" in gene.qualifiers or "pseudogene" in gene.qualifiers:
                    continue

                # The character "<" or ">" in `location `means the gene call runs off a contig, so it
                # is best to mark it as impartial:
                if ('<' in location or '>' in location):
                    partial = 1
                    num_partial_genes += 1
                else:
                    partial = 0

                # cleaning up gene coordinates to more easily parse:
                location = location.replace("[", "").replace('>', '').replace('<', '')
                location = re.sub('](.*)', '', location)
                location = location.split(":")

                start = location[0] # start coordinate
                end = location[1] # end coordinate

                # setting direction to "f" or "r":
                if gene.strand == 1:
                    direction="f"
                else:
                    direction="r"

                # for accession, storing protein id if it has one, else the locus tag, else "None"
                if "protein_id" in gene.qualifiers:
                    accession = gene.qualifiers["protein_id"][0]
                elif "locus_tag" in gene.qualifiers:
                    accession = gene.qualifiers["locus_tag"][0]
                else:
                    accession = "None"

                # storing gene product annotation if present
                if "product" in gene.qualifiers:
                    function = gene.qualifiers["product"][0]
                    # trying to capture all different ways proteins are listed as hypothetical and
                    # setting to same thing so can prevent from adding to output functions table below
                    if function in ["hypothetical", "hypothetical protein", "conserved hypothetical",
                                    "conserved hypotheticals", "Conserved hypothetical protein"]:
                        function = None
                else:
                    function = None

                # if present, adding gene name to product annotation (so long as not a hypothetical,
                # sometimes these names are useful, sometimes they are not):
                if "gene" in gene.qualifiers:
                    if function:
                        gene_name=str(gene.qualifiers["gene"][0])
                        function = function + " (" + gene_name + ")"

                output_gene_calls[self.gene_callers_id] = {'contig': genbank_record.name,
                                                           'start': start,
                                                           'stop': end,
                                                           'direction': direction,
                                                           'partial': partial,
                                                           'call_type': 1,
                                                           'source': self.source,
                                                           'version': self.version}

                # let's keep the amino acid sequences if present
                if aa_sequences_present:
                    if 'translation' in gene.qualifiers and not partial:
                        output_gene_calls[self.gene_callers_id]['aa_sequence'] = gene.qualifiers["translation"][0]
                        num_genes_with_AA_sequences += 1
                    else:
                        output_gene_calls[self.gene_callers_id]['aa_sequence'] = ''

                num_genes_reported += 1

                # not writing gene out to functions table if no annotation
                if function:
                    output_functions.append({'gene_callers_id': self.gene_callers_id,
                                             'source': self.source,
                                             'accession': accession,
                                             'function': function,
                                             'e_value': 0})
                    num_genes_with_functions += 1

                if self.include_locus_tags_as_functions and "locus_tag" in gene.qualifiers:
                    locus_tag = gene.qualifiers["locus_tag"][0]
                    output_functions.append({'gene_callers_id': self.gene_callers_id,
                                             'source': 'GENBANK_LOCUS_TAG',
                                             'accession': locus_tag,
                                             'function': locus_tag,
                                             'e_value': 0})

                # increment the gene callers id for the next
                self.gene_callers_id += 1

        if num_genbank_records_processed == 0:
            raise ConfigError("It seems there was no records in your input genbank file :/ Are you sure you "
                              "gave the right file path that actually resolves to a genbank formatted "
                              "text file?")

        self.run.info('Num GenBank entries processed', num_genbank_records_processed)
        self.run.info('Num gene records found', num_genes_found)
        self.run.info('Num genes reported', num_genes_reported, mc='green')
        self.run.info('Num genes with AA sequences', num_genes_with_AA_sequences, mc='green')
        self.run.info('Num genes with functions', num_genes_with_functions, mc='green')
        self.run.info('Locus tags included in functions output?', "Yes" if self.include_locus_tags_as_functions else "No", mc='green')
        self.run.info('Num partial genes', num_partial_genes, mc='cyan')
        self.run.info('Num genes excluded', num_genes_excluded, mc='red', nl_after=1)

        if num_genes_excluded:
            msg = (f"A total of {num_genes_excluded} in this file were excluded during the processing of "
                   f"the GFF file because they were marked as 'spanning multiple contigs'.")
                # dumping gene if "location" section contains any of these terms set above: "join" means
                # the gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig
            if len(genes_excluded):
                self.run.warning(msg + "Here are some protein IDs from the file that were associated with "
                                 f"them: {', '.join(genes_excluded)}. But please note that such genes that "
                                 f"span multiple contigs are often associated with multiple 'protein IDs' in "
                                 f"the GFF description, and this list includes only the first protein ID if "
                                 f"there are multiple. In addition, not every gene that spans across multiple "
                                 f"contigs will have a protein ID, so you will not see them in this list :/ "
                                 f"If you think anvi'o should try to include them in its report anyway, let us "
                                 f"know.", header="WEIRD GENE CALLS EXCLUDED")
            else:
                self.run.warning(msg, header="WERID GENE CALLS EXCLUDED")

        if aa_sequences_present and num_partial_genes > 1:
            self.run.warning(f"Anvi'o found out that aminio acid seqeunces were present in the GFF file, so it "
                             f"reported them in the external gene calls file. But the {num_partial_genes} that "
                             f"were marked as partial in the GFF file, there will not be any amino acid sequences. "
                             f"Alternatively, you could use the flag `--omit-aa-seqeunces-column` to instruct "
                             f"this program to report an external gene calls file WITHOUT an amino acid sequences, "
                             f"column. When it is missing, anvi'o would do its best to translate your gene calls "
                             f"and may be able to find out what to do with your partial gene calls. But it usually "
                             f"is a MUCH better idea to learn translated sequences from a GFF file since they may "
                             f"include careful consideration of organism-specific genetic code.",
                             header="JUSTICE FOR PARTIAL GENES: NOT FOUND")

        # time to write these down:
        utils.store_dict_as_FASTA_file(output_fasta,
                                       self.output_fasta_path,
                                       wrap_from=None)
        self.run.info('FASTA file path', self.output_fasta_path)

        if len(output_gene_calls):
            header_for_external_gene_calls = ["gene_callers_id", "contig", "start", "stop", "direction", "partial", "call_type", "source", "version"]
            header_for_functions = ['gene_callers_id', 'source', 'accession', 'function', 'e_value']

            if aa_sequences_present:
                header_for_external_gene_calls.append('aa_sequence')

            utils.store_dict_as_TAB_delimited_file(output_gene_calls,
                                                   self.output_gene_calls_path,
                                                   headers=header_for_external_gene_calls)

            with open(self.output_functions_path, 'w') as output_functions_file:
                header_text = '\t'.join(header_for_functions)
                output_functions_file.write(f"{header_text}\n")
                for entry in output_functions:
                    entry_text = '\t'.join([f"{entry[k]}" for k in header_for_functions])
                    output_functions_file.write(f"{entry_text}\n")

            self.run.info('External gene calls file', self.output_gene_calls_path)
            self.run.info('TAB-delimited functions', self.output_functions_path)
        else:
            self.output_gene_calls_path = None
            self.output_functions_path = None
            self.run.warning("Anvi'o couldn't find any gene calles in the GenBank file, hence you will get "
                             "no output files for external gene calls or functions :/ We hope you can "
                             "survive this terrible terrible news :(")

        self.run.info_single('Mmmmm ☘ ', nl_before=1, nl_after=1)

        return {'external_gene_calls': self.output_gene_calls_path,
                'gene_functional_annotation': self.output_functions_path,
                'path': self.output_fasta_path}

