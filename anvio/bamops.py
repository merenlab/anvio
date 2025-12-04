# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes for anything BAM-related"""


import os
import sys
import pysam
import numpy as np
import hashlib

from numba import jit
from collections import Counter
from dataclasses import dataclass, field

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class PairedEndTemplateDist:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bam_file_path = A('bam_file')
        self.min_frequency_of_tlen = A('min_template_length_frequency') or 10
        self.max_tlen_to_consider = A('max_template_length_to_consider') or 500000
        self.output_file_path = A('output_file') or 'TEMPLATE-LENGTH-STATS.txt'
        self.plot_data = A('plot_data')

        self.run.info('BAM file', self.bam_file_path)

        self.template_lengths = {}


    def process(self):
        filesnpaths.is_file_bam_file(self.bam_file_path)
        filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=False)

        bam = BAMFileObject(self.bam_file_path, 'rb')
        self.contig_names = bam.references
        self.contig_lengths = bam.lengths

        self.progress.new('Initializing')
        self.progress.update('Counting num contigs ...')
        num_contigs = len(self.contig_names)
        self.progress.update('Counting num reads ...')
        num_reads = bam.count()
        self.progress.end()

        self.run.info('Number of contigs', pp(num_contigs))
        self.run.info('Number of reads', pp(num_reads))
        self.run.info('Minimum template length frequency', pp(self.min_frequency_of_tlen))
        self.run.info('Maximum template length to consider', pp(self.max_tlen_to_consider), nl_after=1)

        self.progress.new("Bleep bloop", progress_total_items=num_reads)
        self.progress.update('...')

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        for j in range(0, num_contigs):
            contig_name = self.contig_names[j]

            template_lengths_counter = Counter()

            num_reads_processed = 0
            for read in bam.fetch_only(contig_name):
                num_reads_processed += 1

                if read.is_paired and not read.mate_is_unmapped and read.tlen != 0:
                    template_lengths_counter[abs(read.tlen)] += 1

                if num_reads_processed % 100000 == 0:
                    if mem_tracker.measure():
                        mem_usage = mem_tracker.get_last()
                        mem_diff = mem_tracker.get_last_diff()

                    self.progress.increment(increment_to=num_reads_processed)
                    self.progress.update(f"READ {pp(num_reads_processed)}/{pp(num_reads)} :: CONTIG {pp(j+1)}/{pp(num_contigs)} :: MEMORY 🧠 {mem_usage} ({mem_diff})")

            # remove rarely occuring template lengths
            tlens_to_discard = [tlen for tlen in template_lengths_counter if template_lengths_counter[tlen] < self.min_frequency_of_tlen or tlen > self.max_tlen_to_consider]
            [template_lengths_counter.pop(tlen) for tlen in tlens_to_discard]

            # code below looks dumb, b/c we first count frequencies using a Counter,
            # and then turn it into a huge Python array by adding each tlen as many
            # times as it was observed. but in fact this is faster than `.append()`
            # for an array, or adding numbers to a numpy array. LOOK I HAM AS CONFUSE
            # AS U R, BUT I HAS TESTED IT.
            self.template_lengths[contig_name] = []
            for tlen in template_lengths_counter:
                self.template_lengths[contig_name].extend([tlen] * template_lengths_counter[tlen])

        self.progress.end()

        contigs_without_tlen_data = [c for c in self.template_lengths if not len(self.template_lengths[c])]
        if len(contigs_without_tlen_data) == num_contigs:
            raise ConfigError("So. None of your contigs seem to have any template length data, which suggests that "
                              "either the reads that were mapped to your sequences to generate this BAM file were "
                              "single reads rather than paired-end reads, or you set a crazy high `--min-tlen-frequency` "
                              "value even though you didn't have a large number of reads in your BAM file. WHICH ONE "
                              "IS IT, FLORIAN, WHICH ONE?")
        elif len(contigs_without_tlen_data) < num_contigs:
            self.run.warning(f"Some of your contigs, {pp(len(contigs_without_tlen_data))} of {pp(num_contigs)} to "
                             f"be precise, did not seem to have any template lenght data. There are many reasons "
                             f"this could happen, including a very high `--min-tlen-frequency` parameter for BAM files "
                             f"with small number of reads. But since there are some contigs that seem to have proper paired-"
                             f"end reads with template lengths, anvi'o will continue reporting and put zeros for "
                             f"contigs that have no data in output files.")
        else:
            # we're golden
            pass

        if self.output_file_path:
            self.report()

        if self.plot_data:
            self.plot_histogram()


    def report(self):
        stats = {}

        self.progress.new("Summarizing data for reporting", progress_total_items=len(self.contig_names))
        for contig_name in self.template_lengths:
            self.progress.update(f"{contig_name}", increment=True)

            if len(self.template_lengths[contig_name]):
                # the following class was designed for coverage stats, but it does a large part
                # of what we need to learn for these np arrays, so why not:
                arr = np.array(self.template_lengths[contig_name])
                C = utils.CoverageStats(arr, skip_outliers=True)

                stats[contig_name] = {'length': self.contig_lengths[self.contig_names.index(contig_name)],
                                      'num_reads_considered': len(arr),
                                      'mean'     : f"{C.mean:.4}",
                                      'mean_Q2Q3': f"{C.mean_Q2Q3:.4}",
                                      'median'   : f"{C.median}",
                                      'min'      : f"{C.min}",
                                      'max'      : f"{C.max}",
                                      'std'      : f"{C.std:.4}"}
            else:
                stats[contig_name] = {'length': self.contig_lengths[self.contig_names.index(contig_name)],
                                      'num_reads_considered': 0,
                                      'mean'     : "",
                                      'mean_Q2Q3': "",
                                      'median'   : "",
                                      'min'      : "",
                                      'max'      : "",
                                      'std'      : ""}

        self.progress.end()

        headers = ['contig', 'length', 'num_reads_considered', 'mean', 'mean_Q2Q3', 'median', 'min', 'max', 'std']
        utils.store_dict_as_TAB_delimited_file(stats, self.output_file_path, headers=headers)

        self.run.info('Output file', self.output_file_path, nl_after=1)


    def plot_histogram(self, num_bins=50):
        try:
            import plotext as plt
        except:
            self.run.warning("You don't have the `plotext` library to plot data :/ You can "
                             "install it by running `pip install plotext` in your anvi'o "
                             "environment.", header="NO PLOT FOR YOU :(")
            return

        self.progress.new("Bleep bloop")
        self.progress.update('Merging all tlens from all contigs ...')
        tlens_across_all_contigs = []
        for contig_name in self.template_lengths:
            tlens_across_all_contigs.extend(self.template_lengths[contig_name])
        self.progress.end()

        self.run.info_single("All the data are ready for plotting, if you are lucky, you "
                             "should see a histogram of tlen distributions below.", nl_after=1)

        plt.clc()
        plt.title(f"The overall distribution of pair1 and pair2 distances in {os.path.basename(self.bam_file_path)}")
        plt.xlabel("Template length")
        plt.ylabel("Frequency")
        plt.hist(tlens_across_all_contigs, num_bins)
        plt.plotsize(progress.terminal_width, 40)
        plt.show()


class BAMFileObject(pysam.AlignmentFile):
    def __init__(self, *args):
        """A class that is essentially pysam.AlignmentFile, with some added bonuses

        This class inherits pysam.AlignmentFile and adds a little flair. Init such an object the
        way you would an AlignmentFile, i.e. bam = bamops.BAMFileObject(path_to_bam)
        """

        self.input_bam_path = args[0]

        filesnpaths.is_file_bam_file(self.input_bam_path)

        self.fetch_filter = None

        pysam.AlignmentFile.__init__(self)


    def fetch_only(self, contig_name, start=None, end=None, *args, **kwargs):
        """A wrapper function for `bam.fetch()`.

        Having a single function enables global application of `fetch_filters` to reads
        reported from a BAM file.
        """

        for read in self.fetch(contig_name, start, end):
            if self.fetch_filter:
                if constants.fetch_filters[self.fetch_filter](read):
                    yield read
            else:
                yield read


    def fetch_and_trim(self, contig_name, start, end, *args, **kwargs):
        """Returns an read iterator that trims overhanging reads

        Like pysam.AlignmentFile.fetch(), except trims reads that overhang the start and end of the
        defined region so that they fit inside the start and stop.
        """

        for read in self.fetch_only(contig_name, start, end, *args, **kwargs):
            if read.is_unmapped or read.cigartuples is None or read.query_sequence is None:
                # This read either has no associated cigar string or no query sequence. If cigar
                # string is None, this means it did not align but is in the BAM file anyways, or the
                # mapping software decided it did not want to include a cigar string for this read.
                # The origin of why query sequence would be None is unkown. Alternatively, this read
                # has a query sequence, has a cigar string (and therefore mapped), but for some
                # reason pysam has given the read the attribute read.is_unmapped == True. We choose
                # not to work with such a read (see https://github.com/merenlab/anvio/issues/1821)
                continue

            read = Read(read)

            if start - read.reference_start > 0:
                read.trim(trim_by=(start - read.reference_start), side='left')

                if read.reference_start >= end:
                    # You may find this line unnecessary. However in extreme fringe cases, this
                    # condition *can* be met, please see https://github.com/merenlab/anvio/issues/1436
                    continue

            if read.reference_end - end > 0:
                read.trim(trim_by=(read.reference_end - end), side='right')

                if read.reference_end <= start:
                    # You may find this line unnecessary. However in extreme fringe cases, this
                    # condition *can* be met, please see https://github.com/merenlab/anvio/issues/1436
                    continue

            yield read


    def fetch_filter_and_trim(self, contig_name, start, end, percent_id_cutoff=0, *args, **kwargs):
        """Returns an read iterator that trims and also filters reads based on their characteristics

        Trims just as self.fetch_and_trim does. Filtering is applied before trimming.

        Parameters
        ==========
        start : int
            Pass None if you don't want a start position cutoff

        end : int
            Pass None if you don't want an end cutoff

        percent_id_cutoff : float, 0
            Reads with a percentage ID relative to their mapped segments of the reference genome
            less than this amount will be excluded. E.g. if you want 95% as a cutoff, use
            percent_id_cutoff=95.0
        """

        for read in self.fetch_only(contig_name, start, end, *args, **kwargs):
            if read.is_unmapped or read.cigartuples is None or read.query_sequence is None:
                # This read either has no associated cigar string or no query sequence. If cigar
                # string is None, this means it did not align but is in the BAM file anyways, or the
                # mapping software decided it did not want to include a cigar string for this read.
                # The origin of why query sequence would be None is unkown. Alternatively, this read
                # has a query sequence, has a cigar string (and therefore mapped), but for some
                # reason pysam has given the read the attribute read.is_unmapped == True. We choose
                # not to work with such a read (see https://github.com/merenlab/anvio/issues/1821)
                continue

            read = Read(read)

            # The read must be vectorized before calling read.get_percent_id
            read.vectorize()

            if read.get_percent_id() < percent_id_cutoff:
                continue

            # Trimming takes place _after_ filtering based on percent ID

            if start - read.reference_start > 0:
                read.trim(trim_by=(start - read.reference_start), side='left')

                if read.reference_start >= end:
                    # You may find this line unnecessary. However in extreme fringe cases, this
                    # condition *can* be met, please see https://github.com/merenlab/anvio/issues/1436
                    continue

            if read.reference_end - end > 0:
                read.trim(trim_by=(read.reference_end - end), side='right')

                if read.reference_end <= start:
                    # You may find this line unnecessary. However in extreme fringe cases, this
                    # condition *can* be met, please see https://github.com/merenlab/anvio/issues/1436
                    continue

            yield read


    def reads_have_MD_tags(self, num_reads=10000, require=all):
        """A holistic approach to testing if reads in the BAM file have MD tags

        Tests the first `num_reads` in the file to see if they contain the MD tag.  The MD tag
        provides information about the reference bases in the aligned segment and are essential for
        some read operations, such as calculating percent ID of a read to its reference.

        Parameters
        ==========
        num_reads : int, 10000
            How many reads should be sampled?

        require : function, all
            function that takes list of bool values and returns single bool. E.g. built-ins all or
            any. The default is all, which means all of the sampled reads are required to have the
            MD tag

        Returns
        =======
        output : bool
            require(has_MD), where has_MD is a list of booleans (length has_MD) that are True if the
            read has the MD tag and False otherwise
        """

        has_MD = []
        reads_checked = 0

        for contig_name in self.references:
            for read in self.fetch_only(contig_name):
                if read.cigartuples is None:
                    # if it didn't map we forgive it for not having an MD tag
                    continue

                has_MD.append(read.has_tag('MD'))

                reads_checked += 1
                if reads_checked == num_reads:
                    return require(has_MD)

        if reads_checked == 0:
            return True

        return require(has_MD)


class Read:
    def __init__(self, read):
        """Class for manipulating reads

        Parameters
        ==========
        read : pysam.AlignedSegment
        """

        # redefine all properties of interest explicitly from pysam.AlignedSegment object as
        # attributes of this class. The reason for this is that some of the AlignedSegment
        # attributes have no __set__ methods, so are read only. Since this class is designed to
        # modify some of these attributes, and since we want to maintain consistency across
        # attributes, all attributes of interest are redefined here
        self.cigartuples = np.array(read.cigartuples)
        self.query_sequence = np.frombuffer(read.query_sequence.encode('ascii'), np.uint8)
        self.reference_start = read.reference_start
        self.reference_end = read.reference_end

        if read.has_tag('MD'):
            self.reference_sequence = np.frombuffer(read.get_reference_sequence().upper().encode('ascii'), np.uint8)
        else:
            # Calculate reference sequence length excluding N operations (skipped regions)
            ref_seq_length = self.reference_end - self.reference_start
            for op, length in read.cigartuples:
                if op == 3:  # N operation (skipped region)
                    ref_seq_length -= length
            self.reference_sequence = np.array([ord('N')] * ref_seq_length)

        # See self.vectorize
        self.v = None


    def vectorize(self):
        """Set the self.v attribute to provide array-like access to the read

        0th column = reference position
        1st column = query sequence
        2nd column = mapping type (0=mapped, 1=read insertion, or 2=read deletion)
        3rd column = reference sequence
        """

        self.v = _vectorize_read(
            self.cigartuples,
            self.query_sequence,
            self.reference_sequence,
            self.reference_start,
            constants.cigar_consumption
        )


    def iterate_blocks_by_mapping_type(self, mapping_type, array=None):
        """Iterate through slices of array that contain blocks of a given mapping type

        Parameters
        ==========
        mapping_type : int
            Any of 0, 1, 2, or -1. 0 = mapping segment, 1 = read insertion segment, 2 = read
            deletion segment, -1 = gap in read and reference

        array : numpy array, None
            If None, self.v will be used

        Yields
        ======
        output : numpy arrays
            Each numpy array corresponds to a section of self.v that contained consecutive
            mapping_types.
        """
        if array is None:
            array = self.v

        for start, stop in utils.get_constant_value_blocks(array[:, 2], mapping_type):
            yield array[start:stop, :]


    def __getitem__(self, key):
        """Used to access the vectorized form of the read, self.v"""

        return self.v.__getitem__(key)


    def get_percent_id(self):
        """Get the percent id of the read to the reference

        Returns
        =======
        output : float
            The percent ID. E.g. 95.2

        Notes
        =====
        - Only mapped segments (i.e. consumed by read and ref) are considered. E.g if 50% of a read
          maps with 100% identity and the other 50% is composed of indels, this function returns 100.0
        - Delegates to the just-in-time compiled function _get_percent_id
        """

        return _get_percent_id(self.v)


    def get_aligned_sequence_and_reference_positions(self):
        """Get the aligned sequence at each mapped position, and the positions themselves

        Notes
        =====
        - Delegates to the just-in-time compiled function
          _get_aligned_sequence_and_reference_positions
        """

        return _get_aligned_sequence_and_reference_positions(
            self.cigartuples,
            self.query_sequence,
            self.reference_start,
            constants.cigar_consumption,
        )


    def get_blocks(self):
        """Mimic the get_blocks function from AlignedSegment.

        Notes
        =====
        - Takes roughly 200us
        """

        blocks = []
        block_start = self.reference_start
        block_length = 0

        for _, length, consumes_read, consumes_ref in iterate_cigartuples(self.cigartuples, constants.cigar_consumption):
            if consumes_read and consumes_ref:
                block_length += length

            elif consumes_read and not consumes_ref:
                if block_length:
                    blocks.append((block_start, block_start + block_length))

                block_start = block_start + block_length
                block_length = 0

            elif not consumes_read and consumes_ref:
                if block_length:
                    blocks.append((block_start, block_start + block_length))

                block_start = block_start + block_length + length
                block_length = 0

            else:
                pass

        if block_length:
            blocks.append((block_start, block_start + block_length))

        return blocks


    def __repr__(self):
        """Fancy output for viewing a read's alignment in relation to the reference"""

        read, ref = self.alignment_as_string_pair()

        lines = [
            '<%s.%s object at %s>' % (self.__class__.__module__, self.__class__.__name__, hex(id(self))),
            ' ├── start, end : [%s, %s)' % (self.reference_start, self.reference_end),
            ' ├── cigartuple : %s' % [tuple(row) for row in self.cigartuples],
            ' ├── read       : %s' % ''.join(read),
            ' └── reference  : %s' % ''.join(ref),
        ]

        return '\n'.join(lines)


    def alignment_as_string_pair(self):
        """Only for visualization purposes--very slow"""
        ref, read, pos_ref, pos_read = [], [], 0, 0
        for _, length, consumes_read, consumes_ref in iterate_cigartuples(self.cigartuples, constants.cigar_consumption):
            if consumes_read:
                read.extend([chr(x) for x in self.query_sequence[pos_read:(pos_read + length)]])
                pos_read += length
            else:
                read.extend(['-'] * length)

            if consumes_ref:
                ref.extend([chr(x) for x in self.reference_sequence[pos_ref:(pos_ref + length)]])
                pos_ref += length
            else:
                ref.extend(['-'] * length)

        return ''.join(read), ''.join(ref)


    def trim(self, trim_by, side='left'):
        """Trims self.read by either the left or right

        Modifies the attributes:

            query_sequence
            reference_sequence
            cigartuples
            reference_start
            reference_end

        Do not expect more than this!

        Parameters
        ==========
        trim_by : int
            The number of REFERENCE bases you would like to trim the read by. If the trim leaves
            operations that are consumed by the reference but not the read, or the read but not the
            reference, these are trimmed AS WELL. For example, if after trimming by `trim_by`, the
            final cigar string is [(2,2),(0,4)], this will be further trimmed to [(0,4)], since
            there is no useful information held in a terminal read gap.

        side : str, 'left'
            Either 'left' or 'right' side.
        """
        if trim_by == 0:
            return

        elif trim_by < 0:
            raise ConfigError("Read.trim :: Requesting to trim an amount %d, which is negative." % trim_by)

        elif trim_by > self.reference_end - self.reference_start:
            raise ConfigError("Read.trim :: Requesting to trim an amount %d that exceeds the alignment"
                              " range of %d" % (trim_by, self.reference_end - self.reference_start))

        if self.cigartuples.shape[0] == 1:
            # There contains only a pure mapping segment, i.e. no indels. This clause accounts for
            # the majority of reads and exists to speed up the code.
            self.cigartuples[0, 1] -= trim_by

            if side == 'left':
                self.query_sequence = self.query_sequence[trim_by:]
                self.reference_sequence = self.reference_sequence[trim_by:]
                self.reference_start += trim_by

            else:
                self.query_sequence = self.query_sequence[:-trim_by]
                self.reference_sequence = self.reference_sequence[:-trim_by]
                self.reference_end -= trim_by

            return

        # We are here because the read was not a simple mapping. There are indels and so we need to
        # parse cigartuples. We delegate to a just-in-time compiled function for a 4X speed gain

        (self.cigartuples,
         self.query_sequence,
         self.reference_sequence,
         self.reference_start,
         self.reference_end) = _trim(self.cigartuples,
                                     constants.cigar_consumption,
                                     self.query_sequence,
                                     self.reference_sequence,
                                     self.reference_start,
                                     self.reference_end,
                                     trim_by,
                                     0 if side == 'left' else 1)


class Coverage:
    def __init__(self):
        self.c = None # becomes a numpy array of coverage values
        self.min = 0
        self.max = 0
        self.std = 0.0
        self.mean = 0.0
        self.median = 0.0
        self.detection = 0.0
        self.mean_Q2Q3 = 0.0
        self.num_reads = 0

        self.read_iterator_dict = {
            'fetch': self._fetch_iterator,
            'fetch_and_trim': self._fetch_and_trim_iterator,
            'fetch_filter_and_trim': self._fetch_filter_and_trim_iterator,
        }


    def run(self, bam, contig_or_split, start=None, end=None, read_iterator='fetch',
            max_coverage=None, skip_coverage_stats=False, **kwargs):
        """Loop through the bam pileup and calculate coverage over a defined region of a contig or split

        Parameters
        ==========
        bam : bamops.BAMFileObject
            Init such an object the way you would a pysam.AlignmentFile, i.e. bam =
            bamops.BAMFileObject(path_to_bam)

        contig_or_split : anvio.contigops.Split or anvio.contigops.Contig or str
            If Split object is passed, and `start` or `end` are None, they are automatically set to
            contig_or_split.start and contig_or_split.end. If str object is passed, it is assumed to
            be a contig name

        start : int
            The index start of where coverage is calculated. Relative to the contig, even when
            `contig_or_split` is a Split object.

        end : int
            The index end of where coverage is calculated. Relative to the contig, even when
            `contig_or_split` is a Split object.

        read_iterator : string
            How do you want to iterate through reads? Options: see self.read_iterator_dict. Some
            notes: If you supply a start and end, you _must_ supply a read iterator that trims
            reads, i.e, not 'fetch', which returns the entirety of a read, even if it only partially
            lands in [start, end)

        skip_coverage_stats : bool, False
            Should the call to process_c be skipped?

        **kwargs : **kwargs
            kwargs are passed to the method chosen
        """

        if isinstance(contig_or_split, anvio.contigops.Split):
            contig_name = contig_or_split.parent
            start = contig_or_split.start if not start else start
            end = contig_or_split.end if not end else end

        elif isinstance(contig_or_split, anvio.contigops.Contig):
            contig_name = contig_or_split.name
            start = 0 if not start else start
            end = contig_or_split.length if not end else end

        elif isinstance(contig_or_split, str):
            contig_name = contig_or_split
            start = 0 if not start else start
            end = bam.get_reference_length(contig_name) if not end else end

        else:
            raise ConfigError("Coverage.run :: You can't pass an object of type %s as contig_or_split" % type(contig_or_split))

        # a coverage array the size of the defined range is allocated in memory
        c = np.zeros(end - start).astype(int)

        try:
            read_iterator = self.read_iterator_dict[read_iterator]
        except KeyError:
            raise ConfigError("Coverage :: %s is not a valid read iterator." % read_iterator)

        for read in read_iterator(bam, contig_name, start, end, **kwargs):
            if len(read.cigartuples) == 1:
                c[read.reference_start:(read.reference_start + read.cigartuples[0][1])] += 1
            else:
                for start, end in read.get_blocks():
                    c[start:end] += 1

            self.num_reads += 1

        self.c = c

        if max_coverage is not None:
            if np.max(self.c) > max_coverage:
                self.c[self.c > max_coverage] = max_coverage

        if len(self.c):
            try:
                contig_or_split.explicit_length = len(self.c)
            except AttributeError:
                pass

            if not skip_coverage_stats:
                self.process_c(self.c)


    def _fetch_iterator(self, bam, contig_name, start, end):
        """Uses standard pysam fetch iterator from AlignmentFile objects, ignores unmapped reads"""

        for read in bam.fetch_only(contig_name, start, end):
            if read.cigartuples is None or read.query_sequence is None:
                # This read either has no associated cigar string or no query sequence. If cigar
                # string is None, this means it did not align but is in the BAM file anyways, or the
                # mapping software decided it did not want to include a cigar string for this read.
                # The origin of why query sequence would be None is unkown.
                continue

            yield read


    def _fetch_and_trim_iterator(self, bam, contig_name, start, end):
        """Uses BAMFileObject.fetch_and_trim iterator"""

        for read in bam.fetch_and_trim(contig_name, start, end):
            yield read


    def _fetch_filter_and_trim_iterator(self, bam, contig_name, start, end, percent_id_cutoff):
        """Uses BAMFileObject.fetch_and_filter iterator"""

        for read in bam.fetch_filter_and_trim(contig_name, start, end, percent_id_cutoff):
            yield read


    def process_c(self, c):
        self.min = np.amin(c)
        self.max = np.amax(c)
        self.median = np.median(c)
        self.mean = np.mean(c)
        self.std = np.std(c)
        self.detection = np.sum(c > 0) / len(c)

        self.is_outlier = utils.get_list_of_outliers(c, median=self.median) # this is an array not a list

        if c.size < 4:
            self.mean_Q2Q3 = self.mean
        else:
            sorted_c = sorted(c)
            Q = int(c.size * 0.25)
            Q2Q3 = sorted_c[Q:-Q]
            self.mean_Q2Q3 = np.mean(Q2Q3)


class LinkMerDatum:
    def __init__(self, sample_id, read_id, is_read1):
        self.read_id = read_id
        self.sample_id = sample_id
        self.read_X = 'read-1' if is_read1 else 'read-2'
        self.read_unique_id = hashlib.sha224((sample_id + read_id + self.read_X).encode('utf-8')).hexdigest()
        self.contig_name = None
        self.request_id = None
        self.pos_in_contig = None
        self.pos_in_read = None
        self.base = None
        self.reverse = None
        self.sequence = None


    def __str__(self):
        return 'Contig: %s, C_pos: %d, R_pos: %d, Base: %s, seq: %s' % (self.contig_name,
                                                                        self.pos_in_contig,
                                                                        self.pos_in_read,
                                                                        self.base,
                                                                        self.sequence)


class LinkMersData:
    def __init__(self, run=run, progress=progress, quiet=False):
        self.data = []
        self.quiet = quiet

        self.run = run
        self.progress = progress


    def append(self, bam_file_object, sample_id, request_id, contig_name, positions, only_complete_links=False):
        data = []

        if not self.quiet:
            self.run.warning('', header="Working on '%s'" % sample_id, lc='cyan')

            self.progress.new('Processing "%s" in "%s"' % (contig_name, sample_id))
            self.progress.update('Analyzing %d positions (stretching %d nts) ...' % (len(positions), max(positions) - min(positions)))

        for pileupcolumn in bam_file_object.pileup(contig_name, min(positions), max(positions) + 1):
            if pileupcolumn.pos not in positions:
                continue

            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del:
                    L = LinkMerDatum(sample_id, pileupread.alignment.qname, pileupread.alignment.is_read1)
                    L.request_id = request_id
                    L.contig_name = contig_name
                    L.pos_in_contig = pileupcolumn.pos
                    L.pos_in_read = pileupread.query_position
                    L.base = pileupread.alignment.seq[pileupread.query_position]
                    L.reverse = pileupread.alignment.is_reverse
                    L.sequence = pileupread.alignment.query
                    # ['aend', 'alen', 'aligned_pairs', 'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples', 'compare', 'flag', 'get_aligned_pairs', 'get_blocks', 'get_overlap',
                    #  'get_reference_positions', 'get_reference_sequence', 'get_tag', 'get_tags', 'has_tag', 'infer_query_length', 'inferred_length', 'is_duplicate', 'is_paired',
                    #  'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary', 'is_unmapped', 'isize', 'mapping_quality', 'mapq',
                    #  'mate_is_reverse', 'mate_is_unmapped', 'mpos', 'mrnm', 'next_reference_id', 'next_reference_name', 'next_reference_start', 'opt', 'overlap', 'pnext', 'pos',
                    #  'positions', 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end', 'query_alignment_length', 'query_alignment_qualities',
                    #  'query_alignment_sequence', 'query_alignment_start', 'query_length', 'query_name', 'query_qualities', 'query_sequence', 'reference_end', 'reference_id',
                    #  'reference_length', 'reference_name', 'reference_start', 'rlen', 'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags',
                    #  'template_length', 'tid', 'tlen', 'tostring']
                    data.append(L)

        if not self.quiet:
            self.progress.end()
            self.run.info('data', '%d entries identified mapping at least one of the nucleotide positions for "%s"' % (len(data), contig_name))

        if only_complete_links:
            num_positions = len(positions)
            num_hits_dict = Counter([d.read_unique_id for d in data])
            read_unique_ids_to_keep = set([read_unique_id for read_unique_id in num_hits_dict if num_hits_dict[read_unique_id] == num_positions])

            if not self.quiet:
                self.run.info('data', '%s unique reads that covered all positions were kept' % (len(read_unique_ids_to_keep)), mc='red')

            self.data.append((contig_name, positions, [d for d in data if d.read_unique_id in read_unique_ids_to_keep]))
        else:
            self.data.append((contig_name, positions, data))


class LinkMers:
    """Get linkmers

    This class handles an input BAM file, a list of contigs, and positions within them to
    report bases in reads that contribute to positions of interest following Chris Quince's
    suggestion. Each read is reported with a unique ID, therefore linkage informaiton can
    be followed.

    FIXME This class might produce dubious results, according to the side-by-side comparions between
    how SCVs used to be calculated, and IGV view. For details please see
    https://github.com/merenlab/anvio/pull/1362. Working with reads directly rather than pileups is
    highly preferable to this framework and should be implemented the same way SCV calculations were
    changed
    """

    def __init__(self, args=None):
        self.args = args
        self.input_file_paths = []
        self.contig_and_position_requests_list = []

        self.progress = terminal.Progress()
        self.run = terminal.Run(width=35)

        if args:
            for input_file_path in args.input_files:
                filesnpaths.is_file_exists(input_file_path)

            self.input_file_paths = [os.path.abspath(p.strip()) for p in args.input_files]

            if len(self.input_file_paths) != len(set(self.input_file_paths)):
                raise ConfigError("You can't declared the same BAM file twice :/")

            self.only_complete_links = args.only_complete_links

            if args.list_contigs:
                self.list_contigs()
                sys.exit()

            filesnpaths.is_file_exists(args.contigs_and_positions)
            filesnpaths.is_file_tab_delimited(args.contigs_and_positions, expected_number_of_fields=2)

            request_id = 0
            f = open(args.contigs_and_positions)
            for line in f.readlines():
                request_id += 1

                contig_name, positions = line.split('\t')

                try:
                    positions = [int(pos) for pos in positions.split(',')]
                except ValueError:
                    raise ConfigError('Positions for contig "%s" does not seem to be comma-separated integers...' % contig_name)

                self.contig_and_position_requests_list.append((request_id, contig_name, set(positions)),)

        self.linkmers = None


    def process(self):
        self.sanity_check()

        self.linkmers = LinkMersData(self.run, self.progress)

        for input_file in self.input_file_paths:
            sample_id = filesnpaths.get_name_from_file_path(input_file)
            bam_file_object = BAMFileObject(input_file)

            self.run.info('input_bam_path', input_file)
            self.run.info('sample_id', sample_id)
            self.run.info('total_reads_mapped', pp(int(bam_file_object.mapped)))
            self.run.info('num_contigs_in_bam', pp(len(bam_file_object.references)))

            for request_id, contig_name, positions in self.contig_and_position_requests_list:
                self.linkmers.append(bam_file_object, sample_id, request_id, contig_name, positions, self.only_complete_links)

            bam_file_object.close()

        return self.linkmers.data


    def report(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        output_file = open(output_file_path, 'w')
        output_file.write('\t'.join(['entry_id', 'sample_id', 'request_id', 'contig_name', 'pos_in_contig',\
                                     'pos_in_read', 'base', 'read_unique_id', 'read_X', 'reverse',\
                                     'sequence']) + '\n')
        entry_id = 0
        for contig_name, positions, data in self.linkmers.data:
            for d in data:
                entry_id += 1
                output_file.write('%.9d\t%s\t%.3d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n'
                                % (entry_id, d.sample_id, d.request_id, d.contig_name,\
                                   d.pos_in_contig, d.pos_in_read, d.base, d.read_unique_id,\
                                   d.read_X, d.reverse, d.sequence))
        output_file.close()

        self.run.info('output_file', output_file_path)


    def list_contigs(self):
        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        bam_file_object = BAMFileObject(self.input_file_paths[0])
        self.progress.end()

        contig_names = bam_file_object.references
        contig_lengths = bam_file_object.lengths

        bam_file_object.close()

        for tpl in sorted(zip(contig_lengths, contig_names), reverse=True):
            print('%-40s %s' % (tpl[1], pp(int(tpl[0]))))


    def sanity_check(self):
        if not self.contig_and_position_requests_list:
            raise ConfigError("Entries dictionary is empty")


class GetReadsFromBAM:
    """This class takes contig names (or a collection ID and bins), and provides
       access short reads in one or more BAM files mapping to those contigs
       (expanded in #173)."""

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        self.args = args

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # groupA inputs (all the anvi'o files)
        self.input_bam_files = A('input_bams')
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')

        # grouB inputs (all the ad hoc requests)
        self.target_contig = A('target_contig')
        self.target_region_start = A('target_region_start')
        self.target_region_end = A('target_region_end')

        # we can has fetch filters?
        self.fetch_filter = A('fetch_filter')

        self.output_file_path = A('output_file')
        self.output_file_prefix = A('output_file_prefix')
        self.gzip = A('gzip_output')
        self.split_R1_and_R2 = A('split_R1_and_R2')

        self.bins = set([])
        self.split_names_of_interest = set([])


    def sanity_check(self):
        """A function to ensure everything is order given the `args`"""

        ##############################################################
        # make sure BAM files are legit
        ##############################################################
        bad_bam_files = []
        error_message = None
        for bam_file_path in self.input_bam_files:
            try:
                bam_file_object = BAMFileObject(bam_file_path)
                bam_file_object.close()
            except ConfigError as e:
                bad_bam_files.append(bam_file_path)
                error_message = e

        if len(bad_bam_files):
            raise ConfigError('Samtools is not happy with some of your bam files. The following '
                              'file(s) do not look like proper BAM files [here is the actual '
                              'error: "%s"]: %s.' % (error_message, ','.join(bad_bam_files)))

        ##############################################################
        # make sure the input files look OK.
        ##############################################################
        if self.profile_db_path or self.contigs_db_path:
            ##############################################################
            # the user has chosen the anvi'o files path
            ##############################################################
            if self.target_contig:
                raise ConfigError("If you choose to go with anvi'o databases to determine which short reads "
                                  "to get from your BAM files, then you can't also provide a target contig "
                                  "name.")

            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)
        elif self.target_contig:
            ##############################################################
            # the user has chosen the ad hoc inputs path
            ##############################################################
            if self.profile_db_path or self.collection_name or self.bin_id or self.bin_ids_file_path:
                raise ConfigError("If you choose to go with an ad hoc contig name to get your short reads for, "
                                  "you can't also provide anvi'o databases, collections, or bin names. Stop "
                                  "confusing anvi'o :(")

        ##############################################################
        # make sure output files are OK
        ##############################################################
        if self.output_file_prefix and self.output_file_path:
            raise ConfigError("You must either use the parameter output file name, or output file prefix.")

        if self.output_file_prefix and not self.split_R1_and_R2:
            raise ConfigError("Output file prefix parameter is only relevant when you want to split R1 reads "
                              "from R2 reads and so on.")

        if self.split_R1_and_R2 and not self.output_file_prefix:
            raise ConfigError("If you wish R1 and R2 reads to be reported in separate FASTA files, \
                               you need to provide an output file prefix so anvi'o can generate\
                               multiple output files that start with it (i.e., PREFIX_R1.fa, PREFIX_R2.fa\
                               PREFIX_UNPAIRED.fa).")

        if self.split_R1_and_R2:
            filesnpaths.is_output_file_writable(self.output_file_prefix + '_R1.fa')
        elif self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)
        else:
            filesnpaths.is_output_file_writable('short_reads.fa')


    def get_short_reads_dict(self, contig_start_stops):
        """Gets short reads from BAM files, literally.

        If `contig_start_stops` is None, the function will learn the contig names in the first
        BAM file and generate a `contig_start_stops` list.
        """

        self.run.info('Input BAM file(s)', ', '.join([os.path.basename(f) for f in self.input_bam_files]), nl_before=1)
        self.run.info('Fetch filter', self.fetch_filter, nl_after=1, mc=('red' if not self.fetch_filter else 'green'))

        # if the user has provided no `contig_start_stops`, then we will have to recover all contig names
        # from the BAM file and use them as our target.
        if not contig_start_stops:
            bam_file_path_to_recover_contig_names_from = self.input_bam_files[0]
            self.progress.new('Recovering contig names')
            self.progress.update(f"from {bam_file_path_to_recover_contig_names_from} ...")
            bam_file_object = BAMFileObject(bam_file_path_to_recover_contig_names_from)
            contig_start_stops = [(contig_name, 0, sys.maxsize) for contig_name in bam_file_object.references]
            bam_file_object.close()
            self.progress.end()

        # prepare the dictionary to be returned
        short_reads_dict = {}
        if self.split_R1_and_R2:
            short_reads_dict['R1'] = {}
            short_reads_dict['R2'] = {}
            short_reads_dict['UNPAIRED'] = {}
        else:
            short_reads_dict['all'] = {}

        self.progress.new("Recovering short reads")

        # at this point contig_start_stops knows every contig we are interested in, and
        # their start and stop positions based on what split ids were requested. we
        # shall go through each bam file the user is interested, and get those short reads
        # that map to regions of interest:
        for bam_file_path in self.input_bam_files:
            bam_file_name = filesnpaths.get_name_from_file_path(bam_file_path)

            bam_file_object = BAMFileObject(bam_file_path)
            bam_file_object.fetch_filter = self.fetch_filter

            for contig_id, start, stop in contig_start_stops:
                if contig_id not in bam_file_object.references:
                    raise ConfigError(f"The contig name '{contig_id}' is not in the BAM file at '{bam_file_path}' :(")

            self.progress.update('Creating a dictionary of matching short reads in %s ...' % bam_file_name)

            '''here's what's available in the read objects below:

            ['aend', 'alen', 'aligned_pairs', 'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples', 'compare',
             'flag', 'get_aligned_pairs', 'get_blocks', 'get_overlap', 'get_reference_positions', 'get_tag',
             'get_tags', 'has_tag', 'infer_query_length', 'inferred_length', 'is_duplicate', 'is_paired',
             'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary',
             'is_unmapped', 'isize', 'mapping_quality', 'mapq', 'mate_is_reverse', 'mate_is_unmapped', 'mpos', 'mrnm',
             'next_reference_id', 'next_reference_start', 'opt', 'overlap', 'pnext', 'pos', 'positions', 'qend',
             'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end', 'query_alignment_length',
             'query_alignment_qualities', 'query_alignment_sequence', 'query_alignment_start', 'query_length',
             'query_name', 'query_qualities', 'query_sequence', 'reference_end', 'reference_id', 'reference_length',
             'reference_start', 'rlen', 'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags',
             'template_length', 'tid', 'tlen']'''

            D = lambda sequence: {'seq': sequence, 'contig': contig_id, 'start': str(read.get_reference_positions()[0]), 'stop': str(read.get_reference_positions()[-1]), 'bam': bam_file_name, 'read': read.query_name}

            # this will serve as a unique id per read
            counter = 0
            unknown_mate_tracker_dict = {}
            unknown_mate_defline_to_counter = {}
            if self.split_R1_and_R2:
                for contig_id, start, stop in contig_start_stops:

                    # learn contig length & ensure the stop position is not longer than the contig length
                    contig_length = bam_file_object.lengths[bam_file_object.references.index(contig_id)]
                    stop = contig_length if stop > contig_length else stop

                    for read in bam_file_object.fetch_only(contig_id, start, stop):
                        counter += 1

                        defline = '_'.join([contig_id, str(start), str(stop), read.query_name, bam_file_name])

                        if not read.is_paired:
                            short_reads_dict['UNPAIRED'][counter] = D(read.query_sequence)

                        elif defline in unknown_mate_tracker_dict:
                            # `read`s mate has already been read. so assign the read and the mate
                            # to their respective 'R1' and 'R2' dictionaries, then remove the mate
                            # from unknown_mate_tracker_dict since its mate is now known.
                            read_DIRECTION = 'R1' if read.is_read1 else 'R2'
                            mate_DIRECTION = 'R2' if read_DIRECTION == 'R1' else 'R1'

                            # rev_comp either R1 or R2 to match the original short reads orientation
                            if read.is_reverse:
                                short_reads_dict[mate_DIRECTION][counter] = D(unknown_mate_tracker_dict[defline]['seq'])
                                short_reads_dict[read_DIRECTION][counter] = D(utils.rev_comp(read.query_sequence))
                            else:
                                short_reads_dict[mate_DIRECTION][counter] = D(utils.rev_comp(unknown_mate_tracker_dict[defline]['seq']))
                                short_reads_dict[read_DIRECTION][counter] = D(read.query_sequence)

                            del unknown_mate_tracker_dict[defline]
                            del unknown_mate_defline_to_counter[defline]

                        else:
                            unknown_mate_tracker_dict[defline] = D(read.query_sequence)
                            unknown_mate_defline_to_counter[defline] = counter


                short_reads_dict['UNPAIRED'].update(unknown_mate_tracker_dict)
            else:
                for contig_id, start, stop in contig_start_stops:
                    # sanity check for the stop position
                    contig_length = bam_file_object.lengths[bam_file_object.references.index(contig_id)]
                    stop = contig_length if stop > contig_length else stop

                    for read in bam_file_object.fetch_only(contig_id, start, stop):
                        counter += 1
                        short_reads_dict['all'][counter] = D(read.query_sequence)

            bam_file_object.close()

        self.progress.end()

        return short_reads_dict


    def get_all_short_reads_from_bam(self):
        """A route to get all short reads from BAM files without anvi'o files or specific contigs"""


        self.run.info('Mode', "Nothing is provided -- will report all reads", mc="green")

        return self.get_short_reads_dict(contig_start_stops=None)


    def get_short_reads_for_contig_and_range(self):
        """A route get short reads based on the contig name and/or start/stop positions defined by the user"""

        if not self.target_contig:
            raise ConfigError("You can't be here without a target contig name :/")

        self.run.info('Mode', "Contig name provided by the user", mc="green")
        self.run.info('Contig name', self.target_contig)
        self.run.info('Target region start', self.target_region_start if self.target_region_start else "Start from the beginning")
        self.run.info('Target region end', self.target_region_end if self.target_region_end else "Go all the way to the end of the contig")

        contig_start_stops = [(self.target_contig, self.target_region_start or 0, self.target_region_end or sys.maxsize)]

        return self.get_short_reads_dict(contig_start_stops)


    def get_short_reads_through_anvio_files(self, split_names_of_interest=None):
        """A route get short reads based on anvi'o collection/bin names or for a given splits dict"""

        if split_names_of_interest:
            # programmer passed the split names of interest manually
            self.split_names_of_interest = split_names_of_interest
            self.run.info('Mode', "Split names are passed to the function as a variable", mc="green")
        elif len(self.split_names_of_interest) > 0:
            # programmer passed the split names of interest by setting a member variable
            self.run.info('Mode', "Split names are passed as a class member variable", mc="green")
        else:
            # we will recover split names of interest from anvi'o files
            self.run.info('Mode', "Recovering split names from anvi'o files", mc="green")
            self.run.info('Profile DB', self.profile_db_path)
            self.run.info('Contigs DB', self.contigs_db_path)
            self.run.info('Collection ID', self.collection_name)

            d = ccollections.GetSplitNamesInBins(self.args).get_dict()
            self.bins = list(d.keys())

            for split_names in list(d.values()):
                self.split_names_of_interest.update(split_names)

            self.run.info('Bin(s)', ', '.join(self.bins))

        if not len(self.split_names_of_interest):
            raise ConfigError("The split names of interest set is empty. This should have never happened. Good "
                              "job.")
        else:
            self.run.info('Number of splits', pp(len(self.split_names_of_interest)), nl_after=1)

        self.progress.new('Identifying contig start/stops')
        self.progress.update('Reading splits info from the contigs database ...')
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)
        contigs_db.disconnect()

        self.progress.update('Identifying contigs associated with splits ...')
        contigs_involved = utils.get_contigs_splits_dict(self.split_names_of_interest, splits_basic_info)

        # this variable will hold a list of (contig_id, start, stop) tuples
        # for each contig and the start and stop positions of sequential blocks
        # of splits identified within them
        contig_start_stops = []

        self.progress.update('Computing start/stops positions of interest in %d contigs ...' % (len(contigs_involved)))
        for contig_id in contigs_involved:
            splits_order = list(contigs_involved[contig_id].keys())
            sequential_blocks = ccollections.GetSequentialBlocksOfSplits(splits_order).process()

            for sequential_block in sequential_blocks:
                first_split = contigs_involved[contig_id][sequential_block[0]]
                last_split = contigs_involved[contig_id][sequential_block[-1]]

                contig_start_stops.append((contig_id,
                                           splits_basic_info[first_split]['start'],
                                           splits_basic_info[last_split]['end']),)

        self.progress.end()

        return self.get_short_reads_dict(contig_start_stops)


    def get_short_reads(self):
        """Function that determines different routes to get short reads.

        A broker of sorts. If you have a nicely prepared `args` object, this is the function
        you should be using.
        """

        self.sanity_check()

        if self.profile_db_path or self.contigs_db_path:
            short_reads_dict = self.get_short_reads_through_anvio_files()
        elif self.target_contig:
            short_reads_dict = self.get_short_reads_for_contig_and_range()
        else:
            short_reads_dict = self.get_all_short_reads_from_bam()

        return short_reads_dict


    def store_short_reads(self):
        short_reads_dict = self.get_short_reads()

        self.progress.new("Storing reads")
        self.progress.update("...")

        if self.split_R1_and_R2:
            for read_type in sorted(list(short_reads_dict.keys())):
                output_file_path = f"{self.output_file_prefix}_{read_type}.fa"
                with open(output_file_path, 'w') as output:
                    for entry_id in short_reads_dict[read_type]:
                        e = short_reads_dict[read_type][entry_id]
                        output.write(f">{entry_id} read_id:{e['read']}|bam:{e['bam']}|contig:{e['contig']}|start:{e['start']}|stop:{e['stop']}\n{e['seq']}\n")

                if self.gzip:
                    utils.gzip_compress_file(output_file_path)
                    output_file_path = output_file_path + ".gz"

                self.run.info('Output file for %s' % read_type, output_file_path, progress=self.progress)

            self.progress.end()
            self.run.info('Num paired-end reads stored', pp(len(short_reads_dict['R1'])), mc='green', nl_before=1)
            self.run.info('Num unpaired reads stored', pp(len(short_reads_dict['UNPAIRED'])), mc='green')
        else:
            output_file_path = self.output_file_path or 'short_reads.fa'

            with open(output_file_path, 'w') as output:
                for entry_id in short_reads_dict['all']:
                    e = short_reads_dict['all'][entry_id]
                    output.write(f">{entry_id} read_id:{e['read']}|bam:{e['bam']}|contig:{e['contig']}|start:{e['start']}|stop:{e['stop']}\n{e['seq']}\n")

            if self.gzip:
                utils.gzip_compress_file(output_file_path)
                output_file_path = output_file_path + ".gz"

            self.progress.end()
            self.run.info('Output file for all short reads',output_file_path)
            self.run.info('Num reads stored', pp(len(short_reads_dict['all'])), mc='green')


class ReadsMappingToARange:
    """Returns all reads from BAM that maps to a range in a contig"""

    def __init__(self, run=run, progress=progress):
        self.data = []

        self.run = run
        self.progress = progress


    def process_range(self, input_bam_paths, contig_name, start, end):
        if end <= start:
            progress.reset()
            raise ConfigError("The end of range cannot be equal or smaller than the start of it :/")

        data = []

        for input_bam_path in input_bam_paths:
            bam_file_object = BAMFileObject(input_bam_path)

            if contig_name not in bam_file_object.references:
                progress.reset()
                raise ConfigError(f"Contig name '{contig_name}' is not represented in the BAM file at '{input_bam_path}'. "
                                  f"Anvi'o is unable to make a intelligent guess regarding how did you end up here, but "
                                  f"the only way for this to happen is that either (1) you have a contigs database that "
                                  f"contains more contig names than the FASTA file you have used for the read recruitment "
                                  f"that gave you this BAM file, or (2) the settings you have used with the read recruitment "
                                  f"tool in fact did not report results for some contigs :/ Either way. There is no where "
                                  f"to go from here :(")

            sample_id = filesnpaths.get_name_from_file_path(input_bam_path)

            self.run.warning('', header="Working on '%s'" % sample_id, lc='cyan')

            self.run.info('input_bam_path', input_bam_path)
            self.run.info('sample_id', sample_id)
            self.run.info('total_reads_mapped', pp(int(bam_file_object.mapped)))
            self.run.info('num_contigs_in_bam', pp(len(bam_file_object.references)))

            self.progress.new('Processing "%s" in "%s"' % (contig_name, input_bam_path))
            self.progress.update('Analyzing positions stretching %d nts ...' % (end - start))

            read_ids = set([])

            for pileupcolumn in bam_file_object.pileup(contig_name, start, end):
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del:
                        L = LinkMerDatum(sample_id, pileupread.alignment.qname, pileupread.alignment.is_read1)
                        L.contig_name = contig_name
                        L.pos_in_contig = pileupcolumn.pos
                        L.pos_in_read = pileupread.query_position
                        L.base = pileupread.alignment.seq[pileupread.query_position]
                        L.reverse = pileupread.alignment.is_reverse
                        L.sequence = pileupread.alignment.query

                        if not L.read_unique_id in read_ids:
                            data.append(L)
                            read_ids.add(L.read_unique_id)

            self.progress.end()
            bam_file_object.close()

        self.run.info('data', '%d reads identified mapping between positions %d and %d in "%s"' % (len(data), start, end, contig_name))

        self.data.extend(data)


    def report(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        output_file = open(output_file_path, 'w')
        entry_id = 0
        for d in self.data:
            entry_id += 1
            output_file.write('>%.9d|sample_id:%s|reverse:%s|contig_name:%s\n' % (entry_id, d.sample_id, d.reverse, d.contig_name))
            output_file.write('%s\n' % (d.sequence))
        output_file.close()

        self.run.info('output_file', output_file_path)


@dataclass
class CircularityResult:
    """Results from circularity prediction for a single contig.

    Attributes
    ----------
    contig_name : str
        Name of the analyzed contig.

    status : str
        Classification result: 'circular', 'linear', or 'indeterminate'.

    circularity_support : float
        Score from 0.0 to 1.0 representing the strength of evidence for circularity.
        Computed as supporting_rf_pairs / observed_rf_pairs.

    edge_coherence : float
        Fraction of edge-mapping reads whose mates also map to the same contig.
        High values suggest that the vast majority of paired-end reads mapping
        to the contig map together on it, low values suggest there are many
        pairs where only one of the reads map to this contig (which is a good
        indication of fragmentation).

    observed_rf_pairs : int
        Total number of RF-oriented read pairs found on this contig.

    supporting_rf_pairs : int
        Number of RF pairs whose circular insert size falls within the expected
        insert size distribution (i.e., plausible junction-spanning pairs).

    expected_rf_pairs : float
        Expected number of RF pairs if the contig were circular, derived from
        FR pair count and the geometry of insert size vs contig length.

    num_fr_pairs : int
        Number of FR (forward-reverse) pairs used to estimate coverage.

    num_edge_reads : int
        Number of reads mapping within M bp of contig edges.

    num_edge_reads_incoherent : int
        Number of edge reads whose mate maps to a different contig.

    contig_length : int
        Length of the contig in base pairs.

    flags : list of str
        Warning flags for edge cases (e.g., low coverage, short contig).
    """

    contig_name: str
    status: str
    circularity_support: float
    edge_coherence: float
    observed_rf_pairs: int
    supporting_rf_pairs: int
    expected_rf_pairs: float
    num_fr_pairs: int
    num_edge_reads: int
    num_edge_reads_incoherent: int
    contig_length: int
    flags: list = field(default_factory=list)

    def __repr__(self):
        flag_str = f", flags={self.flags}" if self.flags else ""
        return (f"CircularityResult(contig='{self.contig_name}', status='{self.status}', "
                f"circularity_support={self.circularity_support:.3f}, edge_coherence={self.edge_coherence:.3f}, "
                f"supporting_rf={self.supporting_rf_pairs}, "
                f"expected_rf={self.expected_rf_pairs:.1f}{flag_str})")


class CircularityPredictor:
    """Predicts contig circularity from paired-end read mapping patterns.

    This class analyzes BAM alignments to determine whether contigs represent circular
    DNA molecules. It first computes global insert size statistics from FR (forward-reverse)
    read pairs across the entire BAM file, then uses these to evaluate whether RF
    (reverse-forward) pairs on individual contigs are consistent with a circular junction.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments namespace containing at minimum 'bam_file'. Optional arguments include:
          - max_num_pairs_for_is_est: Maximum FR pairs to sample for insert size estimation (default: 100000)
          - min_pairs_for_stats: Minimum FR pairs required for statistics (default: 1000)
          - insert_tolerance_factor: MADs from median for valid circular insert (default: 3.0)
          - min_supporting_pairs: Absolute minimum RF pairs to call circular (default: 5)
          - expected_fraction_threshold: Minimum fraction of expected RF pairs (default: 0.3)
          - circularity_support_threshold: Minimum support to call circular (default: 0.5)

    run : terminal.Run
        Anvi'o Run object for logging.

    progress : terminal.Progress
        Anvi'o Progress object for progress bars.

    Attributes
    ----------
    median_insert_size : float or None
        Median insert size computed from FR pairs. None until `process()` is called.

    mad_insert_size : float or None
        Median absolute deviation of insert sizes. None until `process()` is called.

    Examples
    --------
    >>> import argparse
    >>> args = argparse.Namespace(bam_file='metagenome.bam')
    >>> predictor = CircularityPredictor(args)
    >>> predictor.process()
    >>> result = predictor.predict('contig_00042')
    >>> print(result.status, result.circularity_support)
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bam_file_path = A('bam_file')

        # Parameters for insert size estimation
        self.max_num_pairs_for_is_est = A('max_num_pairs_for_is_est') or 100000
        self.min_pairs_for_stats = A('min_pairs_for_stats') or 1000

        # Parameters for circularity prediction
        self.insert_tolerance_factor = A('insert_tolerance_factor') or 3.0
        self.min_supporting_pairs = A('min_supporting_pairs') or 5
        self.expected_fraction_threshold = A('expected_fraction_threshold') or 0.3
        self.circularity_support_threshold = A('circularity_support_threshold') or 0.5

        # Will be populated by process()
        self.median_insert_size = None
        self.mad_insert_size = None
        self.bam = None
        self.contig_lengths = {}

        # Sanity checks
        if not self.bam_file_path:
            raise ConfigError("CircularityPredictor :: You must provide a BAM file path via `args.bam_file`.")

        filesnpaths.is_file_bam_file(self.bam_file_path)


    def process(self):
        """Compute global insert size statistics from the BAM file.

        This method iterates through properly-paired FR reads across all contigs and
        collects insert sizes to compute the median and MAD (median absolute deviation).
        These statistics are stored as instance attributes and used by `predict()`.

        This method must be called before `predict()`.

        Raises
        ------
        ConfigError
            If fewer than `min_pairs_for_stats` FR pairs are found in the BAM file.
        """

        self.run.info('BAM file', self.bam_file_path)

        self.bam = BAMFileObject(self.bam_file_path)

        # Store contig lengths for later use
        self.contig_lengths = dict(zip(self.bam.references, self.bam.lengths))

        self.run.info('Number of contigs', len(self.contig_lengths))
        self.run.info('Sample size for insert estimation', self.max_num_pairs_for_is_est)

        self._compute_insert_size_stats()

        self.run.info('Median insert size', f'{self.median_insert_size:.1f}')
        self.run.info('MAD insert size', f'{self.mad_insert_size:.1f}', nl_after=1)


    def _compute_insert_size_stats(self):
        """Compute median and MAD of insert sizes from FR pairs across the BAM.

        Iterates through properly-paired FR reads and collects insert sizes until
        `max_num_pairs_for_is_est` pairs are collected or all reads are exhausted.
        """

        num_pairs_processed = 0
        insert_sizes = []

        self.progress.new('Computing insert size statistics')
        self.progress.update('Sampling FR pairs across all contigs...')

        for contig_name in self.bam.references:

            if num_pairs_processed >= self.max_num_pairs_for_is_est:
                break

            for read in self.bam.fetch_only(contig_name):
                if not self._is_valid_fr_pair(read):
                    continue

                insert_size = abs(read.template_length)
                if insert_size > 0:
                    num_pairs_processed += 1
                    insert_sizes.append(insert_size)

                    if num_pairs_processed >= self.max_num_pairs_for_is_est:
                        break

        self.progress.end()

        if len(insert_sizes) < self.min_pairs_for_stats:
            raise ConfigError(
                f"CircularityPredictor :: Found only {len(insert_sizes):,} FR pairs, but need at least "
                f"{self.min_pairs_for_stats:,} for reliable insert size estimation. This BAM file may not "
                f"contain paired-end reads, or the reads may not be properly paired."
            )

        insert_sizes = np.array(insert_sizes)
        self.median_insert_size = float(np.median(insert_sizes))
        self.mad_insert_size = self._compute_mad(insert_sizes)

        self.run.info('FR pairs sampled', f'{len(insert_sizes):,}')


    def _compute_mad(self, values):
        """Compute the Median Absolute Deviation (MAD) of an array.

        MAD is a robust measure of dispersion that is less sensitive to outliers
        than standard deviation.

        Parameters
        ----------
        values : np.ndarray
            Array of numeric values.

        Returns
        -------
        float
            The MAD value. Returns 1.0 if MAD would be zero (to avoid division issues).
        """

        median = np.median(values)
        mad = float(np.median(np.abs(values - median)))

        # Guard against MAD of zero (highly uniform library)
        if mad == 0:
            return 1.0

        return mad


    def _sanity_check_insert_size(self):
        """Make sure insert size stats have been computed in places they are needed."""
        if self.median_insert_size is None or self.mad_insert_size is None:
            raise ConfigError("CircularityPredictor speaking. Insert size statistics are not computed :/ You must call "
                              "`process()` before `diagnose()`.")


    def predict(self, contig_name):
        """
        Predict whether a contig is circular based on RF pair evidence.

        For each RF (reverse-forward) pair on the contig, computes the "circular
        insert size"—the distance the reads would span if the contig were circular
        and linearized at an arbitrary point. RF pairs with circular insert sizes
        matching the expected distribution (within tolerance) are counted as
        supporting evidence for circularity.

        Additionally computes an edge coherence score for linear contigs based on
        whether reads near contig edges have mates mapping to the same contig
        (complete linear) or elsewhere (fragmented).

        Parameters
        ----------
        contig_name : str
            Name of the contig to analyze (must exist in BAM header).

        Returns
        -------
        CircularityResult
            Dataclass containing status, circularity_support, edge coherence, counts, and any warning flags.

        Raises
        ------
        ConfigError
            If `process()` has not been called, or if the contig is not found.
        """
        self._sanity_check_insert_size()

        if contig_name not in self.contig_lengths:
            raise ConfigError(f"CircularityPredictor :: Contig '{contig_name}' not found in BAM header. "
                              f"Available contigs: {len(self.contig_lengths):,}")

        flags = []
        M = self.median_insert_size
        MAD = self.mad_insert_size
        L = self.contig_lengths[contig_name]

        # Edge case: contig too short for meaningful analysis
        if L < 2 * M:
            flags.append(f"contig_shorter_than_2x_insert_size (L={L:,}, 2M={2*M:,.0f})")
            return CircularityResult(contig_name=contig_name,
                                     status='indeterminate',
                                     circularity_support=0.0,
                                     edge_coherence=0.0,
                                     observed_rf_pairs=0,
                                     supporting_rf_pairs=0,
                                     expected_rf_pairs=0.0,
                                     num_fr_pairs=0,
                                     num_edge_reads=0,
                                     num_edge_reads_incoherent=0,
                                     contig_length=L,
                                     flags=flags)

        # Define edge zones: within M bp of contig start or end
        edge_zone_start = M
        edge_zone_end = L - M

        # Collect read pairs mapping to this contig
        fr_count = 0
        rf_pairs = []  # List of (left_end, right_start) tuples

        # Track edge reads for edge coherence calculation
        num_edge_reads = 0
        num_edge_reads_incoherent = 0

        for read in self.bam.fetch_only(contig_name):
            # Basic validity: mapped, primary alignment (but don't require mate to be mapped yet)
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Check if this read is in an edge zone (for edge coherence calculation)
            read_in_edge_zone = (read.reference_start < edge_zone_start or 
                                read.reference_start >= edge_zone_end)

            if read_in_edge_zone:
                num_edge_reads += 1
                # Mate is incoherent if unmapped OR maps to different contig
                if read.mate_is_unmapped or read.reference_name != read.next_reference_name:
                    num_edge_reads_incoherent += 1

            # For FR/RF classification, we need both reads mapped to same contig
            if read.mate_is_unmapped:
                continue

            # Only process each pair once (when we see the leftmost read)
            if read.reference_start > read.next_reference_start:
                continue

            # Both reads must be on the same contig for FR/RF classification
            if read.reference_name != read.next_reference_name:
                continue

            orientation = self._get_pair_orientation(read)

            if orientation == 'FR':
                fr_count += 1

            elif orientation == 'RF':
                left_end = read.reference_end
                right_start = read.next_reference_start
                rf_pairs.append((left_end, right_start))

        # Compute edge coherence score
        if num_edge_reads > 0:
            edge_coherence = 1.0 - (num_edge_reads_incoherent / num_edge_reads)
        else:
            edge_coherence = 0.0

        # Calculate expected RF pairs if circular
        # E[RF] ≈ N_FR × M / (L - M)
        if fr_count == 0:
            flags.append("no_fr_pairs_on_contig")
            return CircularityResult(contig_name=contig_name,
                                     status='indeterminate',
                                     circularity_support=0.0,
                                     edge_coherence=edge_coherence,
                                     observed_rf_pairs=len(rf_pairs),
                                     supporting_rf_pairs=0,
                                     expected_rf_pairs=0.0,
                                     num_fr_pairs=0,
                                     num_edge_reads=num_edge_reads,
                                     num_edge_reads_incoherent=num_edge_reads_incoherent,
                                     contig_length=L,
                                     flags=flags)

        expected_rf = fr_count * M / (L - M)

        # Count RF pairs with plausible circular insert sizes
        tolerance = self.insert_tolerance_factor * MAD
        supporting_rf = 0

        for left_end, right_start in rf_pairs:
            # Circular insert: distance walking "around" the circle
            circular_insert = left_end + (L - right_start)

            if abs(circular_insert - M) <= tolerance:
                supporting_rf += 1

        # Compute circularity support score
        if len(rf_pairs) > 0:
            circularity_support = supporting_rf / len(rf_pairs)
        else:
            circularity_support = 0.0

        # Compute supporting fraction: of observed RF pairs, how many support circularity?
        if len(rf_pairs) > 0:
            supporting_fraction = supporting_rf / len(rf_pairs)
        else:
            supporting_fraction = 0.0

        # Dynamic minimum threshold: need both absolute and relative evidence
        min_required = max(self.min_supporting_pairs, self.expected_fraction_threshold * expected_rf)

        # Add diagnostic flags for edge cases
        if fr_count < 100:
            flags.append(f"low_fr_coverage (n={fr_count})")

        if expected_rf < self.min_supporting_pairs:
            flags.append(f"low_expected_rf_pairs (expected={expected_rf:.1f})")

        if num_edge_reads < 20:
            flags.append(f"low_edge_coverage (n={num_edge_reads})")

        # Make the classification decision
        if supporting_rf >= min_required and circularity_support >= self.circularity_support_threshold:
            status = 'circular'
        elif supporting_rf >= self.min_supporting_pairs and supporting_fraction >= 0.8:
            # Strong evidence: most observed RF pairs have correct circular insert size
            status = 'circular'
        elif len(rf_pairs) == 0 and fr_count >= 100:
            # No RF pairs at all with decent coverage suggests linear
            status = 'linear'
        elif supporting_rf < min_required and circularity_support < self.circularity_support_threshold:
            status = 'linear'
        else:
            status = 'indeterminate'
            flags.append("ambiguous_evidence")

        return CircularityResult(contig_name=contig_name,
                                 status=status,
                                 circularity_support=circularity_support,
                                 edge_coherence=edge_coherence,
                                 observed_rf_pairs=len(rf_pairs),
                                 supporting_rf_pairs=supporting_rf,
                                 expected_rf_pairs=expected_rf,
                                 num_fr_pairs=fr_count,
                                 num_edge_reads=num_edge_reads,
                                 num_edge_reads_incoherent=num_edge_reads_incoherent,
                                 contig_length=L,
                                 flags=flags)


    def diagnose(self, contig_name):
        """Return diagnostic information about RF pairs for a contig (when you silly shit).

        Useful for understanding why a contig with many RF pairs might have low circularity supprt. Here is an example:

        >>> from anvio.bamops import CircularityPredictor
        >>> import argparse
        >>> args = argparse.Namespace(bam_file='HADS_20210819_H1030_MGX_STO1_035.bam')
        >>> predictor = CircularityPredictor(args)
        >>> predictor.process()
        >>> for contig in ['STO1_000000177209', 'STO1_000000417019', 'STO1_000000433561', 'STO1_000000494581']:
        >>>     print(f"Contig: {contig}")
        >>>     print('-' * 80)
        >>>     diag = predictor.diagnose(i)
        >>>     for k, v in diag.items():
        >>>         print(f"{k}: {v}")
        >>>     print("\n\n")

        """

        self._sanity_check_insert_size()

        M = self.median_insert_size
        MAD = self.mad_insert_size
        L = self.contig_lengths[contig_name]
        tolerance = self.insert_tolerance_factor * MAD

        rf_pairs = []
        for read in self.bam.fetch_only(contig_name):
            if not self._is_valid_primary_read(read):
                continue
            if read.reference_start > read.next_reference_start:
                continue
            if read.reference_name != read.next_reference_name:
                continue

            if self._get_pair_orientation(read) == 'RF':
                left_end = read.reference_end
                right_start = read.next_reference_start
                rf_pairs.append((left_end, right_start))

        circular_inserts = [left_end + (L - right_start) for left_end, right_start in rf_pairs]

        if not circular_inserts:
            return {'contig_name': contig_name, 'num_rf_pairs': 0}

        ci = np.array(circular_inserts)

        d = {'contig_name': contig_name,
             'contig_length': L,
             'median_insert_size': M,
             'mad_insert_size': MAD,
             'tolerance': tolerance,
             'expected_range': (M - tolerance, M + tolerance),
             'num_rf_pairs': len(rf_pairs),
             'circular_insert_min': int(ci.min()),
             'circular_insert_max': int(ci.max()),
             'circular_insert_median': float(np.median(ci)),
             'circular_insert_mean': float(np.mean(ci)),
             'num_below_range': int(np.sum(ci < M - tolerance)),
             'num_above_range': int(np.sum(ci > M + tolerance))}

        return d


    def predict_all(self, min_contig_length=None, output_file=None):
        """Predict circularity for ALL contigs in the BAM file.

        Parameters
        ----------
        min_contig_length : int, optional
            If provided, skip contigs shorter than this length.

        output_file : str, optional
            If provided, write results to this path as a TAB-delimited file.

        Returns
        -------
        list of CircularityResult
            Results for all analyzed contigs.
        """

        self._sanity_check_insert_size()

        if output_file:
            filesnpaths.is_output_file_writable(output_file)

        results = []
        contigs_to_analyze = [name for name, length in self.contig_lengths.items() if min_contig_length is None or length >= min_contig_length]

        self.progress.new('Predicting circularity', progress_total_items=len(contigs_to_analyze))

        for contig_name in contigs_to_analyze:
            self.progress.increment()

            if self.progress.progress_current_item % 250 == 0:
                self.progress.update(f'... {contig_name} ...')

            result = self.predict(contig_name)
            results.append(result)

        self.progress.end()

        # Summary statistics
        circular_count = sum(1 for r in results if r.status == 'circular')
        linear_count = sum(1 for r in results if r.status == 'linear')
        indeterminate_count = sum(1 for r in results if r.status == 'indeterminate')

        self.run.info('Contigs analyzed', len(results))
        self.run.info('Predicted circular', circular_count)
        self.run.info('Predicted linear', linear_count)
        self.run.info('Indeterminate', indeterminate_count)

        if output_file:
            # generate an output dict dict
            output_dict = {}
            for r in results:
                output_dict[r.contig_name] = {'status': r.status,
                                              'circularity_support': f'{r.circularity_support:.4f}',
                                              'edge_coherence': f'{r.edge_coherence:.4f}',
                                              'contig_length': r.contig_length,
                                              'num_fr_pairs': r.num_fr_pairs,
                                              'num_edge_reads': r.num_edge_reads,
                                              'num_edge_reads_incoherent': r.num_edge_reads_incoherent,
                                              'observed_rf_pairs': r.observed_rf_pairs,
                                              'supporting_rf_pairs': r.supporting_rf_pairs,
                                              'expected_rf_pairs': f'{r.expected_rf_pairs:.2f}',
                                              'flags': ';'.join(r.flags) if r.flags else ''}


            # Store the output
            headers = ['contig', 'contig_length', 'status', 'circularity_support', 'edge_coherence',
                       'num_fr_pairs', 'expected_rf_pairs', 'observed_rf_pairs', 'supporting_rf_pairs',
                       'num_edge_reads', 'num_edge_reads_incoherent', 'flags']
            utils.store_dict_as_TAB_delimited_file(output_dict, output_file, headers=headers)
            self.run.info('Output file', output_file, nl_after=1)

        return results


    def _is_valid_fr_pair(self, read):
        """Check if read is the forward member of a valid FR pair."""
        return (not read.is_unmapped
                and not read.mate_is_unmapped
                and not read.is_secondary
                and not read.is_supplementary
                and read.is_proper_pair
                and not read.is_reverse          # This read is forward
                and read.mate_is_reverse         # Mate is reverse
                and read.template_length > 0)    # Forward read has positive TLEN


    def _is_valid_primary_read(self, read):
        """Check if read is a valid primary alignment with a mapped mate."""

        return (not read.is_unmapped and not read.mate_is_unmapped and not read.is_secondary and not read.is_supplementary)


    def _get_pair_orientation(self, read):
        """Determine the orientation of a read pair.

        Assumes `read` is the leftmost of the pair (by reference_start).

        Returns
        -------
        str
            'FR' if forward-reverse (normal), 'RF' if reverse-forward (junction-spanning),
            'FF' or 'RR' for same-strand pairs.
        """

        left_is_reverse = read.is_reverse
        right_is_reverse = read.mate_is_reverse

        if not left_is_reverse and right_is_reverse:
            return 'FR'
        elif left_is_reverse and not right_is_reverse:
            return 'RF'
        elif not left_is_reverse and not right_is_reverse:
            return 'FF'
        else:
            return 'RR'


# The functions below are helpers of the Read class which exist outside the class because they are
# just-in-time compiled (very very fast) with numba, which has poor support for in-class methods

@jit(nopython=True)
def iterate_cigartuples(cigartuples, cigar_consumption):
    """Iterate through cigartuples

    Parameters
    ==========
    cigartuples : Nx2 array

    Yields
    ======
    output : tuple
        (operation, length, consumes_read, consumes_ref) -> (int, int, bool, bool)
    """

    for i in range(cigartuples.shape[0]):
        operation, length = cigartuples[i, :]

        yield np.array([
            operation,
            length,
            cigar_consumption[operation, 0],
            cigar_consumption[operation, 1]
        ])


@jit(nopython=True)
def _vectorize_read(cigartuples, query_sequence, reference_sequence, reference_start, cigar_consumption):
    # init the array
    size = 0
    for i in range(cigartuples.shape[0]):
        size += cigartuples[i, 1]
    v = np.full((size, 4), -1, dtype=np.int32)

    count = 0
    ref_pos = 0       # position on the reference genome (for column 0)
    ref_seq_idx = 0   # index into reference_sequence array (for column 3)
    read_consumed = 0

    for i in range(cigartuples.shape[0]):
        operation, length = cigartuples[i, :]
        consumes_read, consumes_ref = cigar_consumption[operation, :]

        if consumes_read and consumes_ref:
            v[count:(count + length), 0] = np.arange(ref_pos + reference_start, ref_pos + reference_start + length)
            v[count:(count + length), 1] = query_sequence[read_consumed:(read_consumed + length)]
            v[count:(count + length), 3] = reference_sequence[ref_seq_idx:(ref_seq_idx + length)]
            v[count:(count + length), 2] = 0

            read_consumed += length
            ref_pos += length
            ref_seq_idx += length

        elif consumes_read:
            v[count:(count + length), 0] = ref_pos + reference_start - 1
            v[count:(count + length), 1] = query_sequence[read_consumed:(read_consumed + length)]
            v[count:(count + length), 2] = 1

            read_consumed += length

        elif consumes_ref:
            v[count:(count + length), 0] = np.arange(ref_pos + reference_start, ref_pos + reference_start + length)
            v[count:(count + length), 2] = 2

            # Only fetch reference bases for deletions (op=2), not skips (op=3)
            if operation == 2:  # Deletion
                v[count:(count + length), 3] = reference_sequence[ref_seq_idx:(ref_seq_idx + length)]
                ref_seq_idx += length
            # For N (op=3), leave column 3 as -1 (already initialized)

            ref_pos += length

        count += length

    return v


@jit(nopython=True)
def _get_percent_id(vectorized_read):
    # Filter vectorization to only include mapped segments, i.e. mapping_type == 0
    vectorized_read = vectorized_read[vectorized_read[:, 2] == 0]

    mapped_length = vectorized_read.shape[0]

    num_equal = 0
    for i in range(mapped_length):
        if vectorized_read[i, 1] == vectorized_read[i, 3]:
            num_equal += 1

    return 100 * num_equal / mapped_length


@jit(nopython=True)
def _get_aligned_sequence_and_reference_positions(cigartuples, query_sequence, reference_start, cigar_consumption):

    # get size of arrays to init
    size = 0
    for i in range(cigartuples.shape[0]):
        if cigar_consumption[cigartuples[i, 0], 0] and cigar_consumption[cigartuples[i, 0], 1]:
            size += cigartuples[i, 1]

    # init the arrays
    aligned_sequence = np.zeros(size, dtype=np.int64)
    reference_positions = np.zeros(size, dtype=np.int64)

    ref_consumed, read_consumed = 0, 0
    num_mapped = 0
    for i in range(cigartuples.shape[0]):
        operation, length = cigartuples[i, :]
        consumes_read, consumes_ref = cigar_consumption[operation, :]

        if consumes_read and consumes_ref:
            aligned_sequence[num_mapped:num_mapped+length] = query_sequence[read_consumed:(read_consumed + length)]
            reference_positions[num_mapped:num_mapped+length] = np.arange(ref_consumed + reference_start, ref_consumed + reference_start + length)

            num_mapped += length
            read_consumed += length
            ref_consumed += length

        elif consumes_ref:
            ref_consumed += length

        elif consumes_read:
            read_consumed += length

    return aligned_sequence, reference_positions


@jit(nopython=True)
def _trim(cigartuples, cigar_consumption, query_sequence, reference_sequence, reference_start, reference_end, trim_by, side):

    cigartuples = cigartuples[::-1, :] if side == 1 else cigartuples

    ref_positions_trimmed = 0
    ref_bases_trimmed = 0  # only counts bases actually in reference_sequence
    read_positions_trimmed = 0

    count = 0
    for i in range(cigartuples.shape[0]):
        operation, length = cigartuples[i, :]
        consumes_read, consumes_ref = cigar_consumption[operation, :]

        if consumes_ref:
            remaining = trim_by - ref_positions_trimmed

            if length > remaining:
                # the length of the operation exceeds the required trim amount. So we will
                # terminate this iteration. To trim the cigar tuple, we replace it with a
                # truncated length
                cigartuples[count, 1] = length - remaining
                ref_positions_trimmed += remaining
                # Only count bases for non-N operations (N is op 3)
                if operation != 3:
                    ref_bases_trimmed += remaining
                if consumes_read:
                    read_positions_trimmed += remaining
                break

            ref_positions_trimmed += length
            # Only count bases for non-N operations
            if operation != 3:
                ref_bases_trimmed += length
            if consumes_read:
                read_positions_trimmed += length

        elif consumes_read:
            # this operation does not consume ref so it does not contribute to trim_by
            # but if we met our ref trim target, we should stop
            if ref_positions_trimmed >= trim_by:
                break
            read_positions_trimmed += length

        count += 1

    cigartuples = cigartuples[count:, :]

    if side == 1:
        cigartuples = cigartuples[::-1]
        if read_positions_trimmed > 0:
            query_sequence = query_sequence[:-read_positions_trimmed]
        if ref_bases_trimmed > 0:
            reference_sequence = reference_sequence[:-ref_bases_trimmed]
        reference_end -= ref_positions_trimmed
    else:
        cigartuples = cigartuples
        if read_positions_trimmed > 0:
            query_sequence = query_sequence[read_positions_trimmed:]
        if ref_bases_trimmed > 0:
            reference_sequence = reference_sequence[ref_bases_trimmed:]
        reference_start += ref_positions_trimmed

    return cigartuples, query_sequence, reference_sequence, reference_start, reference_end


