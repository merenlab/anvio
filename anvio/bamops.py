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


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"



class BAMFileObject(pysam.AlignmentFile):
    def __init__(self, *args):
        """A class that is essentially pysam.AlignmentFile, with some added bonuses

        This class inherits pysam.AlignmentFile and adds a little flair. Init such an object the
        way you would an AlignmentFile, i.e. bam = bamops.BAMFileObject(path_to_bam)
        """
        self.input_bam_path = args[0]
        filesnpaths.is_file_exists(self.input_bam_path)

        try:
            pysam.AlignmentFile.__init__(self)
        except ValueError as e:
            raise ConfigError('Are you sure "%s" is a BAM file? Because samtools is not happy with it: """%s"""' % (self.input_bam_path, e))

        try:
            self.mapped
        except ValueError:
            raise ConfigError("It seems the BAM file is not indexed. See 'anvi-init-bam' script.")


    def fetch_and_trim(self, contig_name, start, end, *args, **kwargs):
        """Returns an read iterator that trims overhanging reads

        Like pysam.AlignmeFile.fetch(), except trims reads that overhang the start and end of the
        defined region so that they fit inside the start and stop.
        """
        for read in self.fetch(contig_name, start, end, *args, **kwargs):
            if read.cigartuples is None:
                # This read has no associated cigar string. This either means it did not align but
                # is in the BAM file anyways, or the mapping software decided it did not want to
                # include a cigar string for this read.
                continue

            read = Read(read)

            if start - read.reference_start > 0:
                read.trim(trim_by=(start - read.reference_start), side='left')

            if read.reference_end - end > 0:
                read.trim(trim_by=(read.reference_end - end), side='right')

            yield read


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

        # See self.vectorize
        self.v = None


    def vectorize(self):
        """Set the self.v attribute to provide array-like access to the read"""

        self.v = _vectorize_read(
            self.cigartuples,
            self.query_sequence,
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

        ref, read, pos_ref, pos_read = [], [], 0, 0
        for _, length, consumes_read, consumes_ref in iterate_cigartuples(self.cigartuples, constants.cigar_consumption):
            if consumes_read:
                read.extend([chr(x) for x in self.query_sequence[pos_read:(pos_read + length)]])
                pos_read += length
            else:
                read.extend(['-'] * length)

            if consumes_ref:
                ref.extend(['X'] * length)
                pos_ref += length
            else:
                ref.extend(['-'] * length)

        lines = [
            '<%s.%s object at %s>' % (self.__class__.__module__, self.__class__.__name__, hex(id(self))),
            ' ├── start, end : [%s, %s)' % (self.reference_start, self.reference_end),
            ' ├── cigartuple : %s' % [tuple(row) for row in self.cigartuples],
            ' ├── read       : %s' % ''.join(read),
            ' └── reference  : %s' % ''.join(ref),
        ]

        return '\n'.join(lines)


    def trim(self, trim_by, side='left'):
        """Trims self.read by either the left or right

        Modifies the attributes:

            query_sequence
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
                self.reference_start += trim_by

            else:
                self.query_sequence = self.query_sequence[:-trim_by]
                self.reference_end -= trim_by

            return

        # We are here because the read was not a simple mapping. There are indels and so we need to
        # parse cigartuples. We delegate to a just-in-time compiled function for a 4X speed gain

        (self.cigartuples,
         self.query_sequence,
         self.reference_start,
         self.reference_end) = _trim(self.cigartuples,
                                     constants.cigar_consumption,
                                     self.query_sequence,
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

        self.routine_dict = {
            'accurate': self._accurate_routine,
        }


    def run(self, bam, contig_or_split, start=None, end=None, method='accurate', max_coverage=None, skip_coverage_stats=False, **kwargs):
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

        method : string
            How do you want to calculate? Options: see self.routine_dict

        skip_coverage_stats : bool, False
            Should the call to process_c be skipped?

        **kwargs : **kwargs
            kwargs are passed to the method chosen
        """

        # if there are defined start and ends we have to trim reads so their ranges fit inside self.c
        iterator = bam.fetch if (start is None and end is None) else bam.fetch_and_trim

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
            routine = self.routine_dict[method]
        except KeyError:
            raise ConfigError("Coverage :: %s is not a valid method." % method)

        self.c = routine(c, bam, contig_name, start, end, iterator, **kwargs)

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


    def _accurate_routine(self, c, bam, contig_name, start, end, iterator):
        """Routine that accounts for gaps in the alignment

        Notes
        =====
        - There used to be an '_approximate_routine', but its only negligibly faster
        - Should typically not be called explicitly. Use run instead
        - fancy indexing of reference_positions was also considered, but is much slower because it
          uses fancy-indexing
          https://jakevdp.github.io/PythonDataScienceHandbook/02.07-fancy-indexing.html:
        """
        for read in iterator(contig_name, start, end):
            if read.cigartuples is None:
                # This read has no associated cigar string. This either means it did not align but
                # is in the BAM file anyways, or the mapping software decided it did not want to
                # include a cigar string for this read.
                continue

            if len(read.cigartuples) == 1:
                c[read.reference_start:(read.reference_start + read.cigartuples[0][1])] += 1
            else:
                for start, end in read.get_blocks():
                    c[start:end] += 1

        return c


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
        self.input_bam_files = A('input_bams')
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.output_file_path = A('output_file')
        self.output_file_prefix = A('output_file_prefix')
        self.gzip = A('gzip_output')
        self.split_R1_and_R2 = A('split_R1_and_R2')

        self.bins = set([])
        self.split_names_of_interest = set([])

        self.initialized = False


    def init(self):
        utils.is_contigs_db(self.contigs_db_path)

        self.run.info('Input BAM file(s)', ', '.join([os.path.basename(f) for f in self.input_bam_files]))

        d = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.bins = list(d.keys())

        for split_names in list(d.values()):
            self.split_names_of_interest.update(split_names)

        self.run.info('Collection ID', self.collection_name)
        self.run.info('Bin(s)', ', '.join(self.bins))
        self.run.info('Number of splits', pp(len(self.split_names_of_interest)))

        self.initialized = True


    def get_short_reads_for_splits_dict(self):
        if not self.initialized:
            raise ConfigError('The `GetReadsFromBAM` class is not initialized :/ Ad hoc use of this class is '
                              'OK, but in that case you should set `self.initialized` to True, and provide '
                              'the split names of interest manually.')

        if not len(self.split_names_of_interest):
            raise ConfigError("The split names of interest set is empty. This should have never happened. Good "
                              "job.")

        short_reads_for_splits_dict = {}
        if self.split_R1_and_R2:
            short_reads_for_splits_dict['R1'] = {}
            short_reads_for_splits_dict['R2'] = {}
            short_reads_for_splits_dict['UNPAIRED'] = {}
        else:
            short_reads_for_splits_dict['all'] = {}

        self.progress.new('Accessing reads')
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

        # at this point contig_start_stops knows every contig we are interested in, and
        # their start and stop positions based on what split ids were requested. we
        # shall go through each bam file the user is interested, and get those short reads
        # that map to regions of interest:
        for bam_file_path in self.input_bam_files:
            bam_file_name = filesnpaths.get_name_from_file_path(bam_file_path)

            bam_file_object = BAMFileObject(bam_file_path)

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

            has_unknown_mate = {}
            if self.split_R1_and_R2:
                for contig_id, start, stop in contig_start_stops:
                    for read in bam_file_object.fetch(contig_id, start, stop):

                        defline = '_'.join([contig_id, str(start), str(stop), read.query_name, bam_file_name])

                        if not read.is_paired:
                            short_reads_for_splits_dict['UNPAIRED'][defline] = read.query_sequence

                        elif defline in has_unknown_mate:
                            # `read`s mate has already been read. so assign the read and the mate
                            # to their respective 'R1' and 'R2' dictionaries, then remove the mate
                            # from has_unknown_mate since its mate is now known.
                            read_DIRECTION = 'R1' if read.is_read1 else 'R2'
                            mate_DIRECTION = 'R2' if read_DIRECTION == 'R1' else 'R1'
                            short_reads_for_splits_dict[mate_DIRECTION][defline] = has_unknown_mate[defline]
                            short_reads_for_splits_dict[read_DIRECTION][defline] = read.query_sequence
                            del has_unknown_mate[defline]

                        else:
                            has_unknown_mate[defline] = read.query_sequence
                short_reads_for_splits_dict['UNPAIRED'].update(has_unknown_mate)
            else:
                for contig_id, start, stop in contig_start_stops:
                    for read in bam_file_object.fetch(contig_id, start, stop):
                        short_reads_for_splits_dict['all']['_'.join([contig_id, str(start), str(stop), read.query_name, bam_file_name])] = read.query_sequence
            bam_file_object.close()

        self.progress.end()

        return short_reads_for_splits_dict


    def store_short_reads_for_splits(self):
        self.sanity_check()

        if not self.sanity_checked:
            raise ConfigError("store_short_reads_for_splits :: Cannot be called before running sanity_check")

        short_reds_for_splits_dict = self.get_short_reads_for_splits_dict()

        self.progress.new("Storing reads")
        self.progress.update("...")

        if self.split_R1_and_R2:
            for read_type in sorted(list(short_reds_for_splits_dict.keys())):
                output_file_path = '%s_%s.fa' % (self.output_file_prefix, read_type)

                utils.store_dict_as_FASTA_file(short_reds_for_splits_dict[read_type], output_file_path)
                if self.gzip:
                    utils.gzip_compress_file(output_file_path)
                    output_file_path = output_file_path + ".gz"

                self.run.info('Output file for %s' % read_type, output_file_path, progress=self.progress)

            self.progress.end()
            self.run.info('Num paired-end reads stored',pp(len(short_reds_for_splits_dict['R1'])), mc='green', nl_before=1)
            self.run.info('Num unpaired reads stored',pp(len(short_reds_for_splits_dict['UNPAIRED'])), mc='green')
        else:
            output_file_path = self.output_file_path or 'short_reads.fa'
            utils.store_dict_as_FASTA_file(short_reds_for_splits_dict['all'], output_file_path)

            if self.gzip:
                utils.gzip_compress_file(output_file_path)
                output_file_path = output_file_path + ".gz"

            self.progress.end()
            self.run.info('Output file for all short reads',output_file_path)
            self.run.info('Num reads stored', pp(len(short_reds_for_splits_dict['all'])), mc='green')


    def sanity_check(self):
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

        self.sanity_checked = True


class ReadsMappingToARange:
    """Returns all reads from BAM that maps to a range in a contig"""

    def __init__(self, run=run, progress=progress):
        self.data = []

        self.run = run
        self.progress = progress


    def process_range(self, input_bam_paths, contig_name, start, end):
        if end <= start:
            raise ConfigError("The end of range cannot be equal or smaller than the start of it :/")

        data = []

        for input_bam_path in input_bam_paths:
            bam_file_object = BAMFileObject(input_bam_path)

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


# The below functions are helpers of the Read class which exist outside the class because they are
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
def _vectorize_read(cigartuples, query_sequence, reference_start, cigar_consumption):
    # init the array
    size = 0
    for i in range(cigartuples.shape[0]):
        size += cigartuples[i, 1]
    v = np.full((size, 3), -1, dtype=np.int32)

    count = 0
    ref_consumed = 0
    read_consumed = 0
    for i in range(cigartuples.shape[0]):
        operation, length = cigartuples[i, :]
        consumes_read, consumes_ref = cigar_consumption[operation, :]

        if consumes_read and consumes_ref:
            v[count:(count + length), 0] = np.arange(ref_consumed + reference_start, ref_consumed + reference_start + length)
            v[count:(count + length), 1] = query_sequence[read_consumed:(read_consumed + length)]
            v[count:(count + length), 2] = 0

            read_consumed += length
            ref_consumed += length

        elif consumes_read:
            v[count:(count + length), 1] = query_sequence[read_consumed:(read_consumed + length)]
            v[count:(count + length), 2] = 1

            read_consumed += length

        elif consumes_ref:
            v[count:(count + length), 0] = np.arange(ref_consumed + reference_start, ref_consumed + reference_start + length)
            v[count:(count + length), 2] = 2

            ref_consumed += length

        count += length

    return v


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
def _trim(cigartuples, cigar_consumption, query_sequence, reference_start, reference_end, trim_by, side):

    cigartuples = cigartuples[::-1, :] if side == 1 else cigartuples

    ref_positions_trimmed = 0
    read_positions_trimmed = 0
    terminate_next = False

    count = 0
    for i in range(cigartuples.shape[0]):
        operation, length = cigartuples[i, :]
        consumes_read, consumes_ref = cigar_consumption[operation, :]

        if consumes_ref and consumes_read:
            if terminate_next:
                break

            remaining = trim_by - ref_positions_trimmed

            if length > remaining:
                # the length of the operation exceeds the required trim amount. So we will
                # terminate this iteration. To trim the cigar tuple, we replace it with a
                # truncated length
                cigartuples[count, 1] = length - remaining
                ref_positions_trimmed += remaining
                read_positions_trimmed += remaining
                break

            ref_positions_trimmed += length
            read_positions_trimmed += length

        elif consumes_ref:
            ref_positions_trimmed += length

        elif consumes_read:
            read_positions_trimmed += length

        if ref_positions_trimmed >= trim_by:
            terminate_next = True

        count += 1

    cigartuples = cigartuples[count:, :]

    if side == 1:
        cigartuples = cigartuples[::-1]
        query_sequence = query_sequence[:-read_positions_trimmed]
        reference_end -= ref_positions_trimmed
    else:
        cigartuples = cigartuples
        query_sequence = query_sequence[read_positions_trimmed:]
        reference_start += ref_positions_trimmed

    return cigartuples, query_sequence, reference_start, reference_end


