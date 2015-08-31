# -*- coding: utf-8
"""LinkMer reporting classes.

   The default client is `anvi-report-linkmers`"""

import os
import sys
import pysam
import hashlib

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError

run = terminal.Run()
pp = terminal.pretty_print


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = ["Faruk Uzun"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class LinkMerDatum:
    def __init__(self, read_id):
        self.read_id = read_id
        self.read_hash = hashlib.sha224(read_id).hexdigest()
        self.contig_name = None
        self.pos_in_contig = None
        self.pos_in_read = None
        self.base = None


    def __str__(self):
        return 'Contig: %s, C_pos: %d, R_pos: %d, Base: %s, hash: %s' % (self.contig,
                                                                         self.pos_in_contig,
                                                                         self.pos_in_read,
                                                                         self.base,
                                                                         self.read_id)


class LinkMersData:
    def __init__(self, bam):
        self.bam = bam
        self.data = []


    def append(self, contig_name, positions):
        for pileupcolumn in self.bam.pileup(contig_name):
            if pileupcolumn.pos not in positions:
                continue

            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del:
                    L = LinkMerDatum(pileupread.alignment.qname)
                    L.contig_name = contig_name
                    L.pos_in_contig = pileupcolumn.pos
                    L.pos_in_read = pileupread.query_position
                    L.base = pileupread.alignment.seq[pileupread.query_position]
                    self.data.append(L)


class LinkMers:
    """This class handles an input BAM file, a list of contigs, and positions within them to
       report bases in reads that contribute to positions of interest following Chris Quince's
       suggestion. Each read is reported with a unique ID, therefore linkage informaiton can
       be followed."""
    def __init__(self, args = None):
        self.args = args
        self.input_file_path = None 
        self.contigs_and_positions = {}

        self.progress = terminal.Progress()
        self.run = terminal.Run(width=35)

        if args:
            filesnpaths.is_file_exists(args.input_file)
            self.input_file_path = args.input_file

            if args.list_contigs:
                self.list_contigs()
                sys.exit()

            filesnpaths.is_file_exists(args.contigs_and_positions)
            filesnpaths.is_file_tab_delimited(args.contigs_and_positions, expected_number_of_fields = 2)

            f = open(args.contigs_and_positions)
            for line in f.readlines():
                contig_name, positions = line.split('\t')

                try:
                    positions = [int(pos) for pos in positions.split(',')]
                except ValueError:
                    raise ConfigError, 'Positions for contig "%s" does not seem to be comma-separated integers...' % contig_name

                self.contigs_and_positions[contig_name] = set(positions)

        self.bam = None
        self.linkmers = None


    def process(self):
        self.sanity_check()

        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        try:
            self.bam = pysam.Samfile(self.input_file_path, 'rb')
        except ValueError as e:
            self.progress.end()
            raise ConfigError, 'Are you sure "%s" is a BAM file? Because samtools is not happy with it: """%s"""' % (self.input_file_path, e)
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lenghts = self.bam.lengths

        try:
            self.num_reads_mapped = self.bam.mapped
        except ValueError:
            raise ConfigError, "It seems the BAM file is not indexed. See 'anvi-init-bam' script."

        self.run.info('input_bam', self.input_file_path)
        self.run.info('total_reads_mapped', pp(int(self.num_reads_mapped)))
        self.run.info('num_contigs', pp(len(self.contig_names)))
        self.run.info('num_contigs_of_interest', pp(len(self.contigs_and_positions)))
        self.run.info('num_positions', pp(sum([len(p) for p in self.contigs_and_positions.values()])))

        indexes = [self.contig_names.index(r) for r in self.contigs_and_positions if r in self.contig_names]
        self.contig_names = [self.contig_names[i] for i in indexes]
        self.contig_lenghts = [self.contig_lenghts[i] for i in indexes]

        self.linkmers = LinkMersData(self.bam)
        for i in range(0, len(self.contig_names)):
            contig_name = self.contig_names[i]
            contig_positions = self.contigs_and_positions[contig_name]

            self.linkmers.append(contig_name, contig_positions)

        return self.linkmers.data


    def report(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        output_file = open(output_file_path, 'w')
        output_file.write('\t'.join(['entry_id', 'contig_name', 'pos_in_contig', 'pos_in_read', 'base', 'read_id']) + '\n')
        entry_id = 0
        for d in self.linkmers.data:
            entry_id += 1
            output_file.write('%.5d\t%s\t%d\t%d\t%s\t%s\n' % (entry_id, d.contig_name, d.pos_in_contig, \
                                                             d.pos_in_read, d.base, d.read_hash))
        output_file.close()

        self.run.info('output_file', output_file_path)


    def list_contigs(self):
        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = pysam.Samfile(self.input_file_path, 'rb')
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lenghts = self.bam.lengths

        for tpl in sorted(zip(self.contig_lenghts, self.contig_names), reverse = True):
            print '%-40s %s' % (tpl[1], pp(int(tpl[0])))


    def sanity_check(self):
        if not self.contigs_and_positions:
            raise ConfigError, "Contig names and positions dictionary is empty"


class GetReadsFromBAM:
    """This class takes contig names (or a collection ID and bins), and provides
       access short reads in one or more BAM iles mapping to those contigs
       (expanded in #173)."""

    def __init__(self, args = None, run = terminal.Run(), progress = terminal.Progress()):
        self.run = run
        self.progress = progress

        self.args = args

        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        self.input_bam_files = A('input_bams')
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.collection_id = A('collection_id')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.debug = A('debug')
        self.output_file_path = A('output_file')

        self.bins = set([])
        self.split_names_of_interest = set([])


    def init(self):
        self.sanity_check()

        self.run.info('Input BAM file(s)', ', '.join([os.path.basename(f) for f in self.input_bam_files]))

        d = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.bins = d.keys()

        for split_names in d.values():
            self.split_names_of_interest.update(split_names)

        self.run.info('Collection ID', self.collection_id)
        self.run.info('Bin(s)', ', '.join(self.bins))
        self.run.info('Number of splits', pp(len(self.split_names_of_interest)))


    def get_short_reads_for_splits_dict(self):
        short_reads_for_splits_dict = {}

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
            splits_order = contigs_involved[contig_id].keys()
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
            bam_file_name = '.'.join(os.path.basename(bam_file_path).split('.')[:-1])

            self.progress.update('Creating a dictionary of matching short reads in %s ...' % bam_file_name)

            bam_file = pysam.Samfile(bam_file_path, 'rb')
            for contig_id, start, stop in contig_start_stops:
                for entry in bam_file.fetch(contig_id, start, stop):
                    '''
                    here's what's available in the entry object:
                    
                    ['aend', 'alen', 'aligned_pairs', 'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples', 'compare',
                     'flag', 'get_aligned_pairs', 'get_blocks', 'get_overlap', 'get_reference_positions', 'get_tag',
                     'get_tags', 'has_tag', 'infer_query_length', 'inferred_length', 'is_duplicate', 'is_paired', 
                     'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary',
                     'is_unmapped', 'isize', 'mapping_quality', 'mapq', 'mate_is_reverse', 'mate_is_unmapped', 'mpos', 'mrnm',
                     'next_reference_id', 'next_reference_start', 'opt', 'overlap', 'pnext', 'pos', 'positions', 'qend', 
                     'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end', 'query_alignment_length',
                     'query_alignment_qualities', 'query_alignment_sequence', 'query_alignment_start', 'query_length',
                     'query_name', 'query_qualities', 'query_sequence', 'reference_end', 'reference_id', 'reference_length',
                     'reference_start', 'rlen', 'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags', 'template_length', 'tid', 'tlen']'''

                    # we are doing only for 'single reads', but I think this has to take into account the paired-end case as well.
                    short_reads_for_splits_dict['_'.join([contig_id, str(start), str(stop), entry.query_name, bam_file_name])] = entry.query_sequence

        self.progress.end()

        return short_reads_for_splits_dict


    def store_short_reads_for_splits(self):
        short_reds_for_splits_dict = self.get_short_reads_for_splits_dict()

        self.progress.new('Storing reads')
        self.progress.update('...')
        utils.store_dict_as_FASTA_file(short_reds_for_splits_dict, self.output_file_path)
        self.progress.end()

        self.run.info('Num reads stored', pp(len(short_reds_for_splits_dict)))
        self.run.info('FASTA output', self.output_file_path)


    def sanity_check(self):
        bad_bam_files = []
        for bam_file_path in self.input_bam_files:
            try:
                filesnpaths.is_file_exists(bam_file_path)
                pysam.Samfile(bam_file_path, 'rb')
            except ValueError as e:
                bad_bam_files.append(bam_file_path)
        if len(bad_bam_files):
            raise ConfigError, 'Samtools is not happy with some of your bam files. The following\
                                file(s) do not look like proper BAM files [ here is the actual\
                                error: "%s"]: %s.' % (e, ','.join(bad_bam_files))

        if not self.output_file_path:
            self.output_file_path = 'short_reads.fa'

        filesnpaths.is_output_file_writable(self.output_file_path)
