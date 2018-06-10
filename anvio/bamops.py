# -*- coding: utf-8
# pylint: disable=line-too-long
"""LinkMer reporting classes.

   The default client is `anvi-report-linkmers`"""



import os
import sys
import pysam
import hashlib
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



class CodonFrequencies:
    def __init__(self, run=run):
        self.run = run
        self.not_reported_items = Counter({})


    def get_codon_frequencies_dict(self, codon_frequencies):
        pass


    def process_gene_call(self, bam_file_object, gene_call, contig_sequence, codons_to_profile=None, return_AA_frequencies_instead=False):
        if gene_call['partial']:
            return None

        contig_name = gene_call['contig']

        # here we will create a dictionary to translate codons in a gene to nucleotide positions in the context
        # of the contig. thanks to this dict, we will be able to profile only a small number of codons from a
        # gene, if they are specified in `codons_to_profile` variable. for instance, this function is called
        # during profiling only with codon positions in genes that possess nucleotide variation.
        codon_order_to_nt_positions = utils.get_codon_order_to_nt_positions_dict(gene_call)

        # here we generate the actual 'linkmers' information.
        d = {}
        linkmers = LinkMersData()
        linkmers.quiet = True
        for codon_order in codon_order_to_nt_positions:
            if codons_to_profile and codon_order not in codons_to_profile:
                continue

            nt_positions = codon_order_to_nt_positions[codon_order]

            reference_codon_sequence = contig_sequence[nt_positions[0]:nt_positions[2] + 1]

            # if consensus sequence contains shitty characters, we will not continue
            if reference_codon_sequence not in constants.codon_to_AA:
                continue

            linkmers.data = []
            linkmers.append(bam_file_object, 'sample_id', None, contig_name, nt_positions, only_complete_links=True)
            data = linkmers.data[0][2]

            hash_to_oligotype = {}
            unique_hashes = set([datum.read_unique_id for datum in data])
            for unique_hash in unique_hashes:
                hash_to_oligotype[unique_hash] = []

            for datum in data:
                hash_to_oligotype[datum.read_unique_id].append((datum.pos_in_contig, datum.base),)

            for unique_hash in unique_hashes:
                hash_to_oligotype[unique_hash] = ''.join([e[1] for e in sorted(hash_to_oligotype[unique_hash])])

            codon_frequencies = Counter(list(hash_to_oligotype.values()))

            # depending on what we want to return item frequencies will contain frequencies for amino acids or
            # codons.
            item_frequencies = Counter({})

            # our conversion dicts will differ if the user is asking for codons or AAs, and if the gene is
            # reverse or forward.
            conv_dict = None
            gene_is_reverse = gene_call['direction'] == 'r'
            if return_AA_frequencies_instead:
                if gene_is_reverse:
                    conv_dict = constants.codon_to_AA_RC
                else:
                    conv_dict = constants.codon_to_AA
            else:
                if gene_is_reverse:
                    conv_dict = constants.codon_to_codon_RC
                else:
                    # no conversion is necessary, so this is a mock dictionary that
                    # resturns the key.
                    conv_dict = dict(zip(constants.codons, constants.codons))

            # the magic happens here:
            for codon in codon_frequencies:
                if codon in conv_dict:
                    item_frequencies[conv_dict[codon]] += codon_frequencies[codon]
                else:
                    # so there is a way the programmer can learn that some weird
                    # stuff did not get reported
                    self.not_reported_items[codon] += 1

            reference_item = conv_dict[reference_codon_sequence]
            coverage = sum(item_frequencies.values())

            if not coverage:
                # FIXME: there was at least one case where the coverage here in this context was 0,
                #        which crashed the profiling. we never went after this issue, and it is
                #        important to understand how often this happens, and why.
                continue

            # here we quantify the ratio of frequencies of non-reference-aas observed in this codon
            # to the overall overage, and that is our `departure_from_reference`:
            total_frequency_of_all_items_but_the_conensus = sum([item_frequencies[item] for item in item_frequencies if item != reference_item])
            departure_from_reference = total_frequency_of_all_items_but_the_conensus / coverage

            d[codon_order] = {'reference': reference_item,
                              'coverage': coverage,
                              'frequencies': item_frequencies,
                              'departure_from_reference': departure_from_reference}

        return d


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


class BAMFileObject:
    def __init__(self, input_bam_path, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.input_bam_path = input_bam_path

        filesnpaths.is_file_exists(input_bam_path)

    def get(self):
        try:
            bam_file_object = pysam.Samfile(self.input_bam_path, 'rb')
        except ValueError as e:
            raise ConfigError('Are you sure "%s" is a BAM file? Because samtools is not happy with it: """%s"""' % (self.input_bam_path, e))

        try:
            bam_file_object.mapped
        except ValueError:
            self.progress.end()
            raise ConfigError("It seems the BAM file is not indexed. See 'anvi-init-bam' script.")

        return bam_file_object


class LinkMers:
    """This class handles an input BAM file, a list of contigs, and positions within them to
       report bases in reads that contribute to positions of interest following Chris Quince's
       suggestion. Each read is reported with a unique ID, therefore linkage informaiton can
       be followed."""
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
            bam_file_object = BAMFileObject(input_file).get()

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
        bam_file_object = BAMFileObject(self.input_file_path).get()
        self.progress.end()

        contig_names = self.bam.references
        contig_lengths = self.bam.lengths

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
            raise ConfigError('The `GetReadsFromBAM` class is not initialized :/ Ad hoc use of this class is\
                               OK, but in that case you should set `self.initialized` to True, and provide\
                               the split names of interest manually.')

        if not len(self.split_names_of_interest):
            raise ConfigError("The split names of interest set is empty. This should have never happened. Good\
                               job.")

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

            bam_file_object = BAMFileObject(bam_file_path).get()

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
                bam_file_object = BAMFileObject(bam_file_path).get()
                bam_file_object.close()
            except ConfigError as e:
                bad_bam_files.append(bam_file_path)
                error_message = e

        if len(bad_bam_files):
            raise ConfigError('Samtools is not happy with some of your bam files. The following\
                               file(s) do not look like proper BAM files [here is the actual\
                               error: "%s"]: %s.' % (error_message, ','.join(bad_bam_files)))

        if self.output_file_prefix and self.output_file_path:
            raise ConfigError("You must either use the parameter output file name, or output file prefix.")

        if self.output_file_prefix and not self.split_R1_and_R2:
            raise ConfigError("Output file prefix parameter is only relevant when you want to split R1 reads\
                               from R2 reads and so on.")

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
            bam_file_object = BAMFileObject(input_bam_path).get()

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



