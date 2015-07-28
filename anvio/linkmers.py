# -*- coding: utf-8
"""LinkMer reporting classes.

   The default client is `anvi-report-linkmers`"""

import sys
import pysam
import hashlib

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

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



