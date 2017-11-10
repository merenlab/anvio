# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

import os
import sys
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.genecalling as genecalling

from anvio.errors import ConfigError

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class ContigsProfiler(object):
    def __init__(self, args, run=run, progress=progress):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress

        self.output_db_path = A('output_db_path')
        self.contigs_fasta = os.path.abspath(A('contigs_fasta'))
        self.project_name = A('project_name')
        self.description_file_path = A('description')
        self.split_length = A('split_length')
        self.kmer_size = A('kmer_size')
        self.skip_gene_calling = A('skip_gene_calling')
        self.external_gene_calls = A('external_gene_calls')
        self.skip_mindful_splitting = A('skip_mindful_splitting')
        self.ignore_internal_stop_codons = A('ignore_internal_stop_codons')
        self.debug = A('debug')

        self.process()


    def sanity_check(self):
        if self.external_gene_calls:
            filesnpaths.is_file_exists(self.external_gene_calls)

        if self.external_gene_calls and self.skip_gene_calling:
            raise ConfigError("You provided a file for external gene calls, and used requested gene calling to be\
                                skipped. Please make up your mind.")

        if not self.project_name:
            raise ConfigError("Sorry, you must provide a project name for your contigs database :/")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'rU').read()
        else:
            self.description = ''

        self.check_fasta()

        if not self.split_length:
            raise ConfigError("Creating a new contigs database requires split length information to be\
                                provided. But the ContigsDatabase class was called to create one without this\
                                bit of information. Not cool.")

        try:
            self.split_length = int(self.split_length)
        except:
            raise ConfigError("Split size must be an integer.")

        if self.split_length <= 0:
            self.split_length = sys.maxsize

        try:
            self.kmer_size = int(self.kmer_size)
        except:
            raise ConfigError("K-mer size must be an integer.")
        if self.kmer_size < 2 or self.kmer_size > 8:
            raise ConfigError("We like our k-mer sizes between 2 and 8, sorry! (but then you can always change the\
                                source code if you are not happy to be told what you can't do, let us know how it goes!).")

        if self.skip_gene_calling:
            self.skip_mindful_splitting = True


    def check_fasta(self):
        if not os.path.exists(self.contigs_fasta):
            raise ConfigError("Creating a new contigs database requires a FASTA file with contigs to be provided.")

        filesnpaths.is_file_fasta_formatted(self.contigs_fasta)

        # go throught he FASTA file to make sure there are no surprises with deflines and sequence lengths.
        self.progress.new('Checking deflines and contig lengths')
        self.progress.update('tick tock ...')
        fasta = u.SequenceSource(self.contigs_fasta)
        while next(fasta):
            if not utils.check_contig_names(fasta.id, dont_raise=True):
                self.progress.end()
                raise ConfigError("At least one of the deflines in your FASTA File does not comply with the 'simple deflines'\
                                    requirement of anvi'o. You can either use the script `anvi-script-reformat-fasta` to take\
                                    care of this issue, or read this section in the tutorial to understand the reason behind\
                                    this requirement (anvi'o is very upset for making you do this): %s" % \
                                        ('http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file'))

            if len(fasta.seq) < self.kmer_size:
                self.progress.end()
                raise ConfigError("At least one of the contigs in your input FASTA '%s' is shorter than the k-mer size. The k\
                                    is %d, and your contig is like %d :/ Anvi'o will not judge you for whatever you are doing\
                                    with such short contigs, but the length of each contig must be at least as long as your `k` for\
                                    k-mer analyis. You can use the script `anvi-script-reformat-fasta` to get rid of very short\
                                    contigs if you like." % (contigs_fasta, kmer_size, len(fasta.seq)))
        fasta.close()
        self.progress.end()

        all_ids_in_FASTA = utils.get_all_ids_from_fasta(self.contigs_fasta)
        if len(all_ids_in_FASTA) != len(set(all_ids_in_FASTA)):
            raise ConfigError("Every contig in the input FASTA file must have a unique ID. You know...")


    def print_arguments_summary(self):
        self.run.info('Name', self.project_name, mc='green')
        self.run.info('Description', os.path.abspath(self.description_file_path) if self.description_file_path else 'No description is given', mc='green')
        self.run.info('Input FASTA file', self.contigs_fasta)
        self.run.info('Split Length', pp(self.split_length))
        self.run.info('K-mer size', self.kmer_size)
        self.run.info('Skip gene calling?', self.skip_gene_calling)
        self.run.info('External gene calls provided?', self.external_gene_calls)
        self.run.info('Ignoring internal stop codons?', self.ignore_internal_stop_codons)
        self.run.info('Splitting pays attention to gene calls?', (not self.skip_mindful_splitting))


    def print_final_report(self):
        self.run.info('Contigs database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of contigs', contigs_info_table.total_contigs, quiet=self.quiet)
        self.run.info('Number of splits', splits_info_table.total_splits, quiet=self.quiet)
        self.run.info('Total number of nucleotides', contigs_info_table.total_nts, quiet=self.quiet)
        self.run.info('Gene calling step skipped', skip_gene_calling, quiet=self.quiet)
        self.run.info("Splits broke genes (non-mindful mode)", skip_mindful_splitting, quiet=self.quiet)
        self.run.info('Desired split length (what the user wanted)', split_length, quiet=self.quiet)
        self.run.info("Average split length (wnat anvi'o gave back)", (int(round(numpy.mean(recovered_split_lengths)))) \
                                                                        if recovered_split_lengths \
                                                                            else "(Anvi'o did not create any splits)", quiet=self.quiet)


    def do_gene_calling(self):
        genes_in_contigs_dict = {}
        contig_name_to_gene_start_stops = {}
        if not self.skip_gene_calling:
            #gene_calls_tables = TablesForGeneCalls(self.db_path, contigs_fasta, debug=debug)
            gene_caller = None
            if self.external_gene_calls:
                gene_caller = genecalling.ExternalGeneCaller(self.contigs_fasta, self.external_gene_calls, ignore_internal_stop_codons=self.ignore_internal_stop_codons)
            else:
                gene_caller = genecalling.GeneCaller(self.contigs_fasta, gene_caller=gene_caller, debug=self.debug)
            
            gene_calls_dict, protein_sequences = gene_caller.process()

            self.contigs_db.gene_calls_tables.store(gene_calls_dict, protein_sequences)

            # reconnect and learn about what's done
            #self.db = db.DB(self.db_path, anvio.__contigs__version__)
            #genes_in_contigs_dict = self.db.get_table_as_dict(t.genesx_in_contigs_table_name)

            for gene_unique_id in genes_in_contigs_dict:
                e = genes_in_contigs_dict[gene_unique_id]
                if e['contig'] not in contig_name_to_gene_start_stops:
                    contig_name_to_gene_start_stops[e['contig']] = set([])

                contig_name_to_gene_start_stops[e['contig']].add((gene_unique_id, e['start'], e['stop']), )


    def create_new_contigs_db(self):
        meta_values = {
            'split_length': self.split_length,
            'project_name': self.project_name,
            'description': self.description,
            'kmer_size': self.kmer_size
        }

        self.contigs_db = dbops.ContigsDatabase(self.output_db_path, create_new=True, meta_values=meta_values)

    def process(self):
        self.sanity_check()

        self.print_arguments_summary()

        self.create_new_contigs_db()

        # first things first: do the gene calling on contigs. this part is important. we are doing the
        # gene calling first. so we understand where genes start and end. this information will guide the
        # arrangement of the breakpoint of splits
        self.do_gene_calling()

        self.print_final_report()

