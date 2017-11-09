# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

import os
import sys
import anvio.utils as utils
import anvio.fastalib as u
import anvio.filesnpaths as filesnpaths

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

        all_ids_in_FASTA = utils.get_all_ids_from_fasta(contigs_fasta)
        if len(all_ids_in_FASTA) != len(set(all_ids_in_FASTA)):
            raise ConfigError("Every contig in the input FASTA file must have a unique ID. You know...")

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
        fasta = u.SequenceSource(contigs_fasta)
        while next(fasta):
            if not utils.check_contig_names(fasta.id, dont_raise=True):
                self.progress.end()
                raise ConfigError("At least one of the deflines in your FASTA File does not comply with the 'simple deflines'\
                                    requirement of anvi'o. You can either use the script `anvi-script-reformat-fasta` to take\
                                    care of this issue, or read this section in the tutorial to understand the reason behind\
                                    this requirement (anvi'o is very upset for making you do this): %s" % \
                                        ('http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file'))

            if len(fasta.seq) < kmer_size:
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



    def process(self):
        pass
