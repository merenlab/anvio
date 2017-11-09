# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

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
        self.args = args
        self.run = run
        self.progress = progress

        self.contigs_fasta = A('contigs_fasta')
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


    def process(self):
        pass
