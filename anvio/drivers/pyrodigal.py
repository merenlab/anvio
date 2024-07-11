# coding: utf-8
"""Interface for gene calling that uses `pyrodigal`."""

import os
import argparse
import pyrodigal

import anvio
import anvio.fastalib as f
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Pyrodigal:
    def __init__(self, args=None, progress=progress, run=run):
        self.progress = progress
        self.run = run
        self.args = args
        A = lambda x: (args.__dict__[x] if x in args.__dict__ else None) if args else None
        self.pyrodigal_translation_table = A('pyrodigal_translation_table')
        self.pyrodigal_single_mode = A('pyrodigal_single_mode')
        self.num_threads = A('num_threads')

        self.run.info('Num threads for gene calling', self.num_threads)

        self.installed_version = pyrodigal.__version__


    def process(self, fasta_file_path, output_dir):
        """Take the fasta file, run pyrodigal on it, and make sense of the output

        Returns a gene calls dict, and amino acid sequences dict.
        """

        # Set up the output files.
        self.genes_in_contigs = os.path.join(output_dir, 'contigs.genes')
        self.amino_acid_sequences_in_contigs = os.path.join(output_dir, 'contigs.amino_acid_sequences')

        self.run.warning("Anvi'o will use 'pyrodigal' by XXX (doi:XXX), which uses the approach originally implemented by "
                         "Hyatt et al (doi:10.1186/1471-2105-11-119), to identify open reading frames in your data. When "
                         "you publish your findings, please do not forget to properly credit both work.", lc='green', header="CITATION")

        # let's learn the number of sequences we will work with early on and report
        num_sequences_in_fasta_file = utils.get_num_sequences_in_fasta(fasta_file_path)

        # some nice logs.
        self.run.warning('', header='Finding ORFs in contigs using pyrodigal', lc='green')
        self.run.info('Number of sequences', pp(num_sequences_in_fasta_file))
        self.run.info('Procedure', 'single' if self.pyrodigal_single_mode else 'meta')
        self.run.info('Genes', self.genes_in_contigs)
        self.run.info('Amino acid sequences', self.amino_acid_sequences_in_contigs)

        self.progress.new('Processing')
        self.progress.update(f"Identifying ORFs using {terminal.pluralize('thread', self.num_threads)}.")

        # key variables to fill in
        gene_calls_dict = {}
        amino_acid_sequences_dict = {}

        # something something

        return gene_calls_dict, amino_acid_sequences_dict
