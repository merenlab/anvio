# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for gene calling.
"""

import shutil

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

from anvio.drivers.prodigal import Prodigal


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class GeneCaller:
    def __init__(self, fasta_file_path, gene_caller='prodigal', args=None, progress=progress, run=run, debug=False):
        filesnpaths.is_file_exists(fasta_file_path)
        filesnpaths.is_file_fasta_formatted(fasta_file_path)

        self.fasta_file_path = fasta_file_path

        self.run = run
        self.progress = progress

        self.args = args
        self.debug = debug
        self.tmp_dirs = []

        self.gene_callers = {'prodigal': Prodigal}

        self.gene_caller = gene_caller

        if self.gene_caller not in self.gene_callers:
            raise ConfigError("The gene caller you requested ('%s') is not available at this point. "
                               "here is a list of what we have: %s." % (', '.join(self.gene_callers)))


    def process(self):
        output_dir = filesnpaths.get_temp_directory_path()
        self.tmp_dirs.append(output_dir)
        gene_caller = self.gene_callers[self.gene_caller](args=self.args)

        gene_calls_dict, amino_acid_sequences_dict = gene_caller.process(self.fasta_file_path, output_dir)

        if not self.debug:
            self.clean_tmp_dirs()

        return gene_calls_dict, amino_acid_sequences_dict


    def clean_tmp_dirs(self):
        for tmp_dir in self.tmp_dirs:
            shutil.rmtree(tmp_dir)
