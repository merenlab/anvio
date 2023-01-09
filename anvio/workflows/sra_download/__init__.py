# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes to define and work with anvi'o SRA_downloads workflow. """

import os
import anvio
import pandas as pd

import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = ['mschecht']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


run = terminal.Run()

class SRADownloadWorkflow(WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='sra_download')

        # Snakemake rules
        self.rules.extend(['prefetch',
                           'fasterq_dump',
                           'pigz'])

        self.general_params.extend(['SRA_accession_list']) # user needs provide a tsv of SRA accessions

        # Parameters for each rule that are accessible in the config.json file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['prefetch'] = ['--max-size']
        rule_acceptable_params_dict['fasterq-dump'] = ['--split-files', ' --verbose', '--progress']
        rule_acceptable_params_dict['pigz'] = ['--processes']

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Set default values for accessible rules and order of rules in config.json file
        self.default_config.update({
            'SRA_accession_list': 'SRA_accession_list.txt',
            'prefetch': {'--max-size': "40g", 'threads': 2},
            'fasterq_dump': {'threads': 6},
            'pigz': {'threads': 8}})

        # Directory structure for Snakemake workflow
        self.dirs_dict.update({"SRA_prefetch": "01_NCBI_SRA"})
        self.dirs_dict.update({"FASTAS": "02_FASTA"})


    def init(self):
        """Load the SRA_accession_list and creates a list of target files for the Snakemake workflow."""
        
        super().init()

        # Load SRA_accession_list
        self.SRA_accession_list = self.get_param_value_from_config(['SRA_accession_list'])

        if not self.SRA_accession_list:
            raise ConfigError('Please provide a list of SRA accessions')

        if self.SRA_accession_list:
            filesnpaths.is_file_exists(self.SRA_accession_list)
            try:
                self.SRA_accession_list_df = pd.read_csv(self.SRA_accession_list, sep='\t', index_col=False, header=None, names=['accessions'])
                self.accessions_list = self.SRA_accession_list_df['accessions'].tolist()
            except IndexError as e:
                raise ConfigError(f"Looks like your samples_txt file, {self.SRA_accession_list}, is not properly formatted. "
                                  f"This is what we know: {e}")

        self.target_files = self.get_target_files()


    def get_target_files(self):
        """Get list of target files for snakemake target rule"""

        target_files = []

        for accession in self.accessions_list:
            target_files.extend([os.path.join(self.dirs_dict['FASTAS'], f"{accession}_1.fastq.gz"),
                                 os.path.join(self.dirs_dict['FASTAS'], f"{accession}_2.fastq.gz")])
        
        return target_files