# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes to define and work with anvi'o SRA_downloads workflow. """

import os
import anvio
import pandas as pd

import anvio.utils as utils
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

        # check that NCBI SRA Toolkit and other programs are installed
        NCBI_sra_tool_programs = ['prefetch', 'fasterq-dump']
        other_programs = ['pigz']

        for program in NCBI_sra_tool_programs:
            if not utils.is_program_exists(program, dont_raise=True):
                raise ConfigError(f"The program {program} is not installed in your anvi'o conda environment. "
                                  f"'prefetch' and 'fasterq-dump'  are from the NCBI SRA toolkit and must be installed for the "
                                  f"sra_download workflow to work. Please check out the installation instructions here: "
                                  f"https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit")
        for program in other_programs:
            if not utils.is_program_exists(program, dont_raise=True):
                raise ConfigError(f"The program {program} is not installed in your anvi'o conda environment. Please "
                                  f"double check you installed all of the programs listed in the anvio'o installation tutorial: https://anvio.org/install/")

        # Snakemake rules
        self.rules.extend(['prefetch',
                           'fasterq_dump',
                           'pigz'])

        self.general_params.extend(['SRA_accession_list']) # user needs provide a tsv of SRA accessions

        # Parameters for each rule that are accessible in the config.json file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['prefetch'] = ['--max-size']
        rule_acceptable_params_dict['fasterq-dump'] = ['--split-files', ' --verbose', '--progress']

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
                raise ConfigError(f"Looks like your SRA accession list file, {self.SRA_accession_list}, is not properly formatted. "
                                  f"This is what we know: {e}")

        for accession in self.accessions_list:
            if not accession.startswith(('SRR', 'ERR', 'DRR')):
                if accession.startswith('SAMEA'):
                    raise ConfigError(f"anvi'o found an NCBI BioSample in your {self.SRA_accession_list}: {accession}. "
                                      f"The anvi'o sra-download workflow only processes sequencing accessions that start with the prefix: ERR, SRR, or DRR. "
                                      f"Search for the BioSample accession '{accession}' on the [NCBI SRA website](https://www.ncbi.nlm.nih.gov/sra) "
                                      f"and find the sequencing accessions.")
                else:
                    raise ConfigError(f"Looks like one of your \"SRA accessions\", {accession}, is not an SRA accession :( "
                                      f"anvi'o asks that you kindly double check your SRA_accession_list.txt ({self.SRA_accession_list}) to confirm you "
                                      f"are using the correct accessions. Hint: SRA accessions start with the prefix: ERR, SRR, or DRR")

        self.target_files = self.get_target_files()


    def get_target_files(self):
        """Get list of target files for snakemake target rule"""

        target_files = [os.path.join(self.dirs_dict['FASTAS'], f"generate_samples_txt.done")]

        return target_files
