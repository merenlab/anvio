# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import os
import anvio
import anvio.utils as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.workflows import init_workflow_super_class
from anvio.workflows import WorkflowSuperClass


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


class ContigsDBWorkflow(WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        self.group_names = []
        self.contigs_information = {}
        self.fasta_information = {}

        # initialize the base class
        init_workflow_super_class(self, args, workflow_name='contigs')

        self.rules.extend(['anvi_script_reformat_fasta',
                           'anvi_gen_contigs_database', 'export_gene_calls_for_centrifuge', 'centrifuge',
                           'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database', 'anvi_get_sequences_for_gene_calls',
                           'emapper', 'anvi_script_run_eggnog_mapper', 'gunzip_fasta'])

        self.general_params.extend(["fasta_txt"])

        self.dirs_dict.update({"FASTA_DIR": "01_FASTA",
                               "CONTIGS_DIR": "02_CONTIGS"})

        self.default_config.update({"fasta_txt": "fasta.txt",
                                    "anvi_gen_contigs_database": {"--project-name": "{group}"},
                                    "centrifuge": {"threads": 2},
                                    "anvi_run_hmms": {"run": True, "threads": 5},
                                    "anvi_run_ncbi_cogs": {"run": True, "threads": 5},
                                    "anvi_script_reformat_fasta": {"run": True, "--simplify-names": True},
                                    "emapper": {"--database": "bact", "--usemem": True, "--override": True},
                                    "anvi_script_run_eggnog_mapper": {"--use-version": "0.12.6"}})

        self.rule_acceptable_params_dict['anvi_run_ncbi_cogs'] = ['run', '--cog-data-dir', '--sensitive', '--temporary-dir-path', '--search-with']

        self.rule_acceptable_params_dict['anvi_run_hmms'] = ['run', '--installed-hmm-profile', '--hmm-profile-dir']

        self.rule_acceptable_params_dict['centrifuge'] = ['run', 'db']

        self.rule_acceptable_params_dict['emapper'] = ['--database', '--usemem', '--override', 'path_to_emapper_dir']

        self.rule_acceptable_params_dict['anvi_script_run_eggnog_mapper'] = ['run', '--cog-data-dir', '--drop-previous-annotations',
                                         '--use-version']

        self.rule_acceptable_params_dict['anvi_script_reformat_fasta'] = \
                    ['run', '--simplify-names', '--keep-ids', '--exclude-ids', '--min-len']


        gen_contigs_params = ['--description', '--skip-gene-calling', '--external-gene-calls',\
                              '--ignore-internal-stop-codons', '--skip-mindful-splitting',\
                              '--contigs-fasta', '--project-name',\
                              '--description', '--split-length', '--kmer-size',\
                              '--skip-mindful-splitting', '--skip-gene-calling', '--external-gene-calls',\
                              '--ignore-internal-stop-codons']

        self.rule_acceptable_params_dict['anvi_gen_contigs_database'] = gen_contigs_params


    def init(self):
        super().init()

        fasta_txt_file = self.get_param_value_from_config('fasta_txt', repress_default=True)

        if fasta_txt_file:
            filesnpaths.is_file_exists(fasta_txt_file)
            self.contigs_information = u.get_TAB_delimited_file_as_dictionary(fasta_txt_file)
            self.fasta_information.update(self.contigs_information)
            self.group_names = list(self.contigs_information.keys())
            self.references_mode = True


    def get_raw_fasta(self, wildcards, remove_gz_suffix=True):
        '''
            Define the path to the input fasta files.
        '''
        contigs = self.fasta_information[wildcards.group]['path']
        ends_with_gz = contigs.endswith('.gz')
        if remove_gz_suffix and ends_with_gz:
            # we need to gunzip the fasta file
            # we will create a temporary uncompressed fasta file.
            contigs = os.path.join(self.dirs_dict['FASTA_DIR'], \
                                   wildcards.group + '-temp.fa')
        return contigs


    def get_fasta(self, wildcards):
        '''
            Define the path to the input fasta files.
        '''
        # The raw fasta will be used if no formatting is needed
        contigs = self.get_raw_fasta(wildcards)

        if self.get_param_value_from_config(['anvi_script_reformat_fasta','run']):
            # by default, reformat fasta is ran
            contigs = self.dirs_dict["FASTA_DIR"] + "/{group}/{group}-contigs.fa".format(group=wildcards.group)

        return contigs
