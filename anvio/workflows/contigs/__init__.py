# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import anvio
import anvio.terminal as terminal

from anvio.workflows import WorkflowSuperClass
from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


class ContigsDBWorkflow(WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        # if a regular instance of `ContigsDBWorkflow` is being generated, we
        # expect it to have a parameter `args`. if there is no `args` given, we
        # assume the class is being inherited as a base class from within another
        if args:
            if len(self.__dict__):
                raise ConfigError("Something is wrong. You are ineriting `ContigsDBWorkflow` from \
                                   within another class, yet you are providing an `args` parameter.\
                                   This is not alright.")
            self.args = args
            self.name = 'contigs'
        else:
            if not len(self.__dict__):
                raise ConfigError("When you are *not* inheriting `ContigsDBWorkflow` from within\
                                   a super class, you must provide an `args` parameter.")

            if 'name' not in self.__dict__:
                raise ConfigError("The super class trying to inherit `ContigsDBWorkflow` does not\
                                   have a set `self.name`. Which means there may be other things\
                                   wrong with it, hence anvi'o refuses to continue.")

        self.run = run
        self.progress = progress

        # initialize the base class
        WorkflowSuperClass.__init__(self)

        self.rules.extend(['anvi_script_reformat_fasta', 'remove_human_dna_using_centrifuge',
                           'anvi_gen_contigs_database', 'export_gene_calls_for_centrifuge', 'centrifuge',
                           'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database', 'anvi_get_sequences_for_gene_calls',
                           'emapper', 'anvi_script_run_eggnog_mapper'])

        self.general_params.extend(["fasta_txt"])

        self.dirs_dict.update({"FASTA_DIR": "01_FASTA",
                               "CONTIGS_DIR": "02_CONTIGS"})

        self.default_config.update({"fasta_txt": "fasta.txt",
                                    "anvi_gen_contigs_database": {"--project-name": "{group}", "threads": 5},
                                    "centrifuge": {"threads": 5},
                                    "anvi_run_hmms": {"run": True, "threads": 20},
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

        self.rule_acceptable_params_dict['remove_human_dna_using_centrifuge'] = ['run']

        gen_contigs_params = ['--description', '--skip-gene-calling', '--external-gene-calls',\
                              '--ignore-internal-stop-codons', '--skip-mindful-splitting',\
                              '--contigs-fasta', '--project-name',\
                              '--description', '--split-length', '--kmer-size',\
                              '--skip-mindful-splitting', '--skip-gene-calling', '--external-gene-calls',\
                              '--ignore-internal-stop-codons']

        self.rule_acceptable_params_dict['anvi_gen_contigs_database'] = gen_contigs_params
