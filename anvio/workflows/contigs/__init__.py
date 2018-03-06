# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import anvio
import anvio.terminal as terminal

from anvio.workflows import WorkflowSuperClass


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()

class ContigsDBWorkflow(WorkflowSuperClass):
    def __init__(self, config):
        WorkflowSuperClass.__init__(self, config)

        self.rules.extend(['anvi_script_reformat_fasta', 'remove_human_dna_using_centrifuge',
                           'anvi_gen_contigs_database', 'export_gene_calls', 'centrifuge',
                           'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database'])

        self.general_params.extend(["fasta_txt"])

        self.dirs_dict.update(
                {
                    "FASTA_DIR": "01_FASTA",
                    "CONTIGS_DIR": "02_CONTIGS"
                }
                             )

        self.default_config.update(
                {
                    "fasta_txt": "fasta.txt",
                    "anvi_gen_contigs_database": {"--project-name": "{group}", "threads": 5},
                    "centrifuge": {"threads": 5},
                    "anvi_run_hmms": {"run": True, "threads": 20},
                    "anvi_run_ncbi_cogs": {"run": True, "threads": 5},
                    "anvi_script_reformat_fasta": {"run": True, "--simplify-names": True},
                }
                                  )

        self.rule_acceptable_params_dict['anvi_run_ncbi_cogs'] = ['run', '--cogs-data-dir', '--sensitive', '--temporary-dir-path', '--search-with']

        self.rule_acceptable_params_dict['anvi_run_hmms'] = ['run', '--installed-hmm-profile', '--hmm-profile-dir']

        self.rule_acceptable_params_dict['centrifuge'] = ['run']

        self.rule_acceptable_params_dict['anvi_script_reformat_fasta'] = \
                    ['run', '--simplify-names', '--keep-ids', '--exclude-ids', '--min-len']

        self.rule_acceptable_params_dict['remove_human_dna_using_centrifuge'] = ['run']

        gen_contigs_params = ['--description', '--skip-gene-calling', '--external-gene-calls',\
                              '--ignore-internal-stop-codons', '--skip-mindful-splitting',\
                              '--contigs-fasta', '--project-name', '--output-db-path',\
                              '--description', '--split-length', '--kmer-size',\
                              '--skip-mindful-splitting', '--skip-gene-calling', '--external-gene-calls',\
                              '--ignore-internal-stop-codons']

        self.rule_acceptable_params_dict['anvi_gen_contigs_database'] = gen_contigs_params

        self.rules_dependencies.update({'anvi_script_reformat_fasta': "anvi-script-reformat-fasta",
                                        'remove_human_dna_using_centrifuge': "centrifuge",
                                        'anvi_gen_contigs_database': "anvi-gen-contigs-database",
                                        'export_gene_calls': "anvi-export-gene-calls",
                                        'centrifuge': "centrifuge",
                                        'anvi_import_taxonomy': "anvi-import-taxonomy",
                                        'anvi_run_hmms': "anvi-run-hmms",
                                        'anvi_run_ncbi_cogs': "anvi-run-ncbi-cogs",
                                        'annotate_contigs_database': ""})
