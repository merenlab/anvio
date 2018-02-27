# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import anvio
import anvio.terminal as terminal

from anvio.workflows import WorkflowSuperClass
from anvio.workflows.contigs import ContigsDBWorkflow


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class MetagenomicsWorkflow(ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self):
        ContigsDBWorkflow.__init__(self, config)

        self.rules = ['iu_gen_configs', 'iu_filter_quality_minoche', 'gen_qc_report', 'gzip_fastqs',\
                     'fq2fa', 'merge_fastas_for_co_assembly', 'megahit', 'anvi_script_anvi_script_reformat_fasta',\
                     'anvi_gen_contigs_database', 'export_gene_calls', 'centrifuge',\
                     'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam',\
                     'anvi_profile', 'annotate_contigs_database', 'anvi_merge']

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule
        rule_acceptable_params_dict['iu_gen_configs'] = ["--r1-prefix", "--r2-prefix"]
        rule_acceptable_params_dict['iu_filter_quality_minoche'] = ['visualize_quality_curves', 'ignore_deflines', 'limit_num_pairs', 'print_qual_scores', 'store_read_fate']
        rule_acceptable_params_dict['gzip_fastqs'] = ["run"]
        rule_acceptable_params_dict['fq2fa'] = []
        rule_acceptable_params_dict['merge_fastas_for_co_assembly'] = []
        rule_acceptable_params_dict['megahit'] = []
        rule_acceptable_params_dict['anvi_script_reformat_fasta'] = []
        rule_acceptable_params_dict['anvi_gen_contigs_database'] = []
        rule_acceptable_params_dict['export_gene_calls'] = []
        rule_acceptable_params_dict['centrifuge'] = []
        rule_acceptable_params_dict['anvi_import_taxonomy'] = []
        rule_acceptable_params_dict['anvi_run_hmms'] = []
        rule_acceptable_params_dict['anvi_run_ncbi_cogs'] = []
        rule_acceptable_params_dict['bowtie_build'] = []
        rule_acceptable_params_dict['bowtie'] = []
        rule_acceptable_params_dict['samtools_view'] = []
        rule_acceptable_params_dict['anvi_init_bam'] = []
        rule_acceptable_params_dict['anvi_profile'] = []
        rule_acceptable_params_dict['annotate_contigs_database'] = []
        rule_acceptable_params_dict['anvi_merge'] = []

        self.rule_acceptable_params_dict = rule_acceptable_params_dict
