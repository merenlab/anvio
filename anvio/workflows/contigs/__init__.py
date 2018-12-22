# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import os
import anvio
import anvio.utils as u
import anvio.terminal as terminal
import anvio.workflows as w
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
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
        self.init_workflow_super_class(args, workflow_name='contigs')

        self.group_names = []
        self.contigs_information = {}
        self.fasta_txt_file = None
        self.fasta_information = {}
        # we have external_genomes_file defined here for the sake of pangenomics and phylogenomics workflows
        self.external_genomes_file = ''
        # we have references_mode defined here for the sake of the metagenomics workflow (it is only used when this workflow is inherited)
        self.references_mode = None

        self.rules.extend(['gen_external_genome_file',
                           'anvi_script_reformat_fasta',
                           'anvi_gen_contigs_database', 'export_gene_calls_for_centrifuge', 'centrifuge',
                           'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database', 'anvi_get_sequences_for_gene_calls',
                           'emapper', 'anvi_script_run_eggnog_mapper', 'gunzip_fasta'])
                           'emapper', 'anvi_script_run_eggnog_mapper', 'gunzip_fasta',
                           'translate_external_gene_calls_table'])

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

        self.fasta_txt_file = self.get_param_value_from_config('fasta_txt', repress_default=True)

        if self.fasta_txt_file:
            filesnpaths.is_file_exists(self.fasta_txt_file)
            self.contigs_information = u.get_TAB_delimited_file_as_dictionary(self.fasta_txt_file)
            self.fasta_information.update(self.contigs_information)
            self.group_names = list(self.contigs_information.keys())
            self.references_mode = True
            self.sanity_check_for_fasta_txt()

        self.sanity_check_contigs_project_name()


    def sanity_check_contigs_project_name(self):
        contigs_project_name = self.get_param_value_from_config(['anvi_gen_contigs_database', '--project-name'], repress_default=True)
        if contigs_project_name != self.default_config['anvi_gen_contigs_database']['--project-name'] and contigs_project_name is not None:
            self.run.warning('You chose to set the "project_name" for your contigs databases\
                         in the config file to %s. You are welcomed to do that, but at your own\
                         risk. Just so you know, by default the project name would match\
                         the name for each contigs file (as defined either in the samples_txt\
                         or fasta_txt file that you supplied), by choosing to provide\
                         a different name, it means that all your contigs databases would have\
                         the same name, unless you incloded "{group}" in the name you provided\
                         but even then, we did not test that option and we are not sure it would\
                         work...' % contigs_project_name)


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


    def get_input_for_anvi_gen_contigs_database(self, wildcards):
        d = {}
        d['fasta'] = self.get_fasta(wildcards)
        external_gene_calls = self.contigs_information[wildcards.group].get('external_gene_calls', None)
        if external_gene_calls:
            d['external_gene_calls'] = os.path.join(self.dirs_dict["FASTA_DIR"], wildcards.group + '-external-gene-calls-reformatted.txt')

        return d


    def get_external_gene_calls_param(self, wildcards):
        external_gene_calls = self.contigs_information[wildcards.group].get('external_gene_calls', None)
        if external_gene_calls:
            return "--external-gene-calls " + external_gene_calls
        else:
            return ""


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


    def sanity_check_for_fasta_txt(self):
        """ Run sanity checks on the fasta txt file"""
        columns = next(iter(self.contigs_information.values()))
        bad_columns = [c for c in columns if c not in w.get_fields_for_fasta_information()]
        if bad_columns:
            raise ConfigError("Your fasta_txt file contains columns that are \
                               not familiar to us. These are the only columns \
                               that we accept: '%s'. These are the columns that \
                               we don't like in your file: '%s'." % (", ".join(w.get_fields_for_fasta_information()), \
                                                                   ", ".join(bad_columns)))

        contigs_with_external_functions_and_no_external_gene_calls = \
                [c for c in self.contigs_information \
                    if self.contigs_information[c].get('gene_functional_annotation')
                    and not self.contigs_information[c].get('external_gene_calls')]
        for c in self.contigs_information:
            w.D(self.contigs_information[c].get('external_gene_calls'))
        if contigs_with_external_functions_and_no_external_gene_calls:
            raise ConfigError('You can only provide gene_functional_annotation in \
                               your fasta_txt if you also provide external_gene_calls. \
                               The following entries in "%s" only have functions, but no \
                               gene calls: "%s".' % (self.fasta_txt_file, ', '.join(contigs_with_external_functions_and_no_external_gene_calls)))
