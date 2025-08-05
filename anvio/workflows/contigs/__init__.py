# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import os
import anvio
import anvio.terminal as terminal
import anvio.workflows as w
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.utils.files import get_TAB_delimited_file_as_dictionary
from anvio.utils.validation import is_this_name_OK_for_database


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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
        self.import_external_functions_flags =[]

        self.rules.extend(['gen_external_genome_file',
                           'anvi_script_reformat_fasta',
                           'anvi_gen_contigs_database', 'export_gene_calls_for_centrifuge', 'centrifuge',
                           'anvi_import_taxonomy_for_genes', 'anvi_run_scg_taxonomy', 'anvi_run_trna_scan', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database', 'anvi_get_sequences_for_gene_calls', 'emapper',
                           'anvi_script_run_eggnog_mapper', 'gunzip_fasta', 'reformat_external_gene_calls_table',
                           'reformat_external_functions', 'import_external_functions', 'anvi_run_pfams', 'anvi_run_kegg_kofams'])

        self.general_params.extend(["fasta_txt"])

        self.dirs_dict.update({"FASTA_DIR": "01_FASTA",
                               "CONTIGS_DIR": "02_CONTIGS"})

        self.default_config.update({"fasta_txt": "fasta.txt",
                                    "anvi_gen_contigs_database": {"--project-name": "{group}"},
                                    "centrifuge": {"threads": 2},
                                    "anvi_run_hmms": {"run": True, "threads": 5, "--also-scan-trnas": True},
                                    "anvi_run_kegg_kofams": {"run": True, "threads": 4},
                                    "anvi_run_ncbi_cogs": {"run": True, "threads": 5},
                                    "anvi_run_scg_taxonomy": {"run": True, "threads": 6},
                                    "anvi_run_trna_scan": {"run": False, "threads": 6},
                                    "anvi_script_reformat_fasta": {"run": True, "--prefix": "{group}", "--simplify-names": True},
                                    "emapper": {"--database": "bact", "--usemem": True, "--override": True},
                                    "anvi_script_run_eggnog_mapper": {"--use-version": "0.12.6"}})

        self.rule_acceptable_params_dict['anvi_run_ncbi_cogs'] = ['run', '--cog-data-dir', '--temporary-dir-path', '--search-with']

        self.rule_acceptable_params_dict['anvi_run_scg_taxonomy'] = ['run', '--scgs-taxonomy-data-dir']

        self.rule_acceptable_params_dict['anvi_run_trna_scan'] = ['run', '--trna-cutoff-score']

        self.rule_acceptable_params_dict['anvi_run_hmms'] = ['run', '--installed-hmm-profile', '--hmm-profile-dir', '--also-scan-trnas', '--add-to-functions-table']

        self.rule_acceptable_params_dict['anvi_run_pfams'] = ['run', '--pfam-data-dir']

        self.rule_acceptable_params_dict['anvi_run_kegg_kofams'] = ['run', '--kegg-data-dir', '--hmmer-program', '--keep-all-hits', '--log-bitscores', '--just-do-it']

        self.rule_acceptable_params_dict['centrifuge'] = ['run', 'db']

        self.rule_acceptable_params_dict['emapper'] = ['--database', '--usemem', '--override', 'path_to_emapper_dir']

        self.rule_acceptable_params_dict['anvi_script_run_eggnog_mapper'] = ['run', '--cog-data-dir', '--drop-previous-annotations',
                                         '--use-version']

        self.rule_acceptable_params_dict['anvi_script_reformat_fasta'] = \
                    ['run', '--keep-ids', '--exclude-ids', '--min-len', "--prefix", "--simplify-names", "--seq-type"]


        gen_contigs_params = ['--description', '--skip-gene-calling',\
                              '--ignore-internal-stop-codons', '--skip-mindful-splitting',\
                              '--contigs-fasta', '--project-name',\
                              '--description', '--split-length', '--kmer-size',\
                              '--skip-mindful-splitting', '--skip-gene-calling',\
                              '--ignore-internal-stop-codons', '--skip-predict-frame', '--prodigal-translation-table']

        self.rule_acceptable_params_dict['anvi_gen_contigs_database'] = gen_contigs_params


    def init(self):
        super().init()

        self.fasta_txt_file = self.get_param_value_from_config('fasta_txt')

        if self.fasta_txt_file:
            if not filesnpaths.is_file_exists(self.fasta_txt_file, dont_raise=True):
                raise ConfigError('You know the path you have for `fasta_txt` in your config file? There is no such file on your disk :(')

            self.contigs_information = get_TAB_delimited_file_as_dictionary(self.fasta_txt_file)
            self.fasta_information.update(self.contigs_information)
            self.group_names = list(self.contigs_information.keys())
            self.references_mode = True
            self.sanity_check_for_fasta_txt()

        self.sanity_check_contigs_project_name()

        # check and warn user regarding risky change of parameters with wildcards as default values
        self.warn_user_regarding_param_with_wildcard_default_value('anvi_run_ncbi_cogs', '--temporary-dir-path', '{group}')
        self.warn_user_regarding_param_with_wildcard_default_value('anvi_script_reformat_fasta', '--prefix', '{group}')


    def get_contigs_target_files(self):
        # Add the contigs database as a target file.
        # This is the only mandatory target file for this workflow.
        target_files = []
        target_files.append(self.get_contigs_db_path())

        # Add the optional target files according to user configurations
        target_files.extend(self.get_contigs_workflow_optional_targets())

        return target_files


    def get_contigs_workflow_optional_targets(self):
        optional_targets = []

        run_taxonomy_with_centrifuge = self.get_param_value_from_config(["centrifuge", "run"]) == True
        # sanity check for centrifuge db
        if run_taxonomy_with_centrifuge:
            if not self.get_param_value_from_config(["centrifuge", "db"]):
                raise ConfigError("If you plan to run centrifuge, then you must "\
                                  "provide a path for the centrifuge db in the "\
                                  "config file. See documentation for more details.")

        if run_taxonomy_with_centrifuge:
            optional_targets.append(self.dirs_dict["CONTIGS_DIR"] + "/{group}-anvi_anvi_import_taxonomy_for_genes.done")

        run_anvi_run_hmms = self.get_param_value_from_config(["anvi_run_hmms", "run"]) == True
        if run_anvi_run_hmms:
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "anvi_run_hmms-{group}.done"))

        if self.get_param_value_from_config(["anvi_run_ncbi_cogs", "run"]) == True:
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "anvi_run_ncbi_cogs-{group}.done"))

        run_anvi_run_kofams = self.get_param_value_from_config(["anvi_run_kegg_kofams", "run"]) == True
        if run_anvi_run_kofams:
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "anvi_run_kegg_kofams-{group}.done"))

        run_anvi_run_scg_taxonomy = self.get_param_value_from_config(["anvi_run_scg_taxonomy", "run"]) == True
        if run_anvi_run_scg_taxonomy:
            if not run_anvi_run_hmms:
                self.run.warning('You chose to run anvi_run_scg_taxonomy, but you didn\'t choose to run '
                                 'anvi_run_hmms. Continue at your own risk. If your contigs databases '
                                 'don\'t have HMM hits stored already then anvi_run_scg_taxonomy will fail.')
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "anvi_run_scg_taxonomy-{group}.done"))

        run_anvi_run_trna_scan = self.get_param_value_from_config(["anvi_run_trna_scan", "run"]) == True
        if run_anvi_run_trna_scan:
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "anvi_run_trna_scan-{group}.done"))

        if self.get_param_value_from_config(["anvi_script_run_eggnog_mapper", "run"]) == True:
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "{group}-anvi_script_run_eggnog_mapper.done"))

        # import external functions if provided
        import_external_functions_flags = [os.path.join(self.dirs_dict["CONTIGS_DIR"],
                                           group + "-external-functions-imported.done")\
                                           for group in self.contigs_information \
                                           if self.contigs_information[group].get('gene_functional_annotation')]
        if import_external_functions_flags:
            optional_targets.extend(import_external_functions_flags)

        if self.get_param_value_from_config(['anvi_run_pfams', 'run']):
            optional_targets.append(os.path.join(self.dirs_dict["CONTIGS_DIR"], "anvi_run_pfams-{group}.done"))

        return optional_targets


    def sanity_check_contigs_project_name(self):
        contigs_project_name = self.get_param_value_from_config(['anvi_gen_contigs_database', '--project-name'])
        if contigs_project_name != self.default_config['anvi_gen_contigs_database']['--project-name'] and contigs_project_name is not None:
            self.run.warning('You chose to set the "project_name" for your contigs databases '
                        'in the config file to %s. You are welcomed to do that, but at your own '
                        'risk. Just so you know, by default the project name would match '
                        'the name for each contigs file (as defined either in the samples_txt '
                        'or fasta_txt file that you supplied), by choosing to provide '
                        'a different name, it means that all your contigs databases would have '
                        'the same name, unless you incloded "{group}" in the name you provided '
                        'but even then, we did not test that option and we are not sure it would '
                        'work...' % contigs_project_name)


    def get_contigs_db_path(self):
        return os.path.join(self.dirs_dict["CONTIGS_DIR"], "{group}-contigs.db")


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
        external_gene_calls_file = self.get_external_gene_calls_file_name(wildcards)
        if external_gene_calls_file:
            d['external_gene_calls'] = external_gene_calls_file
        return d


    def get_input_for_reformat_external_gene_calls_table(self, wildcards):
        d = {}
        d['contigs'] = os.path.join(self.dirs_dict["FASTA_DIR"], wildcards.group, wildcards.group + "-contigs.fa"),
        d['reformat_report'] = os.path.join(self.dirs_dict["FASTA_DIR"], wildcards.group, wildcards.group + "-reformat-report.txt"),
        d['external_gene_calls'] = self.contigs_information[wildcards.group]['external_gene_calls']
        return d


    def get_input_for_reformat_external_functions(self, wildcards):
        d = {}
        d['gene_functional_annotation'] = self.contigs_information[wildcards.group]['gene_functional_annotation']
        d['external_gene_calls'] = self.get_external_gene_calls_file_name(wildcards)
        return d


    def get_input_for_import_external_functions(self, wildcards):
        d = {}
        d['gene_functional_annotation'] = self.get_external_gene_functions_file_name(wildcards)
        d['contigs'] = os.path.join(self.dirs_dict["CONTIGS_DIR"], wildcards.group + "-contigs.db")
        return d


    def get_external_gene_calls_file_name(self, wildcards):
        '''If the user is running anvi-script-reformat-fasta then we need to also reformat the external gene calls'''
        external_gene_calls = self.contigs_information.get(wildcards.group, {}).get('external_gene_calls', None)
        if external_gene_calls:
            reformat = self.get_rule_param('anvi_script_reformat_fasta', 'run')
            if reformat:
                return os.path.join(self.dirs_dict['FASTA_DIR'], wildcards.group, wildcards.group + "-external-gene-calls.txt")
            else:
                return self.contigs_information[wildcards.group]['external_gene_calls']
        return ''


    def get_external_gene_functions_file_name(self, wildcards):
        '''If the user is running anvi-script-reformat-fasta then we need to also reformat the external gene calls'''
        gene_functional_annotation = self.contigs_information[wildcards.group].get('gene_functional_annotation', None)
        if gene_functional_annotation:
            reformat = self.get_rule_param('anvi_script_reformat_fasta', 'run')
            if reformat:
                return os.path.join(self.dirs_dict['FASTA_DIR'], wildcards.group, wildcards.group + "-gene-functional-annotation.txt")
            else:
                return self.contigs_information[wildcards.group]['gene_functional_annotation']
        return ''


    def get_external_gene_calls_param(self, wildcards):
        external_gene_calls = self.get_external_gene_calls_file_name(wildcards)
        if external_gene_calls:
            return "--external-gene-calls " + external_gene_calls
        else:
            return ""


    def get_prodigal_translation_table_param(self, wildcards):
        external_gene_calls_param = self.get_external_gene_calls_param(wildcards)
        if external_gene_calls_param:
            # If the user supplied external gene calls then we ignore the --prodigal-translation-table parameter
            return ''
        else:
            return self.get_rule_param("anvi_gen_contigs_database", "--prodigal-translation-table")


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

        if not self.contigs_information:
            raise ConfigError('It looks the fasta_txt file you provided "%s" is empty.' % self.fasta_txt_file)

        for name in self.contigs_information.keys():
            is_this_name_OK_for_database('fasta.txt entry name', name, additional_chars_allowed='.')

        columns = next(iter(self.contigs_information.values()))
        bad_columns = [c for c in columns if c not in w.get_fields_for_fasta_information()]
        if bad_columns:
            raise ConfigError("Your fasta_txt file contains columns that are "
                              "not familiar to us. These are the only columns "
                              "that we accept: '%s'. These are the columns that "
                              "we don't like in your file: '%s'." % (", ".join(w.get_fields_for_fasta_information()), \
                                                                   ", ".join(bad_columns)))

        contigs_with_external_functions_and_no_external_gene_calls = \
                [c for c in self.contigs_information \
                    if self.contigs_information[c].get('gene_functional_annotation')
                    and not self.contigs_information[c].get('external_gene_calls')]
        if contigs_with_external_functions_and_no_external_gene_calls:
            raise ConfigError('You can only provide gene_functional_annotation in '
                              'your fasta_txt if you also provide external_gene_calls. '
                              'The following entries in "%s" only have functions, but no '
                              'gene calls: "%s".' % (self.fasta_txt_file, ', '.join(contigs_with_external_functions_and_no_external_gene_calls)))
