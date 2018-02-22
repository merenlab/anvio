import anvio.workflows as w
from anvio.errors import ConfigError
import anvio.filesnpaths as filesnpaths

class MetagenomicsWorkflow:
    def __init__(self):
        # This dictionary defines the parameters that could be changed for each rule
        self.rules = ['gen_configs', 'qc', 'gen_qc_report', 'gzip_fastqs',\
                     'fq2fa', 'merge_fastas_for_co_assembly', 'megahit', 'anvi_script_anvi_script_reformat_fasta',\
                     'anvi_gen_contigs_database', 'export_gene_calls', 'centrifuge',\
                     'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam',\
                     'anvi_profile', 'annotate_contigs_database', 'anvi_merge']

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule
        rule_acceptable_params_dict['gen_configs'] = []
        rule_acceptable_params_dict['qc'] = []
        rule_acceptable_params_dict['gen_qc_report'] = []
        rule_acceptable_params_dict['gzip_fastqs'] = []
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


class WorkflowSuperClass:
    def __init__(self, config):
            self.config = config
            self.rules = []
            self.rule_acceptable_params_dict = {}
            self.dirs_dict = {}
            self.general_params = []

    def init(self):

        for rule in self.rules:
            if rule not in self.rule_acceptable_params_dict:
                self.rule_acceptable_params_dict[rule] = []

            if 'threads' not in self.rule_acceptable_params_dict[rule]:
                # this should be acceptable for any rule
                self.rule_acceptable_params_dict[rule].append('threads')

        self.dirs_dict = w.get_dir_names(self.config)

        # make sure that config file doesn't have garbage
        self.check_config()


    def check_config(self):
        
        acceptable_params = set(self.rules + self.general_params)
        wrong_params = [p for p in self.config if p not in acceptable_params]
        if wrong_params:
            raise ConfigError("some of the parameters in your config file are not familiar to us. \
                        Here is a list of the wrong parameters: %s. This workflow only accepts \
                        the following general parameters: %s. And these are the rules in this \
                        workflow: %s." % (wrong_params, self.general_params, self.rules))
        
        self.check_rule_params()

    def check_rule_params(self):
        for rule in self.rules:
            if rule in self.config:
                wrong_params = [p for p in self.config[rule] if p not in self.rule_acceptable_params_dict[rule]]
                if wrong_params:
                    raise ConfigError("some of the parameters in your config file for rule %s are not familiar to us. \
                                Here is a list of the wrong parameters: %s. The only acceptable \
                                parameters for this rule are %s." % (rule, wrong_params, self.rule_acceptable_params_dict))


    def save_empty_config_in_json_format(self, filename='empty_config.json'):

        filesnpaths.is_output_file_writable(filename)

        empty_config = self.get_empty_config()

        import json
        with open(filename, 'w') as fp:
            json.dump(empty_config, fp)


    def get_empty_config(self):
        ''' This returns a dictionary with all the possible configurables for a workflow'''

        empty_config = {}

        for rule in self.rules:
            empty_config[rule] = {}
            for param in self.rule_acceptable_params_dict[rule]:
                empty_config[rule][param] = ''

        for param in self.general_params:
            empty_config[param] = ''

        return empty_config


class ContigsDBWorkflow(WorkflowSuperClass):
    def __init__(self, config):
        WorkflowSuperClass.__init__(self, config)

        self.rules.extend(['anvi_script_reformat_fasta', 'remove_human_dna_using_centrifuge',
                           'anvi_gen_contigs_database', 'export_gene_calls', 'centrifuge',
                           'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database'])

        self.general_params.extend(["references_txt"])

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

class PangenomicsWorkflow(ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, config):
        ContigsDBWorkflow.__init__(self, config)

        self.rules.extend(['gen_external_genome_file', 'anvi_gen_genomes_storage',\
                      'anvi_pan_genome'])

        self.general_params.extend(["project_name", "samples_txt"])
        
        pan_params = ["--project-name", "--output-dir", "--genome-names", "--skip-alignments",\
                     "--align-with", "--exclude-partial-gene-calls", "--use-ncbi-blast",\
                     "--minbit", "--mcl-inflation", "--min-occurrence",\
                     "--min-percent-identity", "--sensitive", "--description",\
                     "--overwrite-output-destinations", "--skip-hierarchical-clustering",\
                     "--enforce-hierarchical-clustering", "--distance", "--linkage"]
        self.rule_acceptable_params_dict['anvi_pan_genome'] = pan_params

        storage_params = ["--internal-genomes", "--external-genomes", "--gene-caller"]
        self.rule_acceptable_params_dict['anvi_gen_genomes_storage'] = storage_params
