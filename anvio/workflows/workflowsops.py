import anvio.workflows as w
from anvio.errors import ConfigError

class MetagenomicsWorkflow:
    def __init__(self):
        # This dictionary defines the parameters that could be changed for each rule
        self.rules = ['gen_configs', 'qc', 'gen_qc_report', 'gzip_fastqs',\
                     'fq2fa', 'merge_fastas_for_co_assembly', 'megahit', 'reformat_fasta',\
                     'gen_contigs_db', 'export_gene_calls', 'run_centrifuge',\
                     'import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam',\
                     'anvi_profile', 'annotate_contigs_database', 'anvi_merge']

        acceptable_params_dict = {}

        # defining the accesible params per rule
        acceptable_params_dict['gen_configs'] = []
        acceptable_params_dict['qc'] = []
        acceptable_params_dict['gen_qc_report'] = []
        acceptable_params_dict['gzip_fastqs'] = []
        acceptable_params_dict['fq2fa'] = []
        acceptable_params_dict['merge_fastas_for_co_assembly'] = []
        acceptable_params_dict['megahit'] = []
        acceptable_params_dict['reformat_fasta'] = []
        acceptable_params_dict['gen_contigs_db'] = []
        acceptable_params_dict['export_gene_calls'] = []
        acceptable_params_dict['run_centrifuge'] = []
        acceptable_params_dict['import_taxonomy'] = []
        acceptable_params_dict['anvi_run_hmms'] = []
        acceptable_params_dict['anvi_run_ncbi_cogs'] = []
        acceptable_params_dict['bowtie_build'] = []
        acceptable_params_dict['bowtie'] = []
        acceptable_params_dict['samtools_view'] = []
        acceptable_params_dict['anvi_init_bam'] = []
        acceptable_params_dict['anvi_profile'] = []
        acceptable_params_dict['annotate_contigs_database'] = []
        acceptable_params_dict['anvi_merge'] = []

        self.acceptable_params_dict = acceptable_params_dict


class WorkflowSuperClass:
    def __init__(self, config):
            self.config = config
            self.rules = ['all']
            self.acceptable_params_dict = {}
            self.dirs_dict = {}

    def init(self):

        for rule in self.rules:
            if rule == 'all':
                continue
            elif rule not in self.acceptable_params_dict:
                self.acceptable_params_dict[rule] = []

            if 'threads' not in self.acceptable_params_dict[rule]:
                # this should be acceptable for any rule except rule 'all'
                self.acceptable_params_dict[rule].append('threads')

        self.dirs_dict = w.get_dir_names(self.config)

        self.check_rules_params()


    def check_rules_params(self):
        for rule in self.rules:
            if rule in self.config:
                wrong_params = [p for p in self.config[rule] if p not in self.acceptable_params_dict[rule]]
                if wrong_params:
                    ConfigError("some of the parameters in your config file for rule %s are wrong. \
                                Here is a list of the wrong parameters: %s. The only acceptable \
                                parameters for this rule are %s." % (rule, wrong_params, acceptable_params_dict))


class ContigsDBWorkflow(WorkflowSuperClass):
    def __init__(self, config):
        WorkflowSuperClass.__init__(self, config)
        self.rules.extend(['reformat_fasta', 'remove_human_dna_using_centrifuge',
                           'gen_contigs_db', 'export_gene_calls', 'run_centrifuge',
                           'import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',
                           'annotate_contigs_database'])


class PangenomicsWorkflow(ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, config):
        ContigsDBWorkflow.__init__(self, config)

        self.rules.extend(['gen_external_genome_file', 'anvi_gen_genomes_storage',\
                      'anvi_pan_genome'])
        
        pan_params = ["--project-name", "--output-dir", "--genome-names", "--skip-alignments",\
                     "--align-with", "--exclude-partial-gene-calls", "--use-ncbi-blast",\
                     "--minbit", "--mcl-inflation", "--min-occurrence",\
                     "--min-percent-identity", "--sensitive", "--description",\
                     "--overwrite-output-destinations", "--skip-hierarchical-clustering",\
                     "--enforce-hierarchical-clustering", "--distance", "--linkage"]
        self.acceptable_params_dict['anvi_pan_genome'] = pan_params

        storage_params = ["--internal-genomes", "--external-genomes", "--gene-caller"]
        self.acceptable_params_dict['anvi_gen_genomes_storage'] = storage_params

        self.init()
