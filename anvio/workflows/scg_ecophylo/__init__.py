# -*- coding: utf-8
# pylint: disable=line-too-long


import os
import anvio
import pandas as pd
import anvio.utils as u
import anvio.terminal as terminal
import anvio.workflows as w
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.workflows.metagenomics import MetagenomicsWorkflow


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


run = terminal.Run()

class SCGEcoPhyloWorkflow(WorkflowSuperClass):

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='scg_ecophylo')

        # Snakemake rules
        self.rules.extend(['anvi_run_hmms_hmmsearch',
                           'filter_hmm_hits_by_query_coverage',
                           'anvi_estimate_scg_taxonomy_for_SCGs',
                           'add_misc_data_to_taxonomy',
                           'anvi_get_external_gene_calls_file',
                           'rename_and_filter_external_gene_calls_file',
                           'cat_external_gene_calls_file',
                           'subset_external_gene_calls_file_all',
                           'cluster_X_percent_sim_mmseqs',
                           'subset_AA_seqs_with_mmseqs_reps',
                           'align_sequences',
                           'trim_alignment',
                           'remove_sequences_with_X_percent_gaps',
                           'count_num_sequences_filtered',
                           'subset_DNA_reps_with_QCd_AA_reps_for_mapping',
                           'make_fasta_txt',
                           'make_metagenomics_config_file',
                           'run_metagenomics_workflow',
                           'make_anvio_state_file',
                           'add_default_collection',
                           'anvi_summarize',
                           'rename_tree_tips',
                           'anvi_import_state',
                           'fasttree',
                           'iqtree'
                           ])

        self.general_params.extend(['metagenomes']) # user needs to input a metagenomes.txt file
        self.general_params.extend(['external_genomes']) # user can add isolate genomes if needed
        self.general_params.extend(['SCG_protein_list']) # user must input which Reference proteins will be used for workflow

        # Parameters for each rule that are accessible in the config.json file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['anvi_run_hmms_hmmsearch'] = ['--installed-hmm-profile']
        rule_acceptable_params_dict['filter_hmm_hits_by_query_coverage'] = ['--hmm-source', '--query-coverage', 'additional_params']
        rule_acceptable_params_dict['anvi_estimate_scg_taxonomy_for_SCGs'] = ['--metagenome-mode']
        rule_acceptable_params_dict['cluster_X_percent_sim_mmseqs'] = ['--min-seq-id']
        rule_acceptable_params_dict['trim_alignment'] = ['-gt', "-gappyout", 'additional_params']
        rule_acceptable_params_dict['remove_sequences_with_X_percent_gaps'] = ['--max-percentage-gaps']
        rule_acceptable_params_dict['fasttree'] = ['run']
        rule_acceptable_params_dict['iqtree'] = ['run', '-m', 'additional_params']
        rule_acceptable_params_dict['run_metagenomics_workflow'] = ['clusterize', 'cluster_submission_params']

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Set default values for parameters accessible in config.json file
        self.default_config.update({
            'metagenomes': 'metagenomes.txt',
            'external_genomes': 'external-genomes.txt',
            'SCG_protein_list': 'SCG_protein_list.txt',
            'samples_txt': "samples.txt",
            'anvi_run_hmms_hmmsearch': {'threads': 5, '--installed-hmm-profile': 'Bacteria_71'},
            'filter_hmm_hits_by_query_coverage': {'threads': 5, '--query-coverage': 0.8, '--hmm-source': 'Bacteria_71'},
            'anvi_estimate_scg_taxonomy_for_SCGs': {'threads': 5},
            'cluster_X_percent_sim_mmseqs': {'threads': 5, '--min-seq-id': 0.94},
            'align_sequences': {'threads': 5},
            'remove_sequences_with_X_percent_gaps': {'threads': 5, '--max-percentage-gaps': 50},
            'trim_alignment': {'threads': 5, '-gappyout': True},
            'fasttree': {'run': True, 'threads': 5},
            'iqtree': {'threads': 5,'-m': "MFP"},
            'add_misc_data_to_taxonomy': {'threads': 5},
            'subset_AA_seqs_with_mmseqs_reps': {'threads': 5},
            'anvi_get_external_gene_calls_file': {'threads': 5},
            'subset_DNA_reps_with_QCd_AA_reps_for_mapping': {'threads': 5},
            'count_num_sequences_filtered': {'threads': 5},
            'rename_and_filter_external_gene_calls_file': {'threads': 2},
            'cat_external_gene_calls_file': {'threads': 2},
            'subset_external_gene_calls_file_all': {'threads': 5},
            'trim_alignment': {'threads': 5},
            'make_fasta_txt': {'threads': 2},
            'make_metagenomics_config_file': {'threads': 5},
            'make_anvio_state_file': {'threads': 2},
            'add_default_collection': {'threads': 2},
            'anvi_summarize': {'threads': 10},
            'rename_tree_tips': {'threads': 2},
            'anvi_import_state': {'threads': 2},
            'run_metagenomics_workflow': {'threads': 10, 'clusterize': False}
            })

        # Directory structure for Snakemake workflow
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_DIR": "SCG_ECO_PHYLO_WORKFLOW/01_REFERENCE_PROTEIN_DATA"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS": "SCG_ECO_PHYLO_WORKFLOW/02_NR_FASTAS"})
        self.dirs_dict.update({"MSA": "SCG_ECO_PHYLO_WORKFLOW/03_MSA"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_STATS": "SCG_ECO_PHYLO_WORKFLOW/04_SEQUENCE_STATS"})
        self.dirs_dict.update({"TREES": "SCG_ECO_PHYLO_WORKFLOW/05_TREES"})
        self.dirs_dict.update({"MISC_DATA": "SCG_ECO_PHYLO_WORKFLOW/06_MISC_DATA"})
        self.dirs_dict.update({"SCG_NT_FASTAS": "SCG_ECO_PHYLO_WORKFLOW/07_SCG_NT_FASTAS"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS_RENAMED": "SCG_ECO_PHYLO_WORKFLOW/08_RIBOSOMAL_PROTEIN_FASTAS_RENAMED"})


    def init(self):
        """This function is called from within the snakefile to initialize parameters."""

        super().init()

        # Re-assigning LOGS/ dir to inside SCG_ECO_PHYLO_WORKFLOW/ dir
        self.dirs_dict.update({"LOGS_DIR": "SCG_ECO_PHYLO_WORKFLOW/00_LOGS"})

        self.names_list = []
        self.names_dirs = []

        # Load metagenomes.txt
        self.metagenomes = self.get_param_value_from_config(['metagenomes'])
        self.input_dirs_dict = {}

        if self.metagenomes:
            filesnpaths.is_file_exists(self.metagenomes)
            try:
                self.metagenomes_df = pd.read_csv(self.metagenomes, sep='\t', index_col=False)
                self.metagenomes_name_list = self.metagenomes_df.name.to_list()
                self.metagenomes_path_list = self.metagenomes_df.contigs_db_path.to_list()
                self.metagenomes_dirname_list = [os.path.dirname(x) for x in self.metagenomes_path_list]
                self.input_dirs_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_dirname_list)))
                self.names_dirs.extend(self.metagenomes_dirname_list)
                self.names_list.extend(self.metagenomes_name_list)

            except IndexError as e:
                raise ConfigError("The metagenomes.txt file, '%s', does not appear to be properly formatted. "
                                  "This is the error from trying to load it: '%s'" % (self.metagenomes_df, e))

        # Load external-genomes.txt
        self.external_genomes = self.get_param_value_from_config(['external_genomes'])
        if self.external_genomes:
            filesnpaths.is_file_exists(self.external_genomes)
            try:
                self.external_genomes_df = pd.read_csv(self.external_genomes, sep='\t', index_col=False)
                self.external_genomes_names_list = self.external_genomes_df.name.to_list()
                self.external_genomes_path_list = self.external_genomes_df.contigs_db_path.to_list()
                self.external_genomes_dirname_list = [os.path.dirname(x) for x in self.external_genomes_path_list]
                self.input_dirs_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_dirname_list)))
                self.names_dirs.extend(self.external_genomes_dirname_list)
                self.names_list.extend(self.external_genomes_names_list)

            except IndexError as e:
                raise ConfigError("The external-genomes.txt file, '%s', does not appear to be properly formatted. "
                                  "This is the error from trying to load it: '%s'" % (self.external_genomes_df, e))
        else:
            self.external_genomes_names_list = []

        # Make a unique list
        self.names_dirs = list(set(self.names_dirs))

        # Make variables that tells whether we have metagenomes.txt, external-genomes.txt, or both
        if self.metagenomes and not self.external_genomes:
            self.mode = 'metagenomes'
        if not self.metagenomes and self.external_genomes:
            self.mode = 'external_genomes'
        if self.metagenomes and self.external_genomes:
            self.mode = 'both'

        # Load Reference protein list
        self.SCG_protein_list_path = self.get_param_value_from_config(['SCG_protein_list'])
        filesnpaths.is_file_exists(self.SCG_protein_list_path)
        try:
            self.Reference_protein_df = pd.read_csv(self.SCG_protein_list_path, sep='\t', index_col=False)
            self.SCG_protein_list = self.Reference_protein_df.iloc[:, 0].to_list() 

        except IndexError as e:
            raise ConfigError("The SCG_protein_list.txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.Ribosomal_protein_df, e))

        # Load samples.txt
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])

        if not self.samples_txt_file:
            if os.path.exists('samples.txt'):
                self.samples_txt_file = 'samples.txt'
            else:
                raise ConfigError("Ehem. Your config file does not include a `samples_txt` directive. "
                                  "Anvi'o tried to assume that your `samples.txt` may be in your work "
                                  "directory, but you don't seem to have a `samples.txt` file anywhere "
                                  "around either. So please add a `samples.txt` directive.")

        filesnpaths.is_file_tab_delimited(self.samples_txt_file)
        try:
            # getting the samples information (names, [group], path to r1, path to r2) from samples.txt
            self.samples_information = pd.read_csv(self.samples_txt_file, sep='\t', index_col=False)
        except IndexError as e:
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. "
                              "This is what we know: '%s'" % (self.samples_txt_file, e))
        if 'sample' not in list(self.samples_information.columns):
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. "
                              "We are not sure what's wrong, but we can't find a column with title 'sample'." % self.samples_txt_file)

        self.sample_names_for_mapping_list = self.samples_information['sample'].to_list()

        # Pick which tree algorithm
        self.run_iqtree = self.get_param_value_from_config(['iqtree', 'run'])
        self.run_fasttree = self.get_param_value_from_config(['fasttree', 'run'])

        if not self.run_iqtree and not self.run_fasttree:
            raise ConfigError("Please choose either iqtree or fasttree in your config file to run your phylogenetic tree.")

        # Decide to clusterize metagenomic workflow
        self.clusterize_metagenomics_workflow = self.get_param_value_from_config(['run_metagenomics_workflow', 'clusterize'])
        self.metagenomics_workflow_HPC_string = self.get_param_value_from_config(['run_metagenomics_workflow', 'cluster_submission_params'])

        self.target_files = self.get_target_files()

    def get_target_files(self):
        target_files = []

        for reference_protein_name in self.SCG_protein_list:

            # Count num sequences removed per step
            tail_path = "%s_stats.tsv" % (reference_protein_name)
            target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], reference_protein_name, tail_path)
            target_files.append(target_file)

            # Import state file to customize interactive interface
            target_file = os.path.join(f"{reference_protein_name}_state_imported.done")
            target_files.append(target_file)

            # Reformat names to match splits in interactive interface
            target_file = os.path.join(f"{reference_protein_name}_combined.done")
            target_files.append(target_file)

            target_file = os.path.join(self.dirs_dict['MISC_DATA'], f"{reference_protein_name}_misc.tsv")
            target_files.append(target_file)

            # The FINAL trees :)
            # iq-tree
            if self.run_iqtree == True:
                tail_path = "%s.iqtree" % (reference_protein_name)
                target_file = os.path.join(self.dirs_dict['TREES'], reference_protein_name, tail_path)
                target_files.append(target_file)
            # fasttree
            elif self.run_fasttree == True:
                tail_path = "%s.nwk" % (reference_protein_name)
                target_file = os.path.join(self.dirs_dict['TREES'], reference_protein_name, tail_path)
                target_files.append(target_file)

        return target_files
