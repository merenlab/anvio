# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes to define and work with anvi'o ecophylo workflows. """

from distutils.command.config import config
import os
import anvio
import argparse
import pandas as pd

import anvio
import anvio.data.hmm
import anvio.utils as u
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.genomedescriptions import GenomeDescriptions
from anvio.genomedescriptions import MetagenomeDescriptions


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = ['mschecht']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


run = terminal.Run()

class EcoPhyloWorkflow(WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='ecophylo')

        # Snakemake rules
        self.rules.extend(['anvi_run_hmms_hmmsearch',
                           'filter_hmm_hits_by_model_coverage',
                           'process_hmm_hits',
                           'combine_sequence_data',
                           'anvi_get_external_gene_calls_file',
                           'cat_external_gene_calls_file',
                           'cluster_X_percent_sim_mmseqs',
                           'subset_AA_seqs_with_mmseqs_reps',
                           'align_sequences',
                           'trim_alignment',
                           'remove_sequences_with_X_percent_gaps',
                           'count_num_sequences_filtered',
                           'subset_DNA_reps_with_QCd_AA_reps_for_mapping',
                           'subset_external_gene_calls_file_all',
                           'make_fasta_txt',
                           'fasttree',
                           'iqtree',
                           'make_metagenomics_config_file',
                           'run_metagenomics_workflow',
                           'add_default_collection',
                           'anvi_summarize',
                           'rename_tree_tips',
                           'make_misc_data',
                           'anvi_scg_taxonomy',
                           'make_anvio_state_file',
                           'anvi_import_everything'
                           ])

        self.general_params.extend(['metagenomes']) # user needs to input a metagenomes.txt file
        self.general_params.extend(['external_genomes']) # user can add isolate genomes if needed
        self.general_params.extend(['hmm_list']) # user must input which Reference proteins will be used for workflow
        self.general_params.extend(['samples_txt']) # user must input which Reference proteins will be used for workflow
        self.general_params.extend(['cluster_representative_method']) # pick cluster rep based on single profile coverage values
        self.general_params.extend(['gene_caller_to_use']) # designate gene-caller for all contig-dbs if not default Prodigal
        self.general_params.extend(['run_genomes_sanity_check']) # run GenomeDescriptions and MetagenomesDescriptions

        # Parameters for each rule that are accessible in the config.json file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['anvi_run_hmms_hmmsearch'] = ['threads_genomes', 'threads_metagenomes', 'additional_params']
        rule_acceptable_params_dict['filter_hmm_hits_by_model_coverage'] = ['--model-coverage', '--filter-out-partial-gene-calls', 'additional_params']
        rule_acceptable_params_dict['cluster_X_percent_sim_mmseqs'] = ['--min-seq-id', '--cov-mode', 'clustering_threshold_for_OTUs', 'AA_mode']
        rule_acceptable_params_dict['align_sequences'] = ['additional_params']
        rule_acceptable_params_dict['trim_alignment'] = ['-gt', "-gappyout", 'additional_params']
        rule_acceptable_params_dict['remove_sequences_with_X_percent_gaps'] = ['--max-percentage-gaps']
        rule_acceptable_params_dict['fasttree'] = ['run']
        rule_acceptable_params_dict['iqtree'] = ['run', '-m', 'additional_params']
        rule_acceptable_params_dict['run_metagenomics_workflow'] = ['clusterize', 'clusterize_submission_params', 'HPC_string', 'snakemake_additional_params', 'bowtie2_additional_params', 'anvi_profile_min_percent_identity']

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Set default values for accessible rules and order of rules in config.json file
        self.default_config.update({
            'metagenomes': 'metagenomes.txt',
            'external_genomes': 'external-genomes.txt',
            'hmm_list': 'hmm_list.txt',
            'samples_txt': 'samples.txt',
            'cluster_representative_method': {'method': 'mmseqs'},
            'anvi_run_hmms_hmmsearch': {'threads_genomes': 1, 'threads_metagenomes': 5},
            'filter_hmm_hits_by_model_coverage': {'threads': 5, '--model-coverage': 0.8, '--filter-out-partial-gene-calls': True},
            'process_hmm_hits': {'threads': 2},
            'combine_sequence_data': {'threads': 2},
            'anvi_get_external_gene_calls_file': {'threads': 2},
            'cat_external_gene_calls_file': {'threads': 2},
            'cluster_X_percent_sim_mmseqs': {'threads': 5, '--min-seq-id': 0.94, '--cov-mode': 0, 'clustering_threshold_for_OTUs': [0.99, 0.98, 0.97], 'AA_mode': False},
            'subset_AA_seqs_with_mmseqs_reps': {'threads': 2},
            'align_sequences': {'threads': 5, 'additional_params': '-maxiters 1 -diags -sv -distance1 kbit20_3'},
            'trim_alignment': {'threads': 5, '-gappyout': True},
            'remove_sequences_with_X_percent_gaps': {'threads': 5, '--max-percentage-gaps': 50},
            'count_num_sequences_filtered': {'threads': 5},
            'subset_DNA_reps_with_QCd_AA_reps_for_mapping': {'threads': 2},
            'subset_external_gene_calls_file_all': {'threads': 2},
            'make_fasta_txt': {'threads': 2},
            'fasttree': {'run': True, 'threads': 5},
            'iqtree': {'threads': 5,'-m': "MFP"},
            'make_metagenomics_config_file': {'threads': 1},
            'run_metagenomics_workflow': {'threads': 2, 'clusterize': False},
            'add_default_collection': {'threads': 2},
            'anvi_summarize': {'threads': 5},
            'rename_tree_tips': {'threads': 1},
            'make_misc_data': {'threads': 2},
            'anvi_scg_taxonomy': {'threads': 5},
            'make_anvio_state_file': {'threads': 2},
            'anvi_import_everything': {'threads': 2},
            'run_genomes_sanity_check': True
            })

        # Directory structure for Snakemake workflow
        self.dirs_dict.update({"HOME": "ECOPHYLO_WORKFLOW"})
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_DIR": os.path.join(self.dirs_dict['HOME'],"01_REFERENCE_PROTEIN_DATA")})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS": os.path.join(self.dirs_dict['HOME'],"02_NR_FASTAS")})
        self.dirs_dict.update({"MSA": os.path.join(self.dirs_dict['HOME'],"03_MSA")})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_STATS": os.path.join(self.dirs_dict['HOME'],"04_SEQUENCE_STATS")})
        self.dirs_dict.update({"TREES": os.path.join(self.dirs_dict['HOME'],"05_TREES")})
        self.dirs_dict.update({"MISC_DATA": os.path.join(self.dirs_dict['HOME'],"06_MISC_DATA")})
        self.dirs_dict.update({"SCG_NT_FASTAS": os.path.join(self.dirs_dict['HOME'],"07_SCG_NT_FASTAS")})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS_RENAMED": os.path.join(self.dirs_dict['HOME'],"08_RIBOSOMAL_PROTEIN_FASTAS_RENAMED")})


    def init(self):
        """This function is called from within the Snakefile to initialize parameters."""
        
        super().init()
        #FIXME: Because 00_LOGS is hardcoded in the base class I need to reassign it
        self.dirs_dict.update({"LOGS_DIR": os.path.join(self.dirs_dict['HOME'],"00_LOGS")})

        # Make log directories
        if not os.path.exists(os.path.join(self.dirs_dict['HOME'], '00_LOGS/')):
            os.makedirs(os.path.join(self.dirs_dict['HOME'], '00_LOGS/'))
        if not os.path.exists(os.path.join(self.dirs_dict['HOME'],'METAGENOMICS_WORKFLOW/00_LOGS/')):
            os.makedirs(os.path.join(self.dirs_dict['HOME'],'METAGENOMICS_WORKFLOW/00_LOGS/'))

        self.names_list = []
        self.contigs_db_name_path_dict = {}
        self.contigs_db_name_bam_dict = {}

        # Load input files
        self.metagenomes = self.get_param_value_from_config(['metagenomes'])
        self.external_genomes = self.get_param_value_from_config(['external_genomes'])
        self.hmm_list_path = self.get_param_value_from_config(['hmm_list'])
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        self.run_genomes_sanity_check = self.get_param_value_from_config(['run_genomes_sanity_check'])

        if not self.metagenomes and not self.external_genomes:
            raise ConfigError('Please provide at least a metagenomes.txt or external-genomes.txt in your '
                              'EcoPhylo config file.')

        if not self.hmm_list_path:
            raise ConfigError('Please provide a path to an hmm_list.txt')

        self.init_hmm_list_txt()

        gene_caller_to_use = self.get_param_value_from_config(['gene_caller_to_use'])
        if not gene_caller_to_use:
            gene_caller_to_use = constants.default_gene_caller

        sanity_checked_metagenomes_file = os.path.join(self.dirs_dict['HOME'], "sanity_checked_metagenomes.txt")
        sanity_checked_genomes_file = os.path.join(self.dirs_dict['HOME'], "sanity_checked_genomes.txt")

        if self.metagenomes:
            filesnpaths.is_file_exists(self.metagenomes)
            self.metagenomes_df = pd.read_csv(self.metagenomes, sep='\t', index_col=False)

            if self.run_genomes_sanity_check:
                if not os.path.exists(sanity_checked_metagenomes_file):
                    args = argparse.Namespace(metagenomes=self.get_param_value_from_config(['metagenomes']), gene_caller = gene_caller_to_use)
                    g = MetagenomeDescriptions(args)
                    g.load_metagenome_descriptions(init=False)
                    self.metagenomes_dict = g.metagenomes_dict
                    self.metagenomes_name_list = list(self.metagenomes_dict.keys())
                    self.metagenomes_path_list = [value['contigs_db_path'] for key,value in self.metagenomes_dict.items()]

                    with open(sanity_checked_metagenomes_file, 'w') as fp:
                        pass
                else:
                    self.run.warning(f"You have declared run_genomes_sanity_check == false. anvi'o takes no responsibility "
                                     f"for any genomes or metagenomes that cause issues downstream in ecophylo.")
                    self.metagenomes_name_list = self.metagenomes_df.name.to_list()
                    self.metagenomes_path_list = self.metagenomes_df.contigs_db_path.to_list()
            else:
                self.metagenomes_name_list = self.metagenomes_df.name.to_list()
                self.metagenomes_path_list = self.metagenomes_df.contigs_db_path.to_list()

            self.contigs_db_name_path_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_path_list)))

            if 'bam' in self.metagenomes_df.columns:
                self.contigs_db_name_bam_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_df.bam)))
                self.metagenomes_profiles_list = self.metagenomes_df.bam.to_list()

            self.names_list.extend(self.metagenomes_name_list)

        else:
            self.metagenomes_name_list = []
        
        if self.external_genomes:
            filesnpaths.is_file_exists(self.external_genomes)
            self.external_genomes_df = pd.read_csv(self.external_genomes, sep='\t', index_col=False)
            
            if self.run_genomes_sanity_check:
                if not os.path.exists(sanity_checked_genomes_file):
                    # FIXME: metagenomes.txt or external-genomes.txt with multiple gene-callers will break
                    # here. Users should only have one type of gene-caller e.g. "NCBI_PGAP".

                    args = argparse.Namespace(external_genomes=self.external_genomes, gene_caller = gene_caller_to_use)
                    genome_descriptions = GenomeDescriptions(args)
                    genome_descriptions.load_genomes_descriptions(init=False)
                    self.external_genomes_dict = genome_descriptions.external_genomes_dict
                    self.external_genomes_names_list = list(self.external_genomes_dict.keys())
                    self.external_genomes_path_list = [value['contigs_db_path'] for key,value in self.external_genomes_dict.items()]

                    with open(sanity_checked_genomes_file, 'w') as fp:
                        pass
                else:
                    self.external_genomes_names_list = self.external_genomes_df.name.to_list()
                    self.external_genomes_path_list = self.external_genomes_df.contigs_db_path.to_list()
            else:
                self.external_genomes_names_list = self.external_genomes_df.name.to_list()
                self.external_genomes_path_list = self.external_genomes_df.contigs_db_path.to_list()

            self.contigs_db_name_path_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_path_list)))

            if 'bam' in self.external_genomes_df.columns:
                self.contigs_db_name_bam_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_df.bam)))
                self.external_genomes_profiles_list = self.external_genomes_df.bam.to_list()

            self.names_list.extend(self.external_genomes_names_list)

        else:
            self.external_genomes_names_list = []

        # Make variables that tells whether we have metagenomes.txt, external-genomes.txt, or both
        if self.metagenomes and not self.external_genomes:
            self.mode = 'metagenomes'
        if not self.metagenomes and self.external_genomes:
            self.mode = 'external_genomes'
        if self.metagenomes and self.external_genomes:
            self.mode = 'both'

        self.AA_mode = self.get_param_value_from_config(['cluster_X_percent_sim_mmseqs', 'AA_mode'])

        if self.samples_txt_file:
            filesnpaths.is_file_tab_delimited(self.samples_txt_file)

            try:
                # getting the samples information (names, [group], path to r1, path to r2) from samples.txt
                self.samples_information = pd.read_csv(self.samples_txt_file, sep='\t', index_col=False)
            except IndexError as e:
                raise ConfigError(f"Looks like your samples_txt file, '%s', is not properly formatted. "
                                  f"This is what we know: {self.samples_txt_file}")
            if 'sample' not in list(self.samples_information.columns):
                raise ConfigError(f"Looks like your samples_txt file, {self.samples_txt_file}, is not properly formatted. "
                                  f"We are not sure what's wrong, but we can't find a column with title 'sample'.")

            self.sample_names_for_mapping_list = self.samples_information['sample'].to_list()

            if self.AA_mode == True:
                raise ConfigError("You provided a samples.txt so you're in profile mode! Please change AA_mode to false.")

        else:
            self.run.warning(f"Since you did not provide a samples.txt, EcoPhylo will assume you do not want "
                             f"to profile the ecology of your proteins and will just be making trees for now!")
        
        # Pick which tree algorithm
        self.run_iqtree = self.get_param_value_from_config(['iqtree', 'run'])
        self.run_fasttree = self.get_param_value_from_config(['fasttree', 'run'])

        if not self.run_iqtree and not self.run_fasttree:
            raise ConfigError("Please choose either iqtree or fasttree in your config file to run your phylogenetic tree.")

        # HPC submission of metagenomics workflow of EcoPhylo
        self.clusterize_metagenomics_workflow = self.get_param_value_from_config(['run_metagenomics_workflow', 'clusterize'])
        self.clusterize_metagenomics_submission_params = self.get_param_value_from_config(['run_metagenomics_workflow', 'clusterize_submission_params'])
        self.metagenomics_workflow_HPC_string = self.get_param_value_from_config(['run_metagenomics_workflow', 'HPC_string'])
        self.metagenomics_workflow_snakemake_additional_params = self.get_param_value_from_config(['run_metagenomics_workflow', 'snakemake_additional_params'])

        self.bowtie2_additional_params = self.get_param_value_from_config(['run_metagenomics_workflow','bowtie2_additional_params'])
        self.anvi_profile_min_percent_identity = self.get_param_value_from_config(['run_metagenomics_workflow','anvi_profile_min_percent_identity'])

        metagenomics_workflow_snakemake_additional_params_list = self.metagenomics_workflow_snakemake_additional_params.split(' ')

        if self.clusterize_metagenomics_workflow:
            if not self.metagenomics_workflow_snakemake_additional_params:
                raise ConfigError("If you are going to use Evan's 'clusterize' (https://github.com/ekiefl/clusterize) with  "
                                  "the metagenomics workflow aspect of EcoPhylo you must provide the '--jobs' in snakemake_additional_params "
                                  "so Snakemake knows how many jobs to run at the same time on the HPC.")

            jobs_param = any("--jobs" in param for param in metagenomics_workflow_snakemake_additional_params_list)
            if jobs_param == False:
                raise ConfigError("The EcoPhylo workflow did not detect the parameter '--jobs' in `snakemake_additional_params`. "
                                  "Please include '--jobs'. You can read about it with snakemake -h")
            
            if self.metagenomics_workflow_HPC_string:
                raise ConfigError("You can't clusterize and provide an HPC_string for the metagenomics workflow at the same time. "
                                  "Please choose one or the other. ")

        # Pick clustering method
        self.cluster_representative_method = self.get_param_value_from_config(['cluster_representative_method', 'method'])

        if self.cluster_representative_method not in ['mmseqs', 'cluster_rep_with_coverages']:
            raise ConfigError(f"anvi'o has never heard of this method to pick a cluster representative: {self.cluster_representative_method} "
                              f"Please check your config file {self.config_file} and change cluster_representative_method to one of the following: 'mmseqs' and 'cluster_rep_with_coverages'")

        if self.cluster_representative_method == 'cluster_rep_with_coverages' and len(self.contigs_db_name_bam_dict) == 0:
            raise ConfigError(f"The EcoPhylo workflow can't use the cluster representative method cluster_rep_with_coverages without BAM files..."
                              f"Please edit your metagenomes.txt or external-genomes.txt and add BAM files.")

        if self.cluster_representative_method == 'cluster_rep_with_coverages' and self.AA_mode == True:
            raise ConfigError(f"The EcoPhylo workflow can't use the cluster representative method cluster_rep_with_coverages in AA_mode")

        # Parse clustering parameter space
        self.clustering_param_space = self.get_param_value_from_config(['cluster_X_percent_sim_mmseqs', 'clustering_threshold_for_OTUs'])
        for count, number in enumerate(self.clustering_param_space):
            if type(number) != float:
                raise ConfigError(f"Element number {count} ({number}) in your clustering_threshold_for_OTUs argument does not appear to be an float. Please provide "
                                  f"a clustering threshold in decimal format i.e. 90% as 0.90")
            if number > 1:
                raise ConfigError(f"The number {number} in your clustering_threshold_for_OTUs argument is not less that one. Please provide "
                                  f"a clustering threshold in decimal format i.e. 90% as 0.90")
        self.clustering_param_space_list_strings = [str(format(clustering_threshold, '.2f')).split(".")[1] + "_percent" for clustering_threshold in self.clustering_param_space]
        self.clustering_threshold_dict = dict(zip(self.clustering_param_space_list_strings, self.clustering_param_space))

        self.target_files = self.get_target_files()


    def get_target_files(self):
        """This function creates a list of target files for Snakemake
        
        RETURNS
        =======
        target_files: list
            list of target files for snakemake
        """

        target_files = []

        for hmm in self.hmm_dict.keys():

            # Clustering parameter space
            for clustering_threshold in self.clustering_param_space_list_strings:
                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], f"{hmm}", f"{clustering_threshold}", f"{hmm}-{clustering_threshold}-mmseqs_NR_rep_seq.fasta")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], f"{hmm}", f"{clustering_threshold}", f"{hmm}-{clustering_threshold}-mmseqs_NR_cluster.tsv")
                target_files.append(target_file)

            if not self.samples_txt_file:
                # TREE-MODE
                target_file = os.path.join(self.dirs_dict['TREES'], f"{hmm}", f"{hmm}_renamed.nwk")
                target_files.append(target_file)
                
                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], f"{hmm}", f"{hmm}_stats.tsv")
                target_files.append(target_file)
                
                target_file = os.path.join(self.dirs_dict['HOME'], f"{hmm}_anvi_estimate_scg_taxonomy_for_SCGs.done")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['HOME'], f"{hmm}_state_imported_tree.done")
                target_files.append(target_file)
            
            else:
                # PROFILE-MODE
                target_file = os.path.join(self.dirs_dict['HOME'], f"{hmm}_state_imported_profile.done")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['TREES'], f"{hmm}", f"{hmm}_renamed.nwk")
                target_files.append(target_file)
                
                target_file = os.path.join(self.dirs_dict['HOME'], f"{hmm}_anvi_estimate_scg_taxonomy_for_SCGs.done")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], f"{hmm}", f"{hmm}_stats.tsv")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "07_SUMMARY", f"{hmm}_summarize.done")
                target_files.append(target_file)
        
        return target_files
    
    def init_hmm_list_txt(self):
        """This function will sanity check hmm-list.txt

        PARAMETERS
        ==========
        self.hmm_list_path : str
            Path to hmm_list.txt

        RETURNS
        =======
        self.hmm_dict : dict
            Dict with hmm as primary key and values: hmm_source, PATH
        """
        filesnpaths.is_file_exists(self.hmm_list_path)
        filesnpaths.is_file_tab_delimited(self.hmm_list_path)
        
        try:
            hmm_df = pd.read_csv(self.hmm_list_path, sep='\t', index_col=False)
        except AttributeError as e:
            raise ConfigError(f"The hmm_list.txt file, {self.hmm_list_path}, does not appear to be properly formatted. "
                              f"This is the error from trying to load it: {self.hmm_list_path}")

        hmm_list_txt_columns = ['name', 'source', 'path']

        for column_name in hmm_list_txt_columns:
            if column_name not in list(hmm_df.columns):
                raise ConfigError(f"Looks like your hmm-list.txt file, {self.hmm_list_path}, is not properly formatted. "
                                  f"We are not sure what's wrong, but we can't find a column with title '{column_name}'."
                                  f"Please make sure you have a tsv with the column names: {hmm_list_txt_columns}")

        self.hmm_dict = hmm_df.set_index('name').to_dict('index')

        if any("-" in s for s in self.hmm_dict.keys()):
            raise ConfigError(f"Please do not use "-" in your external hmm names in: "
                              f"{self.hmm_list_path}. It will make our lives "
                              f"easier with Snakemake wildcards :)")

        # FIXME: this line prints the list of hmm_sources to stdout and I don't want that
        self.internal_hmm_sources = list(anvio.data.hmm.sources.keys())

        for hmm in self.hmm_dict.keys():
            hmm_source = self.hmm_dict[hmm]['source']
            hmm_path = self.hmm_dict[hmm]['path']

            if hmm_path == "INTERNAL":
                if hmm_source not in self.internal_hmm_sources:
                    raise ConfigError(f"{hmm_source} is not an 'INTERNAL' hmm source for anvi'o. "
                                      f"Please double check {self.hmm_list_path} to see if you spelled it right or "
                                      f"please checkout the default internal hmms here: https://merenlab.org/software/anvio/help/7/artifacts/hmm-source/#default-hmm-sources")

            if not filesnpaths.is_file_exists(hmm_path, dont_raise=True):
                if hmm_path == 'INTERNAL':
                    pass
                else:
                    raise ConfigError(f"The path to your hmm {hmm} does not exist: {hmm_path}. "
                                      f"Please double check the paths in our hmm-list.txt: {self.hmm_list_path} "
                                      f"If the hmm you want to use is in an internal anvi'o hmm collection e.g. Bacteria_71 "
                                      f"please put 'INTERNAL' for the path.")
            
            if hmm_path != "INTERNAL":
                sources = u.get_HMM_sources_dictionary([hmm_path])
                
                for source,value in sources.items():
                    gene = value['genes']
                    if hmm_source != source:
                        raise ConfigError(f"In your {self.hmm_list_path}, please change the source for gene {hmm} to this: {source}")
                    if len(gene) > 1:
                        raise ConfigError("EcoPhylo can only work with one gene at a time in a hmm directory (at the moment)")
                    if hmm != gene[0]:
                        raise ConfigError(f"In your {self.hmm_list_path}, please change the gene name {hmm} to this: {gene[0]}")