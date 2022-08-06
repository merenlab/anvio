# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes to define and work with anvi'o ecophylo workflows. """

import os
import sys
import anvio
import argparse
import pandas as pd

import anvio
import anvio.data.hmm
import anvio.utils as u
import anvio.genomedescriptions as gd
import anvio.terminal as terminal
import anvio.workflows as w
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.genomedescriptions import GenomeDescriptions
from anvio.genomedescriptions import MetagenomeDescriptions
from anvio.workflows.metagenomics import MetagenomicsWorkflow


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
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


        # Parameters for each rule that are accessible in the config.json file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['anvi_run_hmms_hmmsearch'] = ['threads_genomes', 'threads_metagenomes', 'additional_params']
        rule_acceptable_params_dict['filter_hmm_hits_by_model_coverage'] = ['--model-coverage', 'clustering_threshold_for_OTUs', 'additional_params']
        rule_acceptable_params_dict['cluster_X_percent_sim_mmseqs'] = ['--min-seq-id', 'clustering_threshold_for_OTUs']
        rule_acceptable_params_dict['align_sequences'] = ['additional_params']
        rule_acceptable_params_dict['trim_alignment'] = ['-gt', "-gappyout", 'additional_params']
        rule_acceptable_params_dict['remove_sequences_with_X_percent_gaps'] = ['--max-percentage-gaps']
        rule_acceptable_params_dict['fasttree'] = ['run']
        rule_acceptable_params_dict['iqtree'] = ['run', '-m', 'additional_params']
        rule_acceptable_params_dict['run_metagenomics_workflow'] = ['clusterize', 'cluster_submission_params', 'snakemake_additional_params']

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Set default values for accessible rules and order of rules in config.json file
        self.default_config.update({
            'metagenomes': 'metagenomes.txt',
            'external_genomes': 'external-genomes.txt',
            'hmm_list': 'hmm_list.txt',
            'samples_txt': 'samples.txt',
            'cluster_representative_method': {'method': 'mmseqs'},
            'anvi_run_hmms_hmmsearch': {'threads_genomes': 1, 'threads_metagenomes': 5},
            'filter_hmm_hits_by_model_coverage': {'threads': 5, '--model-coverage': 0.8},
            'process_hmm_hits': {'threads': 2},
            'combine_sequence_data': {'threads': 2},
            'anvi_get_external_gene_calls_file': {'threads': 5},
            'cat_external_gene_calls_file': {'threads': 2},
            'cluster_X_percent_sim_mmseqs': {'threads': 5, '--min-seq-id': 0.94, 'clustering_threshold_for_OTUs': [0.99, 0.98, 0.97]},
            'subset_AA_seqs_with_mmseqs_reps': {'threads': 2},
            'align_sequences': {'threads': 5, 'additional_params': '-maxiters 1 -diags -sv -distance1 kbit20_3'},
            'trim_alignment': {'threads': 5, '-gappyout': True},
            'remove_sequences_with_X_percent_gaps': {'threads': 5, '--max-percentage-gaps': 50},
            'count_num_sequences_filtered': {'threads': 5},
            'subset_DNA_reps_with_QCd_AA_reps_for_mapping': {'threads': 2},
            'subset_external_gene_calls_file_all': {'threads': 2},
            'make_fasta_txt': {'threads': 5},
            'fasttree': {'run': True, 'threads': 5},
            'iqtree': {'threads': 5,'-m': "MFP"},
            'make_metagenomics_config_file': {'threads': 1},
            'run_metagenomics_workflow': {'threads': 5, 'clusterize': False},
            'add_default_collection': {'threads': 5},
            'anvi_summarize': {'threads': 5},
            'rename_tree_tips': {'threads': 1},
            'make_misc_data': {'threads': 2},
            'anvi_scg_taxonomy': {'threads': 5},
            'make_anvio_state_file': {'threads': 5},
            'anvi_import_everything': {'threads': 5},
            })

        # Directory structure for Snakemake workflow
        self.dirs_dict.update({"HOME": "ECOPHYLO_WORKFLOW"})
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_DIR": "ECOPHYLO_WORKFLOW/01_REFERENCE_PROTEIN_DATA"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS": "ECOPHYLO_WORKFLOW/02_NR_FASTAS"})
        self.dirs_dict.update({"MSA": "ECOPHYLO_WORKFLOW/03_MSA"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_STATS": "ECOPHYLO_WORKFLOW/04_SEQUENCE_STATS"})
        self.dirs_dict.update({"TREES": "ECOPHYLO_WORKFLOW/05_TREES"})
        self.dirs_dict.update({"MISC_DATA": "ECOPHYLO_WORKFLOW/06_MISC_DATA"})
        self.dirs_dict.update({"SCG_NT_FASTAS": "ECOPHYLO_WORKFLOW/07_SCG_NT_FASTAS"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS_RENAMED": "ECOPHYLO_WORKFLOW/08_RIBOSOMAL_PROTEIN_FASTAS_RENAMED"})

        # Make log directories
        if not os.path.exists('ECOPHYLO_WORKFLOW/00_LOGS/'):
            os.makedirs('ECOPHYLO_WORKFLOW/00_LOGS/')
        if not os.path.exists('ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/00_LOGS/'):
            os.makedirs('ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/00_LOGS/')


    def init(self):
        """This function is called from within the snakefile to initialize parameters."""
        
        super().init()
        #FIXME: Because 00_LOGS is hardcoded in the base class I need to reassign it
        self.dirs_dict.update({"LOGS_DIR": "ECOPHYLO_WORKFLOW/00_LOGS"})

        self.names_list = []
        self.contigsDB_name_path_dict = {}
        self.contigsDB_name_bam_dict = {}

        # Load metagenomes.txt
        self.metagenomes = self.get_param_value_from_config(['metagenomes'])

        if self.metagenomes:
            args = argparse.Namespace(metagenomes=self.get_param_value_from_config(['metagenomes']))
            g = MetagenomeDescriptions(args)
            g.load_metagenome_descriptions(skip_sanity_check=True)
            self.metagenomes_dict = g.metagenomes_dict
            self.metagenomes_name_list = list(self.metagenomes_dict.keys())
            self.metagenomes_path_list = [value['contigs_db_path'] for key,value in self.metagenomes_dict.items()]
            self.contigsDB_name_path_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_path_list)))

            self.metagenomes_df = pd.read_csv(self.metagenomes, sep='\t', index_col=False)
            if 'bam' in self.metagenomes_df.columns:
                self.contigsDB_name_bam_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_df.bam)))
                self.metagenomes_profiles_list = self.metagenomes_df.bam.to_list()
            self.names_list.extend(self.metagenomes_name_list)

        # Load external-genomes.txt
        self.external_genomes = self.get_param_value_from_config(['external_genomes'])
        
        if self.external_genomes:
            
            args = argparse.Namespace(external_genomes=self.external_genomes)
            genome_descriptions = GenomeDescriptions(args)
            genome_descriptions.load_genomes_descriptions()
            self.external_genomes_dict = genome_descriptions.external_genomes_dict
            self.external_genomes_names_list = list(self.external_genomes_dict.keys())
            self.external_genomes_path_list = [value['contigs_db_path'] for key,value in self.external_genomes_dict.items()]
            self.contigsDB_name_path_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_path_list)))


            self.external_genomes_df = pd.read_csv(self.external_genomes, sep='\t', index_col=False)
            if 'bam' in self.external_genomes_df.columns:
                self.contigsDB_name_bam_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_df.bam)))
                self.external_genomes_profiles_list = self.external_genomes_df.bam.to_list()
            self.names_list.extend(self.external_genomes_names_list)

        else:
            self.external_genomes_names_list = []

        # Concatenate metagenomes.txt and external-genomes.txt
        contigsDB_name_path_list = list(self.contigsDB_name_path_dict.items())
        contigsDB_name_path_df = pd.DataFrame(contigsDB_name_path_list, columns=['name', 'contigs_db_path'])
        self.combined_genomes_df_path = "ECOPHYLO_WORKFLOW/combined_genomes.txt"
        
        contigsDB_name_path_df.to_csv(self.combined_genomes_df_path, \
                                      sep="\t", \
                                      index=False, \
                                      header=True)

        # Make variables that tells whether we have metagenomes.txt, external-genomes.txt, or both
        if self.metagenomes and not self.external_genomes:
            self.mode = 'metagenomes'
        if not self.metagenomes and self.external_genomes:
            self.mode = 'external_genomes'
        if self.metagenomes and self.external_genomes:
            self.mode = 'both'

        # Load HMM list
        self.hmm_list_path = self.get_param_value_from_config(['hmm_list'])
        if self.hmm_list_path: 
            filesnpaths.is_file_exists(self.hmm_list_path)
            filesnpaths.is_file_tab_delimited(self.hmm_list_path)
            try:
                HMM_df = pd.read_csv(self.hmm_list_path, sep='\t', index_col=False)
                self.HMM_source_dict = dict(zip(HMM_df.name, HMM_df.source))
                self.HMM_path_dict = dict(zip(HMM_df.name, HMM_df.path))

            except AttributeError as e:
                raise ConfigError("The hmm_list.txt file, '%s', does not appear to be properly formatted. "
                                  "This is the error from trying to load it: %s" % (self.hmm_list_path, e))

            if any("-" in s for s in self.HMM_source_dict.keys()):
                raise ConfigError(f"Please do not use "-" in your external HMM names in: "
                                  f"{self.hmm_list_path}. It will make our lives "
                                  f"easier with Snakemake wildcards :)")

        # FIXME: this line prints the list of HMM_sources to stdout and I don't want that
        self.internal_HMM_sources = list(anvio.data.hmm.sources.keys())

        for HMM, HMM_path in self.HMM_source_dict.items():
            if HMM_path == "INTERNAL":
                if self.HMM_source_dict[HMM] not in self.internal_HMM_sources:
                    HMM_source = self.HMM_source_dict[HMM]
                    raise ConfigError(f"{HMM_source} is not an 'INTERNAL' HMM source for anvi'o. "
                                      f"Please double check {self.hmm_list_path} to see if you spelled it right or "
                                      f"please checkout the default internal HMMs here: https://merenlab.org/software/anvio/help/7/artifacts/hmm-source/#default-hmm-sources")

        # Load samples.txt
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        if self.samples_txt_file:
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
        else:
            self.run.warning(f"Since you did not provide a samples.txt, EcoPhylo will assume you do not want "
                             f"to profile the ecology and will just be making trees for now!")
        
        # Pick which tree algorithm
        self.run_iqtree = self.get_param_value_from_config(['iqtree', 'run'])
        self.run_fasttree = self.get_param_value_from_config(['fasttree', 'run'])

        if not self.run_iqtree and not self.run_fasttree:
            raise ConfigError("Please choose either iqtree or fasttree in your config file to run your phylogenetic tree.")

        # Decide to clusterize metagenomic workflow
        self.clusterize_metagenomics_workflow = self.get_param_value_from_config(['run_metagenomics_workflow', 'clusterize'])
        self.metagenomics_workflow_HPC_string = self.get_param_value_from_config(['run_metagenomics_workflow', 'cluster_submission_params'])
        self.snakemake_additional_params = self.get_param_value_from_config(['run_metagenomics_workflow', 'snakemake_additional_params'])

        # Pick clustering method
        self.cluster_representative_method = self.get_param_value_from_config(['cluster_representative_method', 'method'])

        if self.cluster_representative_method not in ['mmseqs', 'cluster_rep_with_coverages']:
            raise ConfigError(f"anvi'o has never heard of this method to pick a cluster representative: {self.cluster_representative_method} "
                              f"Please check your config file {self.config_file} and change cluster_representative_method to one of the following: 'mmseqs' and 'cluster_rep_with_coverages'")

        if self.cluster_representative_method == 'cluster_rep_with_coverages' and len(self.contigsDB_name_bam_dict) == 0:
            raise ConfigError(f"The EcoPhylo workflow can't use the cluster representative method cluster_rep_with_coverages without BAM files..."
                              f"Please edit your metagenomes.txt or external-genomes.txt and add BAM files.")

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
        target_files = []

        for HMM in self.HMM_source_dict.keys():

            # Clustering parameter space
            for clustering_threshold in self.clustering_param_space_list_strings:
                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], f"{HMM}", f"{clustering_threshold}", f"{HMM}-{clustering_threshold}-mmseqs_NR_rep_seq.fasta")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], f"{HMM}", f"{clustering_threshold}", f"{HMM}-{clustering_threshold}-mmseqs_NR_cluster.tsv")
                target_files.append(target_file)

            if not self.samples_txt_file:
                # TREE-MODE
                target_file = os.path.join(self.dirs_dict['TREES'], f"{HMM}", f"{HMM}_renamed.nwk")
                target_files.append(target_file)
                
                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], f"{HMM}", f"{HMM}_stats.tsv")
                target_files.append(target_file)
                
                target_file = os.path.join("ECOPHYLO_WORKFLOW", f"{HMM}_anvi_estimate_scg_taxonomy_for_SCGs.done")
                target_files.append(target_file)

                target_file = os.path.join("ECOPHYLO_WORKFLOW", f"{HMM}_state_imported_tree.done")
                target_files.append(target_file)
            
            else:
                # PROFILE-MODE
                target_file = os.path.join("ECOPHYLO_WORKFLOW", f"{HMM}_state_imported_profile.done")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['TREES'], f"{HMM}", f"{HMM}_renamed.nwk")
                target_files.append(target_file)
                
                target_file = os.path.join("ECOPHYLO_WORKFLOW", f"{HMM}_anvi_estimate_scg_taxonomy_for_SCGs.done")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], f"{HMM}", f"{HMM}_stats.tsv")
                target_files.append(target_file)
        
        return target_files