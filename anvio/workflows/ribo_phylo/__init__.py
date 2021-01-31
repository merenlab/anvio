# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes to define and work with anvi'o trnaseq workflows. """


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


# class RibosomalPhylogeneticsWorkflow(MetagenomicsWorkflow, WorkflowSuperClass):
class RibosomalPhylogeneticsWorkflow(WorkflowSuperClass):

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='ribo_phylo')

        # MetagenomicsWorkflow.__init__(self)
        # ribo_phylo Snakemake rules
        self.rules.extend(['anvi_get_sequences_for_hmm_hits_ribosomal_proteins',
                           'anvi_estimate_scg_taxonomy_for_ribosomal_proteins',
                           'filter_for_scg_sequences_and_metadata',
                           'cat_ribo_proteins_to_one_fasta',
                           'anvi_reformat_fasta_ribosomal_protein_file',
                           'cat_misc_data_to_one_file',
                           'join_renamed_fasta_with_misc_data',
                           'remove_redundant_sequences_mmseqs',
                           'align_muscle',
                           'trim_alignment',
                           'remove_sequences_with_X_percent_gaps',
                           'get_gap_count_distribution',
                           'filter_out_outlier_sequences',
                           'align_muscle_2',
                           'trim_alignment_2',
                           'calculate_tree',
                           'anvi_get_sequences_for_gene_calls',
                           'cluster_90_mmseqs'
                           ])

        self.general_params.extend(['metagenomes']) # user needs to input a metagenomes.txt file
        self.general_params.extend(['external_genomes']) # user can add isolate genomes if needed
        self.general_params.extend(['Ribosomal_protein_list']) # user must input which Ribosomal proteins will be used for workflow
        self.general_params.extend(['MSA_gap_threshold']) # user can input a num gaps threshold to filter the SCG MSA


        # Parameters for each rule that are accessible in the config file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['anvi_get_sequences_for_hmm_hits_ribosomal_proteins'] = ['--hmm-source']
        rule_acceptable_params_dict['anvi_estimate_scg_taxonomy_for_ribosomal_proteins'] = ['--metagenome-mode']
        rule_acceptable_params_dict['anvi_reformat_fasta_ribosomal_protein_file'] = ['--gzip-output',
                                                                                     '--simplify-names',
                                                                                     '--keep-ids',
                                                                                     '--exclude-ids',
                                                                                     '--min-len',
                                                                                     "--prefix"]
        rule_acceptable_params_dict['remove_redundant_sequences_mmseqs'] = ['--min-seq-id']
        rule_acceptable_params_dict['cluster_90_mmseqs'] = ['--min-seq-id']
        rule_acceptable_params_dict['trim_alignment'] = ['-gt', "-gappyout", 'additional_params']
        rule_acceptable_params_dict['trim_alignment_2'] = ['-gt', "-gappyout", 'additional_params']
        rule_acceptable_params_dict['remove_sequences_with_X_percent_gaps'] = ['--max-percentage-gaps']
        rule_acceptable_params_dict['filter_out_outlier_sequences'] = ['-M']
        rule_acceptable_params_dict['calculate_tree'] = ['run',
                                                         '-m',
                                                         'additional_params']

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Set default values for certain accessible parameters
        self.default_config.update({
            'metagenomes': 'metagenomes.txt',
            'external_genomes': 'external-genomes.txt',
            'anvi_reformat_fasta_ribosomal_protein_file': {'--simplify-names': True, 'threads': 5},
            'Ribosomal_protein_list': 'Ribosomal_protein_list.txt',
            'MSA_gap_threshold': '',
            'anvi_estimate_scg_taxonomy_for_ribosomal_proteins': {'threads': 5, '--metagenome-mode': True},
            'filter_for_scg_sequences_and_metadata': {'threads': 5},
            'cat_ribo_proteins_to_one_fasta': {'threads': 5},
            'anvi_get_sequences_for_hmm_hits_ribosomal_proteins': {'threads': 5, '--hmm-source': 'Bacteria_71'},
            'join_renamed_fasta_with_misc_data': {'threads': 5},
            'remove_redundant_sequences_mmseqs': {'threads': 5, '--min-seq-id': 1},
            'cluster_90_mmseqs': {'threads': 5, '--min-seq-id': 0.9},
            'align_muscle': {'threads': 5},
            'remove_sequences_with_X_percent_gaps': {'threads': 5, '--max-percentage-gaps': 50},
            'get_gap_count_distribution': {'threads': 5},
            'filter_out_outlier_sequences': {'threads': 5},
            'align_muscle_2': {'threads': 5},
            'trim_alignment': {'threads': 5, '-gappyout': True},
            'trim_alignment_2': {'threads': 5, '-gappyout': True},
            'calculate_tree': {'run': True, 'threads': 5,'-m': "MFP"}
            })

        # Workflow directory structure

        # The magical line that will put the whole metagenomics workflow into a dir to keep things organized, thanks Sam!
        self.dirs_dict = {directory: 'METAGENOMICS_WORKFLOW/' + dir_path for directory,dir_path in self.dirs_dict.items()}

        # Adding directories specific to Ribo_phylo workflow
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_DIR": 'RIBO_PHYLO_WORKFLOW/01_SCG_HMM_HITS'})
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_TAXONOMY_DIR": "RIBO_PHYLO_WORKFLOW/02_SCG_TAXONOMY"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS": "RIBO_PHYLO_WORKFLOW/03_NR_FASTAS"})
        self.dirs_dict.update({"MSA": "RIBO_PHYLO_WORKFLOW/04_MSA"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_STATS": "RIBO_PHYLO_WORKFLOW/05_SEQUENCE_STATS"})
        self.dirs_dict.update({"TREES": "RIBO_PHYLO_WORKFLOW/06_TREES"})
        self.dirs_dict.update({"MISC_DATA": "RIBO_PHYLO_WORKFLOW/07_MISC_DATA"})
        self.dirs_dict.update({"SCG_NT_FASTAS": "RIBO_PHYLO_WORKFLOW/08_SCG_NT_FASTAS"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS_RENAMED": "RIBO_PHYLO_WORKFLOW/9_RIBOSOMAL_PROTEIN_FASTAS_RENAMED"})


    def init(self):
        """This function is called from within the snakefile to initialize parameters."""

        super().init()
        # WorkflowSuperclass().init()

        # initiating a list to fill with names of contigDBs
        self.names_list = []
        self.names_dirs = []

        # Load metagenomes.txt
        self.metagenomes = self.get_param_value_from_config(['metagenomes'])
        if self.metagenomes:
            filesnpaths.is_file_exists(self.metagenomes)
            try:
                self.metagenomes_df = pd.read_csv(self.metagenomes, sep='\t', index_col=False)
                self.metagenomes_name_list = self.metagenomes_df.name.to_list()
                self.metagenomes_path_list = self.metagenomes_df.contigs_db_path.to_list()
                self.metagenomes_dirname_list = [os.path.dirname(x) for x in self.metagenomes_path_list]
                self.names_dirs.extend(self.metagenomes_dirname_list)
                self.names_list.extend(self.metagenomes_name_list)

            except IndexError as e:
                raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
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
                self.names_dirs.extend(self.external_genomes_dirname_list)
                self.names_list.extend(self.external_genomes_names_list)
            except IndexError as e:
                raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                                  "This is the error from trying to load it: '%s'" % (self.external_genomes_df, e))

        # Make a unique list
        self.names_dirs = list(set(self.names_dirs))

        # Make variables that tells whether we have metagenomes.txt, external-genomes.txt, or both
        if self.metagenomes and not self.external_genomes:
            self.mode = 'metagenomes'
        if not self.metagenomes and self.external_genomes:
            self.mode = 'external_genomes'
        if self.metagenomes and self.external_genomes:
            self.mode = 'both'
        # Load Ribosomal protein list
        self.Ribosomal_protein_list_path = self.get_param_value_from_config(['Ribosomal_protein_list'])
        filesnpaths.is_file_exists(self.Ribosomal_protein_list_path)
        try:
            self.Ribosomal_protein_df = pd.read_csv(self.Ribosomal_protein_list_path, sep='\t', index_col=False)
            self.Ribosomal_protein_list = self.Ribosomal_protein_df['Ribosomal_protein'].to_list()
        except IndexError as e:
            raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.Ribosomal_protein_df, e))

        self.target_files = self.get_target_files()

    def get_target_files(self):
        target_files = []

        # FIXME: need to add target files for metagenomics workflow here!
        # target_files.extend(list(self.profile_databases.values()))

        for ribosomal_protein_name in self.Ribosomal_protein_list:


            # Count num sequences removed per step
            tail_path = "%s_stats.tsv" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

            # tail_path = "%s_all.fna" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], ribosomal_protein_name, tail_path)
            # target_files.append(target_file)

            # # Get SCG taxonomy for each metagenome or external-genome 
            # for external_genome_name in self.names_list:
            #     tail_path = "%s_%s_estimate_scg_taxonomy_results.tsv" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['EXTRACTED_RIBO_PROTEINS_TAXONOMY_DIR'], external_genome_name, tail_path)
            #     target_files.append(target_file)

            # Get final misc data for anvi-interactive display of tree
            tail_path = "%s_all_misc_data_final.tsv" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['MISC_DATA'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

            # Get fasta of nt SCGs for mapping
            tail_path = "%s_scgs_for_mapping.fna" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

            # Get external-gene-calls file for fasta of nt SCGs for mapping
            tail_path = "%s_external_gene_calls_all_renamed.tsv" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

            #########################
            # The FINAL trees :)
            # For iq-tree
            # tail_path = "%s.iqtree" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['TREES'], ribosomal_protein_name, tail_path)
            # target_files.append(target_file)
            # for fasttree
            tail_path = "%s.nwk" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['TREES'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)
            #########################

            # print(target_files)
        return target_files
