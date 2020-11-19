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


class RibosomalPhylogeneticsWorkflow(MetagenomicsWorkflow, WorkflowSuperClass):
# class RibosomalPhylogeneticsWorkflow(WorkflowSuperClass):

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='ribo_phylo')

        MetagenomicsWorkflow.__init__(self)
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

        self.general_params.extend(['external_genomes']) # general section of config file
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
            'external_genomes': 'external_genomes.txt',
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

        # Added directories in the workflow

        # 01_SCG_HMM_HITS
        # 02_NR_FASTA
        # 03_MSA
        # 04_TREE

        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_DIR": "01_SCG_HMM_HITS"})
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_TAXONOMY_DIR": "02_SCG_TAXONOMY"})
        self.dirs_dict.update({"FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY_DIR": "03_FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS": "04_NR_FASTAS"})
        self.dirs_dict.update({"MSA": "05_MSA"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_STATS": "06_SEQUENCE_STATS"})
        self.dirs_dict.update({"TREES": "07_TREES"})
        self.dirs_dict.update({"MISC_DATA": "08_MISC_DATA"})
        self.dirs_dict.update({"SCG_NT_FASTAS": "09_SCG_NT_FASTAS"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS_RENAMED": "10_RIBOSOMAL_PROTEIN_FASTAS_RENAMED"})


    def init(self):
        """This function is called from within the snakefile to initialize parameters."""

        super().init()

    #     self.run_iu_merge_pairs = self.get_param_value_from_config(['iu_merge_pairs', 'run'])
    #     self.run_anvi_reformat_fasta = self.get_param_value_from_config(['anvi_reformat_fasta', 'run'])
    #     self.run_anvi_trnaseq = self.get_param_value_from_config(['anvi_trnaseq', 'run'])
    #     self.gzip_iu_merge_pairs_output = self.get_param_value_from_config(['iu_merge_pairs', '--gzip-output'])
    #     self.gzip_anvi_reformat_fasta_output = self.get_param_value_from_config(['anvi_reformat_fasta', '--gzip-output'])

        # Load table of sample info from samples_txt (sample names, split types, paths to r1, r2).
        self.external_genomes = self.get_param_value_from_config(['external_genomes'])
        filesnpaths.is_file_exists(self.external_genomes)
        try:
            # An error will subsequently be raised in `check_samples_txt` if there is no header.
            self.external_genomes_df = pd.read_csv(self.external_genomes, sep='\t', index_col=False)
            self.external_genome_name_list = self.external_genomes_df.name.to_list()
            self.external_genome_path_list = self.external_genomes_df.contigs_db_path.to_list()
            # Dict for when the time comes to tell SCG_taxonomy if it is working with a genome for metagnome
            # add a new column to the external_genomes.txt called type to denote genome or metagenome
            self.external_genome_type_dict = dict(zip(self.external_genomes_df.name, self.external_genomes_df.type))

        except IndexError as e:
            # FIXME: need to make a better ConfigError
            raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.samples_txt_file, e))

        self.contig_dir = os.path.dirname(self.external_genome_path_list[0])

        # Load Ribosomal protein list
        self.Ribosomal_protein_list_path = self.get_param_value_from_config(['Ribosomal_protein_list'])
        filesnpaths.is_file_exists(self.Ribosomal_protein_list_path)

        try:
            # An error will subsequently be raised in `check_samples_txt` if there is no header.
            self.Ribosomal_protein_df = pd.read_csv(self.Ribosomal_protein_list_path, sep='\t', index_col=False)

        except IndexError as e:
            raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.Ribosomal_protein_df, e))
        self.Ribosomal_protein_list = self.Ribosomal_protein_df['Ribosomal_protein'].to_list()

        self.target_files = self.get_target_files()

    def get_target_files(self):
        target_files = []

        # FIXME: need to add target files for metagenomics workflow here!
        # target_files.extend(list(self.profile_databases.values()))

        for ribosomal_protein_name in self.Ribosomal_protein_list:


            # Num sequences removed per step
            tail_path = "%s_stats.tsv" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_MSA_STATS'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

            # Misc metadata files
            tail_path = "%s_all_misc_data_final.tsv" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['MISC_DATA'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

               # Misc metadata files
            # tail_path = "%s_reformat_report_mmseqs_NR.txt" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, tail_path)
            # target_files.append(target_file)

            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_reps_leeway.fna" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, external_genome_name, tail_path)
            #     print(target_file)
            #     target_files.append(target_file)
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_gene_callers_ids_reps.tsv" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, external_genome_name, tail_path)
            #     print(target_file)
            #     target_files.append(target_file)

            # The FINAL trees :)
            tail_path = "%s.iqtree" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['TREES'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)


            # tail_path = "%s_trimmed.fasta" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['MSA'], "MSA_2", ribosomal_protein_name, tail_path)
            # target_files.append(target_file)

            # #
            # tail_path = "%s_gene_callers_ids_reps.tsv" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY_DIR'], ribosomal_protein_name, tail_path)
            # # print(target_file)
            # target_files.append(target_file)

            # rename_gene_calls
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_renamed.fna" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'],ribosomal_protein_name, external_genome_name, tail_path)
            #     target_files.append(target_file)


            # # cluster_90_mmseqs
            # tail_path = "%s_mmseqs_NR_rep_seq.fasta" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, tail_path)
            # target_files.append(target_file)

            # get_filtered_reformat_file
            # tail_path = "%s_reformat_report_mmseqs_NR.txt" % (ribosomal_protein_name)
            # target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, tail_path)
            # target_files.append(target_file)

            # anvi_get_gene_caller_ids_for_reps
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_gene_callers_ids_reps.tsv" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'],external_genome_name, ribosomal_protein_name, tail_path)
            #     target_files.append(target_file)

            # # anvi_get_gene_caller_ids_for_reps
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_gene_callers_ids_reps.tsv" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, external_genome_name, tail_path)
            #     print(target_file)
            #     target_files.append(target_file)

            # # anvi_get_sequences_for_gene_calls_reps_with_leeway
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_reps_leeway.fna" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, external_genome_name, tail_path)
            #     print(target_file)
            #     target_files.append(target_file)

            # # rename_gene_calls_reps  
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_all_reps_leeway.fna" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, external_genome_name, tail_path)
            #     print(target_file)
            #     target_files.append(target_file)

            # # anvi_reformat_fasta_SCG_NT  
            # for external_genome_name in self.external_genome_name_list:
            #     # Nucleotide fasta
            #     tail_path = "%s_%s_all_reps_leeway_renamed.fna" % (external_genome_name, ribosomal_protein_name)
            #     target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, external_genome_name, tail_path)
            #     print(target_file)
            #     target_files.append(target_file)
           # # cat_SCG_NT_to_one_fasta_reps
           #  tail_path = "%s_all_reps_leeway.fna" % (ribosomal_protein_name)
           #  target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, tail_path)
           #  target_files.append(target_file)

            # cat_SCG_NT_to_one_fasta
            tail_path = "%s_all_leeway.fna" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['SCG_NT_FASTAS'], ribosomal_protein_name, tail_path)
            target_files.append(target_file)

            # from contigs workflow
            # contigs_annotated = [os.path.join(self.dirs_dict["CONTIGS_DIR"],\
            # g + "-annotate_contigs_database.done") for g in self.group_names]
            # target_files.extend(contigs_annotated)


        return target_files
