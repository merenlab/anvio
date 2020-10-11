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


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


run = terminal.Run()


class RibosomalPhylogeneticsWorkflow(WorkflowSuperClass):

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='ribo_phylo')

        # ribo_phylo Snakemake rules
        self.rules.extend(['anvi_get_sequences_for_hmm_hits_ribosomal_proteins', 
                            'anvi_estimate_scg_taxonomy_for_ribosomal_proteins'
                           # 'align_ribosomal_sequences',
                           # 'trim_alignment',
                           # 'filter_sequences_max_percen_gaps_50'
                           ])

        self.general_params.extend(['samples_txt']) # general section of config file
        self.general_params.extend(['Ribosomal_protein_list']) # user must input which Ribosomal proteins will be used for workflow


        # Parameters for each rule that are accessible in the config file
        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['anvi_get_sequences_for_hmm_hits_ribosomal_proteins'] = ['run']
        rule_acceptable_params_dict['anvi_estimate_scg_taxonomy_for_ribosomal_proteins'] = ['run']
        rule_acceptable_params_dict['anvi_reformat_fasta'] = ['run',
                                                              '--gzip-output',
                                                              '--simplify-names']
        rule_acceptable_params_dict['anvi_trnaseq'] = ['run',
                                                       '--overwrite-output-destinations',
                                                       '--description',
                                                       '--skip-fasta-check',
                                                       '--write-checkpoints',
                                                       '--load-checkpoint',
                                                       '--write-buffer-size',
                                                       '--alignment-target-chunk-size',
                                                       '--fragment-mapping-query-chunk-length',
                                                       '--feature-param-file',
                                                       '--min-trna-fragment-size',
                                                       '--agglomeration-max-mismatch-freq',
                                                       '--max-deletion-size',
                                                       '--alignment-progress-interval',
                                                       '--agglomeration-progress-interval']
        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Set default values for certain accessible parameters
        self.default_config.update({
            'samples_txt': 'samples.txt',
            'Ribosomal_protein_list': 'Ribosomal_protein_list.txt',
            'anvi_get_sequences_for_hmm_hits_ribosomal_proteins': {'run': True, 'threads': 5},
            'anvi_estimate_scg_taxonomy_for_ribosomal_proteins': {'run': True, 'threads': 5}})

        # Added directories in the workflow
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_DIR": "01_EXTRACTED_RIBO_PROTEINS"})
        self.dirs_dict.update({"EXTRACTED_RIBO_PROTEINS_TAXONOMY_DIR": "02_EXTRACTED_RIBO_PROTEINS_TAXONOMY"})
        self.dirs_dict.update({"FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY_DIR": "03_FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_FASTAS": "04_RIBOSOMAL_PROTEIN_FASTAS"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_1": "05_RIBOSOMAL_PROTEIN_MSA_1"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_MSA_2": "06_RIBOSOMAL_PROTEIN_MSA_2"})
        self.dirs_dict.update({"RIBOSOMAL_PROTEIN_TREES": "07_RIBOSOMAL_PROTEIN_TREES"})


    def init(self):
        """This function is called from within the snakefile to initialize parameters."""

        super().init()

    #     self.run_iu_merge_pairs = self.get_param_value_from_config(['iu_merge_pairs', 'run'])
    #     self.run_anvi_reformat_fasta = self.get_param_value_from_config(['anvi_reformat_fasta', 'run'])
    #     self.run_anvi_trnaseq = self.get_param_value_from_config(['anvi_trnaseq', 'run'])
    #     self.gzip_iu_merge_pairs_output = self.get_param_value_from_config(['iu_merge_pairs', '--gzip-output'])
    #     self.gzip_anvi_reformat_fasta_output = self.get_param_value_from_config(['anvi_reformat_fasta', '--gzip-output'])

        # Load table of sample info from samples_txt (sample names, split types, paths to r1, r2).
        self.samples_txt_path = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(self.samples_txt_path)
        try:
            # An error will subsequently be raised in `check_samples_txt` if there is no header.
            self.sample_txt = pd.read_csv(self.samples_txt_path, sep='\t', index_col=False)

        except IndexError as e:
            raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.samples_txt_file, e))
        self.sample_name_list = self.sample_txt['name'].to_list()
        self.sample_path_list = self.sample_txt['path'].to_list()
        # 03_FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY/ACE_MG_NA_0150_1482_118C/ACE_MG_NA_0150_1482_118C_Ribosomal_L16.fasta    

        self.contig_dir = os.path.dirname(self.sample_path_list[0])
        # self.snakemake_input = os.path.join(contig_dir, "{sample_name}-contigs.db")


        # self.check_samples_txt()

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
    #     self.sample_names = self.sample_info['sample'].tolist()
    #     if self.run_iu_merge_pairs:
    #         self.r1_paths = self.sample_info['r1'].tolist()
    #         self.r2_paths = self.sample_info['r2'].tolist()
    #         self.fasta_paths = None
    #     else:
    #         self.r1_paths = None
    #         self.r2_paths = None
    #         self.fasta_paths = self.sample_info['fasta'].tolist()

    #     # Determine which samples have which splits.
    #     self.sample_splits_dict = dict([(sample_name, []) for sample_name in self.sample_names])
    #     for sample_name, split_type in zip(self.sample_names, self.sample_info['split']):
    #         self.sample_splits_dict[sample_name].append(split_type)

    #     self.sample_split_prefixes = [sample_name + '_' + split_type
    #                                   for sample_name, split_type
    #                                   in zip(self.sample_names, self.sample_info['split'])]

        self.target_files = self.get_target_files()

    def get_target_files(self):
        target_files = []

        # Target files for anvi_get_sequences_for_hmm_hits_ribosomal_proteins
        # for ribosomal_protein_name in self.Ribosomal_protein_list:
        #     for sample_name in self.sample_name_list:
        #         tail_path = "%s_%s_hmm_hits_all.fasta" % (sample_name, ribosomal_protein_name)
        #         target_file = os.path.join(self.dirs_dict['EXTRACTED_RIBO_PROTEINS_DIR'], sample_name, tail_path)
        #         target_files.append(target_file)

        # # Target files for anvi_estimate_scg_taxonomy_for_ribosomal_proteins
        # for ribosomal_protein_name in self.Ribosomal_protein_list:
        #     for sample_name in self.sample_name_list:
        #         tail_path = "%s_%s_estimate_scg_taxonomy_results.tsv" % (sample_name, ribosomal_protein_name)
        #         target_file = os.path.join(self.dirs_dict['EXTRACTED_RIBO_PROTEINS_TAXONOMY_DIR'], sample_name, tail_path)
        #         target_files.append(target_file)

        # # Target files for filter_for_scg_sequences_and_metadata:
        # for ribosomal_protein_name in self.Ribosomal_protein_list:
        #     for sample_name in self.sample_name_list:
        #         tail_path = "%s_%s_misc_data.tsv" % (sample_name, ribosomal_protein_name)
        #         target_file = os.path.join(self.dirs_dict['FILTERED_RIBO_PROTEINS_SEQUENCES_TAXONOMY_DIR'], sample_name, tail_path)
        #         target_files.append(target_file)

        # # Target files for filter_for_scg_sequences_and_metadata:
        # for ribosomal_protein_name in self.Ribosomal_protein_list:

        #     tail_path = "%s_misc_data_final.tsv" % (ribosomal_protein_name)
        #     target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], tail_path)
        #     target_files.append(target_file)

        # Target files for mmseqs:
        for ribosomal_protein_name in self.Ribosomal_protein_list:

            tail_path = "%s.contree" % (ribosomal_protein_name)
            target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_TREES'], tail_path)
            target_files.append(target_file)
        # # Target files for filter_for_scg_sequences_and_metadata:
        # for ribosomal_protein_name in self.Ribosomal_protein_list:

        #     tail_path = "%s_reformat_report_all.txt" % (ribosomal_protein_name)
        #     target_file = os.path.join(self.dirs_dict['RIBOSOMAL_PROTEIN_FASTAS'], tail_path)
        #     target_files.append(target_file)

        return target_files


    # def check_samples_txt(self):

    #     if self.run_iu_merge_pairs:
    #         proper_header = ['sample', 'split', 'r1', 'r2']
    #     else:
    #         proper_header = ['sample', 'split', 'fasta']
    #     missing_columns = []
    #     for column_title in proper_header:
    #         if column_title not in self.sample_info.columns:
    #             missing_columns.append(column_title)
    #     if missing_columns:
    #         raise ConfigError("The samples_txt file, '%s', is not properly formatted, "
    #                           "as the following columns are missing: '%s'."
    #                           % (self.sample_info, ', '.join(missing_columns)))

    #     for sample_name in self.sample_info['sample']:
    #         try:
    #             u.check_sample_id(sample_name)
    #         except ConfigError as e:
    #             raise ConfigError("While processing the samples_txt file, '%s', "
    #                               "Anvi'o ran into the following error: %s" % (self.samples_txt_file, e))

    #     unknown_split_types = []
    #     for split_type in self.sample_info['split']:
    #         if split_type not in TRNASeqWorkflow.known_split_types:
    #             unknown_split_types.append(split_type)
    #     if unknown_split_types:
    #         run.warning("Some of the names of split types in the samples_txt file, '%s', "
    #                     "are not what we were expecting (%s). "
    #                     "That's okay, but Anvi'o decided it should warn you. "
    #                     "Here are the names of split types that are not in our little list: %s. " % (
    #                         self.samples_txt_file,
    #                         ', '.join(TRNASeqWorkflow.known_split_types),
    #                         ', '.join(sorted(set(unknown_split_types)))))

    #     if self.run_iu_merge_pairs:
    #         fastq_paths = self.sample_info['r1'].tolist() + self.sample_info['r2'].tolist()
    #         bad_fastq_paths = [s for s in fastq_paths if not filesnpaths.is_file_exists(s, dont_raise=True)]
    #         if bad_fastq_paths:
    #             raise ConfigError("The following FASTQ files in the samples_txt file, '%s', cannot be found: %s."
    #                               % (self.samples_txt_file, ', '.join(bad_fastq_paths)))
    #         bad_fastq_names = [
    #             s for s in fastq_paths
    #             if (not s.endswith('.fq')
    #                 and not s.endswith('.fq.gz')
    #                 and not s.endswith('.fastq')
    #                 and not s.endswith('.fastq.gz'))]
    #         if bad_fastq_names:
    #             run.warning("Some of the sequence files in the samples_txt file, '%s', "
    #                         "do not end with '.fq', '.fq.gz', 'fastq' or '.fastq.gz'. "
    #                         "That's okay, but Anvi'o decided it should warn you. "
    #                         "Here are the first 5 such files that have unconventional file extensions: %s."
    #                         % (self.samples_txt_file, ', '.join(bad_fastq_names[:5])))
    #     else:
    #         fasta_paths = self.sample_info['fasta'].tolist()

    #         bad_fasta_paths = [s for s in fasta_paths if not filesnpaths.is_file_exists(s, dont_raise=True)]
    #         if bad_fasta_paths:
    #             raise ConfigError("The following FASTA files in the samples_txt file, '%s', cannot be found: %s."
    #                               % (self.samples_txt_file, ', '.join(bad_fasta_paths)))

    #         bad_fasta_names = [
    #             s for s in fasta_paths
    #             if (not s.endswith('.fa')
    #                 and not s.endswith('.fa.gz')
    #                 and not s.endswith('.fasta')
    #                 and not s.endswith('.fasta.gz'))]
    #         if bad_fasta_names:
    #             run.warning("Some of the FASTA files in the samples_txt file, '%s', "
    #                         "do not end with '.fa', '.fa.gz', 'fasta' or '.fasta.gz'. "
    #                         "That's okay, but Anvi'o decided it should warn you. "
    #                         "Here are the first 5 such files that have unconventional file extensions: %s."
    #                         % (self.samples_txt_file, ', '.join(bad_fasta_names[:5])))


    # def get_input_for_anvi_reformat_fasta(self, wildcards):
    #     """Input can come from two possible sources:
    #     a user-supplied FASTA file
    #     or the FASTA file of merged reads generated from user-supplied FASTQ files.
    #     """

    #     if self.fasta_paths:
    #         return self.fasta_paths[self.sample_split_prefixes.index(wildcards.sample_split_prefix)]
    #     return os.path.join(self.dirs_dict['QC_DIR'], wildcards.sample_split_prefix + "_MERGED")


    # def get_input_for_anvi_trnaseq(self, wildcards):
    #     """Input can come from two possible sources:
    #     a FASTA file with Anvi'o-compliant deflines supplied by the user
    #     or the reformatted FASTA file produced by the rule, anvi_reformat_fasta.
    #     """

    #     if self.run_anvi_reformat_fasta:
    #         return os.path.join(self.dirs_dict['QC_DIR'], wildcards.sample_split_prefix + "-reformatted.fasta")
    #     return self.fasta_paths[self.sample_split_prefixes.index(wildcards.sample_split_prefix)]
