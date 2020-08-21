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
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


run = terminal.Run()


class TRNASeqWorkflow(WorkflowSuperClass):
    known_split_types = ['demethylase', 'untreated']

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='trnaseq')

        self.rules.extend(['make_iu_input', # trnaseq Snakemake rules
                           'iu_gen_configs',
                           'iu_merge_pairs',
                           'gen_qc_report',
                           'anvi_reformat_fasta',
                           'anvi_trnaseq'])

        self.general_params.extend(['samples_txt']) # general section of config file

        # Parameters for each rule that are accessible in the config file
        rule_acceptable_params_dict = {}
        rule_acceptable_params_dict['iu_gen_configs'] = ['--r1-prefix',
                                                         '--r2-prefix']
        rule_acceptable_params_dict['iu_merge_pairs'] = ['run',
                                                         '--gzip-output',
                                                         '--marker-gene-stringent',
                                                         '--max-num-mismatches',
                                                         '--report-r1-prefix',
                                                         '--report-r2-prefix']
        rule_acceptable_params_dict['anvi_reformat_fasta'] = ['run',
                                                              '--gzip-output',
                                                              '--simplify-names']
        rule_acceptable_params_dict['anvi_trnaseq'] = ['run',
                                                       '--skip-fasta-check',
                                                       '--write-buffer-size']
        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Default values for certain accessible parameters
        self.default_config.update({
            'samples_txt': 'samples.txt',
            'iu_merge_pairs': {'run': True,
                               '--gzip-output': False,
                               '--marker-gene-stringent': True,
                               '--max-num-mismatches': 0,
                               '--report-r1-prefix': False,
                               '--report-r2-prefix': False,
                               'threads': 1},
            'anvi_reformat_fasta': {'run': True,
                                    '--gzip-output': False,
                                    '--simplify-names': True},
            'anvi_trnaseq': {'run': True,
                             '--skip-fasta-check': True}})

        self.dirs_dict.update({'QC_DIR': '01_QC', 'IDENT_DIR': '02_IDENT'})


    def init(self):
        """ This function is called from within the snakefile to initialize parameters. """
        super().init()

        self.run_iu_merge_pairs = self.get_param_value_from_config(['iu_merge_pairs', 'run'])
        self.run_anvi_reformat_fasta = self.get_param_value_from_config(['anvi_reformat_fasta', 'run'])
        self.run_anvi_trnaseq = self.get_param_value_from_config(['anvi_trnaseq', 'run'])
        self.gzip_iu_merge_pairs_output = self.get_param_value_from_config(['iu_merge_pairs', '--gzip-output'])
        self.gzip_anvi_reformat_fasta_output = self.get_param_value_from_config(['anvi_reformat_fasta', '--gzip-output'])

        # Load table of sample info from samples_txt (sample names, split types, paths to r1, r2)
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(self.samples_txt_file)
        try:
            # An error will subsequently be raised in `check_samples_txt` if there is no header
            self.sample_info = pd.read_csv(self.samples_txt_file, sep='\t', index_col=False)
        except IndexError as e:
            raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.samples_txt_file, e))
        self.check_samples_txt()

        self.sample_names = self.sample_info['sample'].tolist()
        if self.run_iu_merge_pairs:
            self.r1_paths = self.sample_info['r1'].tolist()
            self.r2_paths = self.sample_info['r2'].tolist()
            self.fasta_paths = None
        else:
            self.r1_paths = None
            self.r2_paths = None
            self.fasta_paths = self.sample_info['fasta'].tolist()

        # Determine which samples have which splits
        self.sample_splits_dict = dict([(sample_name, []) for sample_name in self.sample_names])
        for sample_name, split_type in zip(self.sample_names, self.sample_info['split']):
            self.sample_splits_dict[sample_name].append(split_type)

        self.sample_split_prefixes = [sample_name + '_' + split_type
                                      for sample_name, split_type
                                      in zip(self.sample_names, self.sample_info['split'])]

        self.target_files = self.get_target_files()


    def get_target_files(self):
        target_files = []

        for sample_split_prefix in self.sample_split_prefixes:
            output_dir = os.path.join(self.dirs_dict['IDENT_DIR'], sample_split_prefix)
            target_files.append(os.path.join(output_dir, sample_split_prefix + "-TRNASEQ.db"))

        if self.run_iu_merge_pairs:
            target_files.append(os.path.join(self.dirs_dict['QC_DIR'], "qc_report.txt"))

        if self.run_anvi_reformat_fasta:
            for sample_split_prefix in self.sample_split_prefixes:
                target_files.append(os.path.join(self.dirs_dict['QC_DIR'], sample_split_prefix + "-reformat_report.txt"))

        return target_files


    def check_samples_txt(self):

        if self.run_iu_merge_pairs:
            proper_header = ['sample', 'split', 'r1', 'r2']
        else:
            proper_header = ['sample', 'split', 'fasta']
        missing_columns = []
        for column_title in proper_header:
            if column_title not in self.sample_info.columns:
                missing_columns.append(column_title)
        if missing_columns:
            raise ConfigError("The samples_txt file, '%s', is not properly formatted, "
                              "as the following columns are missing: '%s'."
                              % (self.sample_info, ', '.join(missing_columns)))

        for sample_name in self.sample_info['sample']:
            try:
                u.check_sample_id(sample_name)
            except ConfigError as e:
                raise ConfigError("While processing the samples_txt file, '%s', "
                                  "Anvi'o ran into the following error: %s" % (self.samples_txt_file, e))

        unknown_split_types = []
        for split_type in self.sample_info['split']:
            if split_type not in TRNASeqWorkflow.known_split_types:
                unknown_split_types.append(split_type)
        if unknown_split_types:
            run.warning("Some of the names of split types in the samples_txt file, '%s', "
                        "are not what we were expecting (%s). "
                        "That's okay, but Anvi'o decided it should warn you. "
                        "Here are the names of split types that are not in our little list: %s. " % (
                            self.samples_txt_file,
                            ', '.join(TRNASeqWorkflow.known_split_types),
                            ', '.join(sorted(set(unknown_split_types)))))

        if self.run_iu_merge_pairs:
            fastq_paths = self.sample_info['r1'].tolist() + self.sample_info['r2'].tolist()
            bad_fastq_paths = [s for s in fastq_paths if not filesnpaths.is_file_exists(s, dont_raise=True)]
            if bad_fastq_paths:
                raise ConfigError("The following FASTQ files in the samples_txt file, '%s', cannot be found: %s."
                                  % (self.samples_txt_file, ', '.join(bad_fastq_paths)))
            bad_fastq_names = [
                s for s in fastq_paths
                if (not s.endswith('.fq')
                    and not s.endswith('.fq.gz')
                    and not s.endswith('.fastq')
                    and not s.endswith('.fastq.gz'))]
            if bad_fastq_names:
                run.warning("Some of the sequence files in the samples_txt file, '%s', "
                            "do not end with '.fq', '.fq.gz', 'fastq' or '.fastq.gz'. "
                            "That's okay, but Anvi'o decided it should warn you. "
                            "Here are the first 5 such files that have unconventional file extensions: %s."
                            % (self.samples_txt_file, ', '.join(bad_fastq_names[:5])))
        else:
            fasta_paths = self.sample_info['fasta'].tolist()

            bad_fasta_paths = [s for s in fasta_paths if not filesnpaths.is_file_exists(s, dont_raise=True)]
            if bad_fasta_paths:
                raise ConfigError("The following FASTA files in the samples_txt file, '%s', cannot be found: %s."
                                  % (self.samples_txt_file, ', '.join(bad_fasta_paths)))

            bad_fasta_names = [
                s for s in fasta_paths
                if (not s.endswith('.fa')
                    and not s.endswith('.fa.gz')
                    and not s.endswith('.fasta')
                    and not s.endswith('.fasta.gz'))]
            if bad_fasta_names:
                run.warning("Some of the FASTA files in the samples_txt file, '%s', "
                            "do not end with '.fa', '.fa.gz', 'fasta' or '.fasta.gz'. "
                            "That's okay, but Anvi'o decided it should warn you. "
                            "Here are the first 5 such files that have unconventional file extensions: %s."
                            % (self.samples_txt_file, ', '.join(bad_fasta_names[:5])))


    def get_input_for_anvi_reformat_fasta(self, wildcards):
        """
        Input can come from two possible sources:
        a user-supplied FASTA file
        or the FASTA file of merged reads generated from user-supplied FASTQ files.
        """

        if self.fasta_paths:
            return self.fasta_paths[self.sample_split_prefixes.index(wildcards.sample_split_prefix)]
        return os.path.join(self.dirs_dict['QC_DIR'], wildcards.sample_split_prefix + "_MERGED")


    def get_input_for_anvi_trnaseq(self, wildcards):
        """
        Input can come from two possible sources:
        a FASTA file with Anvi'o-compliant deflines supplied by the user
        or the reformatted FASTA file produced by the rule, anvi_reformat_fasta.
        """
        if self.run_anvi_reformat_fasta:
            return os.path.join(self.dirs_dict['QC_DIR'], wildcards.sample_split_prefix + "-reformatted.fasta")
        return self.fasta_paths[self.sample_split_prefixes.index(wildcards.sample_split_prefix)]