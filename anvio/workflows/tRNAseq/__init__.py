# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o tRNA-seq workflows.
"""


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


class tRNASeqWorkflow(WorkflowSuperClass):
    known_split_types = ['demethylase', 'untreated']

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='tRNAseq')

        self.target_files = [] # Snakemake target files
        self.sample_info = pd.DataFrame()
        self.samples_txt_file = None
        self.sample_names = None

        self.rules.extend([ # tRNAseq Snakemake rules
            'make_iu_input',
            'iu_gen_configs',
            'iu_merge_pairs',
            'gen_qc_report',
            'compress_merged_fasta',
            'anvi_reformat_fasta',
            'anvi_gen_tRNAseq_database',
            'compress_reformatted_fasta',
            'get_unique_tRNA_fasta',
            'prefix_derep_tRNA_fasta',
            'map_tRNA',
            'anvi_condense_clusters'])

        self.general_params.extend(['samples_txt']) # general section of config file

        # Parameters for each rule that are accessible in the config file
        rule_acceptable_params_dict = {}
        rule_acceptable_params_dict['iu_gen_configs'] = ['--r1-prefix', '--r2-prefix']
        rule_acceptable_params_dict['iu_merge_pairs'] = [
            '--marker-gene-stringent',
            '--max-num-mismatches',
            '--report-r1-prefix',
            '--report-r2-prefix']
        rule_acceptable_params_dict['compress_merged_fasta'] = ['run']
        rule_acceptable_params_dict['anvi_reformat_fasta'] = ['--simplify-names']
        rule_acceptable_params_dict['anvi_gen_tRNAseq_database'] = [
            '--charging-recorded', '--trust-fasta', '--verbose']
        rule_acceptable_params_dict['compress_reformatted_fasta'] = ['run']
        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Default values for certain accessible parameters
        self.default_config.update({
            'samples_txt': 'samples.txt',
            'iu_merge_pairs': {
                '--marker-gene-stringent': True,
                '--max-num-mismatches': 0,
                '--report-r1-prefix': False,
                '--report-r2-prefix': False,
                'threads': 1},
            'compress_merged_fasta': {'run': True},
            'anvi_reformat_fasta': {'--simplify-names': True},
            'anvi_gen_tRNAseq_database': {
                '--charging-recorded': False,
                '--trust-fasta': False,
                '--verbose': False,
                'threads': 1},
            'compress_reformatted_fasta': {'run': True}})

        # Directories where output beside log files are written
        self.dirs_dict.update({'QC_DIR': '01_QC', 'IDENT_DIR': '02_IDENT'})


    def init(self):
        super().init()

        # Load table of sample info from samples_txt (sample names, split types, paths to r1, r2)
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(self.samples_txt_file)
        try:
            # An error will subsequently be raised in `check_samples_txt` if there is no header
            self.sample_info = pd.read_csv(
                self.samples_txt_file, sep='\t', index_col=False)
        except IndexError as e:
            raise ConfigError(
                "The samples_txt file, '%s', does not appear to be properly formatted. "
                "This is the error from trying to load it: '%s'" % (self.samples_txt_file, e))

        self.check_samples_txt(self.samples_txt_file, self.sample_info)

        self.sample_names = self.sample_info['sample'].tolist()

        # Determine which samples have which splits
        self.sample_splits_dict = dict([(sample_name, []) for sample_name in self.sample_names])
        for sample_name, split_type in zip(self.sample_names, self.sample_info['split']):
            self.sample_splits_dict[sample_name].append(split_type)

        self.sample_split_prefixes = [
            sample_name + '_' + split_type
            for sample_name, split_type in zip(self.sample_names, self.sample_info['split'])]

        # The target files are extended depending on the rules set to run by the config file
        self.init_target_files(
            self.get_param_value_from_config(['compress_merged_fasta', 'run']),
            self.get_param_value_from_config(['compress_reformatted_fasta', 'run']))


    def init_target_files(self, run_compress_merged_fasta, run_compress_reformatted_fasta):

        # The target files of all tRNAseq workflows
        for sample_split_prefix in self.sample_split_prefixes:
            self.target_files.append(
                os.path.join(self.dirs_dict['IDENT_DIR'], sample_split_prefix + "_CONDENSED.bam"))

        # Compress by default unless specified otherwise in the config file
        if run_compress_merged_fasta:
            for sample_split_prefix in self.sample_split_prefixes:
                self.target_files.append(
                    os.path.join(self.dirs_dict['QC_DIR'], sample_split_prefix + "_MERGED.gz"))

        # Compress by default unless specified otherwise in the config file
        if run_compress_reformatted_fasta:
            for sample_split_prefix in self.sample_split_prefixes:
                self.target_files.append(
                    os.path.join(self.dirs_dict['QC_DIR'], sample_split_prefix + ".fasta.gz"))


    @staticmethod
    def check_samples_txt(sample_info_path, sample_info):

        # Must have four columns: sample, split, r1, r2
        missing_headers = []
        for header in ['sample', 'split', 'r1', 'r2']:
            if header not in sample_info.columns:
                missing_headers.append(header)
        if missing_headers:
            raise ConfigError(
                "The samples_txt file, '%s', is not properly formatted, "
                "as the following columns are missing: '%s'."
                % (sample_info, ', '.join(missing_headers)))

        for sample_name in sample_info['sample']:
            try:
                u.check_sample_id(sample_name)
            except ConfigError as e:
                raise ConfigError(
                    "While processing the samples_txt file, '%s', "
                    "Anvi'o ran into the following error: %s" % (samples_txt_file, e))

        unknown_split_types = []
        for split_type in sample_info['split']:
            if split_type not in tRNASeqWorkflow.known_split_types:
                unknown_split_types.append(split_type)
        if unknown_split_types:
            run.warning(
                "Some of the names of split types in the samples_txt file, '%s', "
                "are not what we were expecting (%s). "
                "That's okay, but Anvi'o decided it should warns you. "
                "Here are the names of split types that are not in our little list: %s. " % (
                    samples_txt_file,
                    ', '.join(tRNASeqWorkflow.known_split_types),
                    ', '.join(sorted(set(unknown_split_types)))))

        fastq_names = sample_info['r1'].tolist() + sample_info['r2'].tolist()
        bad_fastq_names = [
            s for s in fastq_names if (not s.endswith('.fastq') and not s.endswith('.fastq.gz'))]
        if bad_fastq_names:
            run.warning(
                "Some of the sequence files in the samples_txt file, '%s', "
                "do not end with either '.fastq' or '.fastq.gz'. "
                "That's okay, but Anvi'o decided it should warn you. "
                "Here are the first 5 such files that have unconventional file extensions: "
                "%s." % (samples_txt_file, ', '.join(bad_fastq_names[:5])))