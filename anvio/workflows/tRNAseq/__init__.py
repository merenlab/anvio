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
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='tRNAseq')

        self.target_files = [] # TODO: Once we update all other workflows then this will be initiated in WorkflowSuperClass
        self.samples_information = pd.DataFrame()
        self.samples_txt_file = None
        self.sample_names = None
        self.valid_splits = ['demethylase', 'untreated']

        self.rules.extend(['make_iu_input',
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

        self.general_params.extend(['samples_txt'])

        rule_acceptable_params_dict = {}

        # defining the accessible params per rule
        rule_acceptable_params_dict['iu_gen_configs'] = ['--r1-prefix',
                                                         '--r2-prefix']
        rule_acceptable_params_dict['iu_merge_pairs'] = ['--marker-gene-stringent',
                                                         '--max-num-mismatches',
                                                         '--report-r1-prefix',
                                                         '--report-r2-prefix']
        rule_acceptable_params_dict['compress_merged_fasta'] = ['run']
        rule_acceptable_params_dict['anvi_reformat_fasta'] = ['--simplify-names']
        rule_acceptable_params_dict['anvi_gen_tRNAseq_database'] = ['--charging-recorded',
                                                                    '--trust-fasta',
                                                                    '--verbose']
        rule_acceptable_params_dict['compress_reformatted_fasta'] = ['run']

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.dirs_dict.update({'QC_DIR': '01_QC',
                               'IDENT_DIR': '02_IDENT'})

        self.default_config.update({'samples_txt': 'samples.txt',
                                    'iu_merge_pairs': {'--marker-gene-stringent': True,
                                                       '--max-num-mismatches': 0,
                                                       '--report-r1-prefix': False,
                                                       '--report-r2-prefix': False,
                                                       'threads': 1},
                                    'compress_merged_fasta': {'run': True},
                                    'anvi_reformat_fasta': {'--simplify-names': True},
                                    'anvi_gen_tRNAseq_database': {'--charging-recorded': False,
                                                                  '--trust-fasta': False,
                                                                  '--verbose': False,
                                                                  'threads': 1},
                                    'compress_reformatted_fasta': {'run': True}})


    def init(self):
        super().init()

        # loading the samples.txt file
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(self.samples_txt_file)
        try:
            # getting the samples information
            # (names, splits, [groups], paths to r1, paths to r2)
            # from samples.txt
            self.samples_information = pd.read_csv(self.samples_txt_file, sep='\t', index_col=False)
        except IndexError as e:
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. \
                               This is what we know: '%s'" % (self.samples_txt_file, e))
        if 'sample' not in self.samples_information.columns:
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. \
                               We aren't seeing a column with title 'sample'." % self.samples_txt_file)
        if 'split' not in self.samples_information.columns:
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. \
                               We aren't seeing a column with title 'split'." % self.samples_txt_file)
        if 'r1' not in self.samples_information.columns:
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. \
                               We aren't seeing a column with title 'r1'." % self.samples_txt_file)
        if 'r2' not in self.samples_information.columns:
            raise ConfigError("Looks like your samples_txt file, '%s', is not properly formatted. \
                               We aren't seeing a column with title 'r2'." % self.samples_txt_file)

        self.check_samples_txt()
        # get a list of the sample names
        self.sample_names = self.samples_information['sample'].tolist()[::len(self.valid_splits)]

        self.init_target_files(self.get_param_value_from_config(['compress_merged_fasta', 'run']),
                               self.get_param_value_from_config(['compress_reformatted_fasta', 'run']))


    def check_samples_txt(self):
        sample_names = self.samples_information['sample']
        num_valid_splits = len(self.valid_splits)

        for i in range(0, len(sample_names), num_valid_splits):
            if len(set(sample_names[i: i + num_valid_splits])) != 1:
                raise ConfigError("There should be %d consecutive lines "
                                  "for each sample in your samples_txt file ('%s'). "
                                  "Each line should be for a different split. "
                                  "This is not the case, for example, with this sample: "
                                  "%s" % (num_valid_splits, self.samples_txt_file, sample_names[i]))
        if len(sample_names) / num_valid_splits != len(set(sample_names)):
            raise ConfigError("Names of samples in your samples_txt file, %s, must be unique. "
                              "Each sample should be represented by %d consecutive lines "
                              "for different splits." % (self.samples_txt_file, len(self.valid_splits)))
        for sample_name in sample_names:
            try:
                u.check_sample_id(sample_name)
            except ConfigError as e:
                raise ConfigError("While processing the samples_txt file ('%s'), "
                                  "Anvi'o ran into the following error: "
                                  "%s" % (self.samples_txt_file, e))

        split_entries = self.samples_information['split']
        valid_splits_set = set(self.valid_splits)
        for i in range(0, len(split_entries), num_valid_splits):
            if set(split_entries[i: i + num_valid_splits]) != valid_splits_set:
                raise ConfigError("Each sample in your samples_txt file ('%s') "
                                  "should have splits corresponding to these valid splits: %s. "
                                  "This is not the case, for example, with this sample: "
                                  "%s" % (self.samples_txt_file, ', '.join(self.valid_splits), sample_names[i]))

        fastq_file_names = self.samples_information['r1'].tolist() + self.samples_information['r2'].tolist()
        bad_fastq_names = [s for s in fastq_file_names if (not s.endswith('.fastq') and not s.endswith('.fastq.gz'))]
        if bad_fastq_names:
            run.warning("We noticed some of your sequence files in your samples_txt file ('%s') "
                        "do not end with either '.fastq' or '.fastq.gz'. "
                        "That's okay, but anvi'o decided it should warn you. "
                        "Here are the first 5 such files that have unconventional file extensions: "
                        "%s." % (self.samples_txt_file, ', '.join(bad_fastq_names[:5])))


    def init_target_files(self, run_compress_merged_fasta, run_compress_reformatted_fasta):

        for sample in self.sample_names:
            for split in self.valid_splits:
                self.target_files.append(os.path.join(self.dirs_dict['IDENT_DIR'], "%s_%s_CONDENSED.bam" % (sample, split)))

        if run_compress_merged_fasta:
            for sample in self.sample_names:
                for split in self.valid_splits:
                    self.target_files.append(os.path.join(self.dirs_dict['QC_DIR'], "%s_%s_MERGED.gz" % (sample, split)))

        if run_compress_reformatted_fasta:
            for sample in self.sample_names:
                for split in self.valid_splits:
                    self.target_files.append(os.path.join(self.dirs_dict['QC_DIR'], "%s_%s.fasta.gz" % (sample, split)))