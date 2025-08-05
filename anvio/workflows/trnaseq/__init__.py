# -*- coding: utf-8
# pylint: disable=line-too-long
"""Sets up an anvi'o tRNA-seq snakemake workflow."""

import os
import pandas as pd

from snakemake.io import ancient

import anvio
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.utils.validation import check_sample_id


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


run = terminal.Run()


class TRNASeqWorkflow(WorkflowSuperClass):
    """Sets up a workflow for the analysis of tRNA-seq reads."""

    KNOWN_TREATMENTS = ['untreated', 'demethylase']

    RECOGNIZED_PREFIX_NTS = constants.nucleotides

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='trnaseq')

        self.rules.extend([
            'iu_merge_pairs',
            'anvi_reformat_fasta',
            'anvi_trnaseq',
            'anvi_merge_trnaseq',
            'anvi_run_trna_taxonomy',
            'anvi_tabulate_trnaseq'
        ])

        # "General" section of the workflow config file.
        self.general_params.extend(['samples_txt'])

        # Parameters for each rule that are accessible in the config file.
        rule_acceptable_params_dict = {}
        rule_acceptable_params_dict['iu_merge_pairs'] = [
            'run',
            '--gzip-output',
            '--marker-gene-stringent',
            '--max-num-mismatches',
            '--report-r1-prefix',
            '--report-r2-prefix'
        ]
        rule_acceptable_params_dict['anvi_reformat_fasta'] = [
            'run',
            '--gzip-output',
            '--simplify-names'
        ]
        rule_acceptable_params_dict['anvi_trnaseq'] = [
            'run',
            '--treatment',
            '--overwrite-output-destinations',
            '--description',
            '--write-checkpoints',
            '--load-checkpoint',
            '--feature-param-file',
            '--threeprime-termini',
            '--min-length-long-fiveprime',
            '--min-trna-fragment-size',
            '--agglomeration-max-mismatch-freq',
            '--skip-INDEL-profiling',
            '--max-indel-freq',
            '--left-indel-buffer',
            '--right-indel-buffer',
            '--skip-fasta-check',
            '--alignment-target-chunk-size',
            '--profiling-chunk-size'
        ]
        rule_acceptable_params_dict['anvi_merge_trnaseq'] = [
            'run',
            '--project-name',
            '--max-reported-trna-seeds',
            '--overwrite-output-destinations',
            '--description',
            '--feature-threshold',
            '--preferred-treatment',
            '--nonspecific-output',
            '--min-variation',
            '--min-third-fourth-nt',
            '--min-indel-fraction',
            '--distance',
            '--linkage'
        ]
        rule_acceptable_params_dict['anvi_run_trna_taxonomy'] = [
            'run',
            '--trna-taxonomy-data-dir',
            '--min-percent-identity',
            '--max-num-target-sequences',
            '--num-parallel-processes',
            '--write-buffer-size'
        ]
        rule_acceptable_params_dict['anvi_tabulate_trnaseq'] = [
            'run',
            '--overwrite-output-destinations'
        ]
        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        # Default values for accessible parameters: all defaults are written to the config file so
        # the user can see them succinctly.

        # Though the workflow superclass automatically adds a threads argument of "" to each
        # workflow, here we make explicit that the default is 1 and the user does not need to
        # enclose the value in quotes. Likewise, the superclass adds mandatory arguments at the end
        # of the list for each rule in the config file, but we explicitly add them here to ensure
        # they appear in the order of each script's help display.
        self.default_config.update({
            'samples_txt': 'samples.txt',
            'iu_merge_pairs': {
                'run': True,
                '--gzip-output': False,
                '--marker-gene-stringent': True,
                '--max-num-mismatches': 0,
                '--report-r1-prefix': False,
                '--report-r2-prefix': False,
                'threads': 1
            },
            'anvi_reformat_fasta': {
                'run': True,
                '--gzip-output': False, # not an argument of anvi-script-reformat-fasta
                '--simplify-names': True, # not the default in anvi-script-reformat-fasta
                'threads': 1
            },
            'anvi_trnaseq': {
                'run': True,
                '--treatment': "", # if provided in the config file, the treatment is assumed to be for all samples
                '--overwrite-output-destinations': anvio.D['overwrite-output-destinations'][1]['default'],
                '--description': "",
                '--write-checkpoints': anvio.D['write-checkpoints'][1]['default'],
                '--load-checkpoint': "",
                '--feature-param-file': "",
                '--threeprime-termini': anvio.D['threeprime-termini'][1]['default'],
                '--min-length-long-fiveprime': anvio.D['min-length-long-fiveprime'][1]['default'],
                '--min-trna-fragment-size': anvio.D['min-trna-fragment-size'][1]['default'],
                '--agglomeration-max-mismatch-freq': anvio.D['agglomeration-max-mismatch-freq'][1]['default'],
                '--skip-INDEL-profiling': anvio.D['skip-INDEL-profiling'][1]['default'],
                '--max-indel-freq': anvio.D['max-indel-freq'][1]['default'],
                '--left-indel-buffer': anvio.D['left-indel-buffer'][1]['default'],
                '--right-indel-buffer': anvio.D['right-indel-buffer'][1]['default'],
                '--skip-fasta-check': True, # not the default in anvi-trnaseq
                '--profiling-chunk-size': anvio.D['profiling-chunk-size'][1]['default'],
                '--alignment-target-chunk-size': anvio.D['alignment-target-chunk-size'][1]['default'],
                'threads': 1
            },
            'anvi_merge_trnaseq': {
                'run': True,
                '--project-name': "",
                '--max-reported-trna-seeds': anvio.D['max-reported-trna-seeds'][1]['default'],
                '--overwrite-output-destinations': anvio.D['overwrite-output-destinations'][1]['default'],
                '--description': "",
                '--feature-threshold': anvio.D['feature-threshold'][1]['default'],
                '--preferred-treatment': "",
                '--nonspecific-output': anvio.D['nonspecific-output'][1]['default'],
                '--min-variation': anvio.D['min-variation'][1]['default'],
                '--min-third-fourth-nt': anvio.D['min-third-fourth-nt'][1]['default'],
                '--min-indel-fraction': anvio.D['min-indel-fraction'][1]['default'],
                '--distance': anvio.D['distance'][1]['default'],
                '--linkage': anvio.D['linkage'][1]['default'],
                'threads': 1
            },
            'anvi_run_trna_taxonomy': {
                'run': True,
                '--trna-taxonomy-data-dir': "",
                '--min-percent-identity': 90, # default in anvi-run-trna-taxonomy
                '--max-num-target-sequences': 100, # default in anvi-run-trna-taxonomy
                '--num-parallel-processes': anvio.D['num-parallel-processes'][1]['default'],
                '--write-buffer-size': anvio.D['write-buffer-size'][1]['default'],
                'threads': 1
            },
            'anvi_tabulate_trnaseq': {
                'run': True,
                '--overwrite-output-destinations': anvio.D['overwrite-output-destinations'][1]['default'],
                'threads': 1
            },
            'output_dirs': {}, # This ensures that output_dirs comes before max_threads in the file
            'max_threads': 1
        })

        self.dirs_dict.update({
            'QC_DIR': '01_QC',
            'IDENT_DIR': '02_IDENT',
            'CONVERT_DIR': '03_CONVERT'})


    def init(self):
        """This function is called from within the snakefile to initialize parameters."""
        super().init()

        self.run_iu_merge_pairs = self.get_param_value_from_config(['iu_merge_pairs', 'run'])
        self.gzip_iu_merge_pairs_output = self.get_param_value_from_config(['iu_merge_pairs', '--gzip-output']) if self.run_iu_merge_pairs else False
        self.run_anvi_reformat_fasta = self.get_param_value_from_config(['anvi_reformat_fasta', 'run'])
        self.gzip_anvi_reformat_output = self.get_param_value_from_config(['anvi_reformat_fasta', '--gzip-output']) if self.run_anvi_reformat_fasta else False
        self.run_anvi_trnaseq = self.get_param_value_from_config(['anvi_trnaseq', 'run'])
        self.run_anvi_merge_trnaseq = self.get_param_value_from_config(['anvi_merge_trnaseq', 'run'])
        self.run_anvi_run_trna_taxonomy = self.get_param_value_from_config(['anvi_run_trna_taxonomy', 'run'])
        self.run_anvi_tabulate_trnaseq = self.get_param_value_from_config(['anvi_tabulate_trnaseq', 'run'])

        # Load table of sample info from samples_txt.
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(self.samples_txt_file)
        try:
            # An error will subsequently be raised in `check_samples_txt` if there is no header.
            self.sample_info = pd.read_csv(self.samples_txt_file, sep='\t', index_col=False)
        except IndexError as e:
            raise ConfigError("The samples_txt file, '%s', does not appear to be properly formatted. "
                              "This is the error from trying to load it: '%s'" % (self.samples_txt_file, e))
        self.check_samples_txt()

        self.sample_names = self.sample_info['sample'].tolist()
        if 'treatment' not in self.sample_info:
            # The treatment is the same for each sample and is set in the config file.
            self.sample_info['treatment'] = [self.get_param_value_from_config(['anvi_trnaseq', '--treatment'])] * len(self.sample_names)
        if self.run_iu_merge_pairs:
            self.treatments = self.sample_info['treatment']
            self.r1_paths = self.sample_info['r1'].tolist()
            self.r2_paths = self.sample_info['r2'].tolist()
            self.r1_prefixes = self.get_r1_prefixes()
            self.r2_prefixes = self.get_r2_prefixes()
            self.fasta_paths = None
        else:
            self.treatments = self.sample_info['treatment']
            self.r1_paths = None
            self.r2_paths = None
            self.r1_prefixes = None
            self.r2_prefixes = None
            self.fasta_paths = self.sample_info['fasta'].tolist()

        self.target_files = self.get_target_files()

        # The `anvi-run-workflow --cluster` option, which submits each rule as a separate job,
        # requires that the rule's log directory exist before running the rule. This workflow
        # differs from others by writing log files for each sample to log directories for each
        # sample.
        for sample_name in self.sample_names:
            filesnpaths.gen_output_directory(os.path.join(self.dirs_dict['LOGS_DIR'], sample_name))


    def check_samples_txt(self):
        """Check the format and content of the file of sample information provided by the user."""
        self.check_samples_txt_header()
        self.check_sample_names()
        self.check_treatments()
        self.check_input_paths()
        self.check_project_name()


    def check_samples_txt_header(self):
        if self.run_iu_merge_pairs:
            min_header = ['sample', 'r1', 'r2']
        else:
            min_header = ['sample', 'fasta']
        missing_columns = []
        for column_title in min_header:
            if column_title not in self.sample_info.columns:
                missing_columns.append(column_title)
        if missing_columns:
            raise ConfigError("The samples_txt file, '%s', is not properly formatted, "
                              "as the following columns are missing: '%s'."
                              % (self.samples_txt_file, ', '.join(missing_columns)))


    def check_sample_names(self):
        """Check that the name of each tRNA-seq library is anvi'o-compliant."""
        for sample_name in self.sample_info['sample']:
            try:
                check_sample_id(sample_name)
            except ConfigError as e:
                raise ConfigError("While processing the samples_txt file, '%s', "
                                  "anvi'o ran into the following error: %s" % (self.samples_txt_file, e))

        if len(set(self.sample_info['sample'])) != len(self.sample_info['sample']):
            raise ConfigError("Sample names in the samples_txt file, '%s', must be unique." % self.samples_txt_file)


    def check_treatments(self):
        """Check the types of sample treatments, which are primarily meant to differentiate
        technical replicates involving different enzymatic treatments to remove nucleotide
        modifications."""
        if 'treatment' in self.sample_info.columns:
            treatments = self.sample_info['treatment'].tolist()
            unknown_treatments = []
            for treatment in treatments:
                if treatment not in TRNASeqWorkflow.KNOWN_TREATMENTS:
                    unknown_treatments.append(treatment)
            if unknown_treatments:
                run.warning("Some of the names of treatments in the samples_txt file, '%s', "
                            "are not known to anvi'o (these are the known treatments: %s). "
                            "That's okay, but anvi'o decided it should warn you. "
                            "Here are the names of treatments that are not in our little list: %s. "
                            % (self.samples_txt_file,
                               ', '.join(TRNASeqWorkflow.KNOWN_TREATMENTS),
                               ', '.join(sorted(set(unknown_treatments)))))
        else:
            treatment = self.get_param_value_from_config(['anvi_trnaseq', '--treatment'])
            if treatment:
                if treatment not in TRNASeqWorkflow.KNOWN_TREATMENTS:
                    run.warning("The treatment, %s, specified in the config file "
                                "is not known to anvi'o (these are the known treatments: %s). "
                                "That's okay, but anvi'o decided it should warn you."
                                % (treatment,
                                   ', '.join(TRNASeqWorkflow.KNOWN_TREATMENTS)))
            else:
                run.warning("Since you did not specify treatments for your samples -- "
                            "either individually for each sample in the samples_txt file, '%s', "
                            "or for all samples in the config file using the anvi-trnaseq parameter, `--treatment` -- "
                            "anvi'o will set the treatment for all samples to \"untreated\", "
                            "assuming that no enzymatic treatment was applied." % self.samples_txt_file)


    def check_input_paths(self):
        """Check FASTQ file paths if running Illumina-utils for read merging, or FASTA file paths if
        considering merged or unpaired reads. Allow both absolute and relative paths in
        samples_txt."""
        if self.run_iu_merge_pairs:
            fastq_paths = self.sample_info['r1'].tolist() + self.sample_info['r2'].tolist()
            bad_fastq_paths = []
            for fastq_path in fastq_paths:
                if os.path.isabs(fastq_path):
                    if not filesnpaths.is_file_exists(fastq_path, dont_raise=True):
                        bad_fastq_paths.append(fastq_path)
                else:
                    if not filesnpaths.is_file_exists(os.path.join(os.getcwd(), fastq_path), dont_raise=True):
                        bad_fastq_paths.append(fastq_path)
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
                            "That's okay, but anvi'o decided it should warn you. "
                            "Here are the first 5 such files that have unconventional file extensions: %s."
                            % (self.samples_txt_file, ', '.join(bad_fastq_names[:5])))
        else:
            fasta_paths = self.sample_info['fasta'].tolist()
            bad_fasta_paths = []
            for fasta_path in fasta_paths:
                if os.path.isabs(fasta_path):
                    if not filesnpaths.is_file_exists(fasta_path, dont_raise=True):
                        bad_fasta_paths.append(fasta_path)
                else:
                    if not filesnpaths.is_file_exists(os.path.join(os.getcwd(), fasta_path), dont_raise=True):
                        bad_fasta_paths.append(fasta_path)
            bad_fasta_paths = [s for s in fasta_paths if not filesnpaths.is_file_exists(os.path.abspath(s), dont_raise=True)]
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
                            "That's okay, but anvi'o decided it should warn you. "
                            "Here are the first 5 such files that have unconventional file extensions: %s."
                            % (self.samples_txt_file, ', '.join(bad_fasta_names[:5])))


    def check_project_name(self):
        """Check the name of the tRNA-seq project."""
        if self.run_anvi_merge_trnaseq:
            project_name = self.get_param_value_from_config(['anvi_merge_trnaseq', '--project-name'])
            if not project_name:
                raise ConfigError("Since you are running anvi-merge-trnaseq, "
                                  "please provide a project name for the sample(s) in the config file.")
            try:
                check_sample_id(project_name)
            except ConfigError as e:
                raise ConfigError("While checking the project name, '%s', "
                                  "anvi'o ran into the following error: %s" % (project_name, e))


    def get_r1_prefixes(self):
        """Subsequences found at the beginning of reads can include sample barcodes and unimolecular
        identifiers (random nucleotides). The user should provide these sequences in samples_txt.
        This function converts read 1 prefixes into regex strings that are recognized by
        Illumina-utils."""
        if 'r1_prefix' not in self.sample_info:
            return None

        r1_prefixes = []
        for r1_prefix_input in self.sample_info['r1_prefix'].tolist():
            r1_prefix = '^'
            for nt in r1_prefix_input:
                if nt not in self.RECOGNIZED_PREFIX_NTS:
                    raise ConfigError("The read prefix sequences provided in the samples_txt file, '%s', "
                                      "do not exclusively contain recognized nucleotide characters (%s). "
                                      "One faulty read 1 prefix that was provided is %s."
                                      % (self.samples_txt_file, ', '.join(self.RECOGNIZED_PREFIX_NTS), r1_prefix_input))
                if nt == 'N':
                    r1_prefix += '.'
                else:
                    r1_prefix += nt
            if r1_prefix == '^':
                r1_prefix == ''
            r1_prefixes.append(r1_prefix)
        return r1_prefixes


    def get_r2_prefixes(self):
        """Subsequences found at the beginning of reads can include sample barcodes and unimolecular
        identifiers (random nucleotides). The user should provide these sequences in samples_txt.
        This function converts read 2 prefixes into regex strings that are recognized by
        Illumina-utils."""
        if 'r2_prefix' not in self.sample_info:
            return None

        r2_prefixes = []
        for r2_prefix_input in self.sample_info['r2_prefix'].tolist():
            r2_prefix = '^'
            for nt in r2_prefix_input:
                if nt not in self.RECOGNIZED_PREFIX_NTS:
                    raise ConfigError("The read prefix sequences provided in the samples_txt file, '%s', "
                                      "do not exclusively contain recognized nucleotide characters (%s). "
                                      "One faulty read 2 prefix that was provided is %s."
                                      % (self.samples_txt_file, ', '.join(self.RECOGNIZED_PREFIX_NTS), r2_prefix_input))
                if nt == 'N':
                    r2_prefix += '.'
                else:
                    r2_prefix += nt
            if r2_prefix == '^':
                r2_prefix = ''
            r2_prefixes.append(r2_prefix)
        return r2_prefixes


    def get_target_files(self):
        """Determine the snakemake target files."""
        target_files = []

        if self.run_iu_merge_pairs:
            # Each library of paired-end reads is run separately through Illumina-utils, as each can
            # have different prefix sequences, and a single run cannot handle multiple prefixes.
            target_files.append(os.path.join(self.dirs_dict['QC_DIR'], "qc_report.txt"))

        if self.run_anvi_reformat_fasta:
            # Do not use `REFORMAT.done` output files produced by rule `run_anvi-reformat_fasta` as
            # workflow targets. If this were the case and these files were last modified before
            # `run_anvi_reformat_fasta` input files (`MERGE.done` when starting with FASTQ files) --
            # which can occur when copying the contents of `01_QC` from one location to another --
            # then `run_anvi_reformat_fasta` would always be rerun, even if attempting to start the
            # workflow after this rule, say at `anvi_trnaseq`. Therefore, touch
            # `ALL_REFORMATTING.done` as the output of dummy rule `all_reformatting_done` in lieu of
            # using `REFORMAT.done` files as workflow targets.
            target_files.append(os.path.join(self.dirs_dict['QC_DIR'], "ALL_REFORMATTING.done"))

        if self.run_anvi_trnaseq:
            for sample_name in self.sample_names:
                out_dir = os.path.join(self.dirs_dict['IDENT_DIR'], sample_name)
                # If anvi-trnaseq output files such as the tRNA-seq database or analysis summary
                # were used as targets, they would be deleted upon workflow failure.
                target_files.append(os.path.join(out_dir, "IDENT.done"))

        if self.run_anvi_merge_trnaseq:
            project_name = self.get_param_value_from_config(['anvi_merge_trnaseq', '--project-name'])
            target_files.append(os.path.join(self.dirs_dict['CONVERT_DIR'], "CONVERT.done"))

        if self.run_anvi_run_trna_taxonomy:
            target_files.append(os.path.join(self.dirs_dict['CONVERT_DIR'], "TAXONOMY.done"))

        if self.run_anvi_tabulate_trnaseq:
            target_files.append(os.path.join(self.dirs_dict['CONVERT_DIR'], "TABULATE.done"))

        return target_files


    def get_input_for_anvi_reformat_fasta(self, wildcards):
        """Input can come from two possible sources: a user-supplied FASTA file or the FASTA file of
        merged reads generated from user-supplied FASTQ files."""
        sample_name = wildcards.sample_name
        if self.fasta_paths:
            return ancient(self.fasta_paths[self.sample_names.index(sample_name)])
        return ancient(os.path.join(os.path.join(self.dirs_dict['QC_DIR'], sample_name), "MERGE.done"))


    def get_input_for_anvi_trnaseq(self, wildcards):
        """Input can come from two possible sources: a FASTA file with Anvi'o-compliant deflines
        supplied by the user or the reformatted FASTA file produced by the rule,
        anvi_reformat_fasta."""
        sample_name = wildcards.sample_name
        if self.run_anvi_reformat_fasta:
            return os.path.join(os.path.join(self.dirs_dict['QC_DIR'], sample_name), "REFORMAT.done")
        elif self.fasta_paths:
            return self.fasta_paths[self.sample_names.index(sample_name)]
        else:
            raise ConfigError("The rule, `anvi_reformat_fasta`, should be run after the rule, `iu_merge_pairs`, "
                              "to produce an anvi'o-compliant FASTA file. "
                              "The `run` parameter for `anvi_reformat_fasta` should be set to `True` in the workflow config file.")
