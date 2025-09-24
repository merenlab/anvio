# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import os
import anvio
import pandas as pd
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio import utils as u
from anvio.drivers import driver_modules
from anvio.workflows import WorkflowSuperClass
from anvio.workflows.contigs import ContigsDBWorkflow
from anvio.errors import ConfigError
from anvio.artifacts.samples_txt import SamplesTxt


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()

min_contig_length_for_assembly = 1000

class MetagenomicsWorkflow(ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='metagenomics')

        self.samples_information = {}
        self.kraken_annotation_dict = {}
        self.run_krakenuniq = None
        self.run_metaspades = None
        self.use_scaffold_from_metaspades = None
        self.use_scaffold_from_idba_ud = None
        self.remove_short_reads_based_on_references = None
        self.references_for_removal_txt = None
        self.references_for_removal = {}
        self.references_mode = None
        self.fasta_txt_file = None
        self.samples_txt_file = None
        self.sample_names = None
        self.group_sizes = None
        self.collections_txt = None
        self.collections = None

        # initialize the base class
        ContigsDBWorkflow.__init__(self)

        self.rules.extend(['iu_gen_configs', 'iu_filter_quality_minoche', 'gen_qc_report', 'gzip_fastqs',\
                     'merge_fastqs_for_co_assembly', 'megahit', 'merge_fastas_for_co_assembly',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam', 'idba_ud',\
                     'anvi_profile', 'anvi_merge', 'import_percent_of_reads_mapped', 'anvi_cluster_contigs',\
                     'krakenuniq', 'krakenuniq_mpa_report', 'import_krakenuniq_taxonomy', 'metaspades',\
                     'flye', 'hifiasm_meta', 'minimap2_index', 'minimap2',\
                     'remove_short_reads_based_on_references', 'anvi_summarize', 'anvi_split'])

        self.general_params.extend(['samples_txt', "references_mode", "all_against_all",\
                                    "kraken_txt", "collections_txt", "read_type_suffix"])

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule. NOTE --threads is a parameter for every rule
        # and is not explicitly provided in what follows
        rule_acceptable_params_dict['iu_gen_configs'] = ["--r1-prefix", "--r2-prefix"]
        rule_acceptable_params_dict['iu_filter_quality_minoche'] = ['run', '--visualize-quality-curves', '--ignore-deflines', '--limit-num-pairs', '--print-qual-scores', '--store-read-fate']
        rule_acceptable_params_dict['gzip_fastqs'] = ["run"]

        # add parameters for modifying binning algorithms
        additional_params_for_anvi_cluster_contigs = [self.get_param_name_for_binning_driver(d) for d in driver_modules['binning'].keys()]
        rule_acceptable_params_dict['anvi_cluster_contigs'] = ["run", "--collection-name", "--driver", "--just-do-it"]
        rule_acceptable_params_dict['anvi_cluster_contigs'].extend(additional_params_for_anvi_cluster_contigs)

        rule_acceptable_params_dict['anvi_summarize'] = ["additional_params", "run"]
        rule_acceptable_params_dict['anvi_split'] = ["additional_params", "run"]
        rule_acceptable_params_dict['metaspades'] = ["run", "conda_yaml", "conda_env", "additional_params", "use_scaffolds"]
        rule_acceptable_params_dict['megahit'] = ["run", "conda_yaml", "conda_env", "--min-contig-len", "--min-count", "--k-min",
                                                  "--k-max", "--k-step", "--k-list",
                                                  "--no-mercy", "--no-bubble", "--merge-level",
                                                  "--prune-level", "--prune-depth", "--low-local-ratio",
                                                  "--max-tip-len", "--no-local", "--kmin-1pass",
                                                  "--presets", "--memory", "--mem-flag",
                                                  "--use-gpu", "--gpu-mem", "--keep-tmp-files",
                                                  "--tmp-dir", "--continue", "--verbose"]
        rule_acceptable_params_dict['idba_ud'] = ["run", "conda_yaml", "conda_env", "--mink", "--maxk", "--step", "--inner_mink",
                                                  "--inner_step", "--prefix", "--min_count",
                                                  "--min_support", "--seed_kmer", "--min_contig",
                                                  "--similar", "--max_mismatch", "--min_pairs",
                                                  "--no_bubble", "--no_local", "--no_coverage",
                                                  "--no_correct", "--pre_correction", "use_scaffolds"]
        rule_acceptable_params_dict['flye'] = ["run", "conda_yaml", "conda_env", "--meta", "--pacbio-raw", "--pacbio-corr",
                                                   "--pacbio-hifi", "--nano-raw", "--nano-corr",
                                                   "--nano-hq", "--genome-size", "--iterations",
                                                   "--min-overlap", "--read-error", "--keep-haplotypes",
                                                   "--no-alt-contigs", "--scaffold", "--polish-target",
                                                   "additional_params", "threads"]
        rule_acceptable_params_dict['hifiasm_meta'] = ["run", "conda_yaml", "conda_env", "additional_params", "threads"]
        rule_acceptable_params_dict['bowtie'] = ["conda_yaml", "conda_env", "additional_params"]
        rule_acceptable_params_dict['bowtie_build'] = ["conda_yaml", "conda_env", "additional_params"]
        rule_acceptable_params_dict['minimap2_index'] = ["conda_yaml", "conda_env", "additional_params"]
        rule_acceptable_params_dict['minimap2'] = ["preset","conda_yaml", "conda_env",  "additional_params", "threads"]
        rule_acceptable_params_dict['samtools_view'] = ["additional_params"]
        rule_acceptable_params_dict['anvi_profile'] = ["--overwrite-output-destinations", "--report-variability-full",
                                                        "--skip-SNV-profiling", "--profile-SCVs", "--description",
                                                        "--skip-hierarchical-clustering", "--distance", "--linkage", "--min-contig-length",
                                                        "--min-mean-coverage", "--min-coverage-for-variability", "--cluster-contigs",
                                                        "--contigs-of-interest", "--queue-size", "--write-buffer-size-per-thread",
                                                        "--fetch-filter", "--min-percent-identity", "--max-contig-length"]
        rule_acceptable_params_dict['merge_fastas_for_co_assembly'] = []
        rule_acceptable_params_dict['merge_fastqs_for_co_assembly'] = []
        rule_acceptable_params_dict['anvi_merge'] = ["--sample-name", "--description", "--skip-hierarchical-clustering",
                                                     "--enforce-hierarchical-clustering", "--distance", "--linkage",
                                                     "--overwrite-output-destinations"]
        rule_acceptable_params_dict['import_percent_of_reads_mapped'] = ["run"]
        rule_acceptable_params_dict['krakenuniq'] = ["additional_params", "run", "--db", "--gzip-compressed"]
        rule_acceptable_params_dict['import_krakenuniq_taxonomy'] = ["--min-abundance"]
        rule_acceptable_params_dict['remove_short_reads_based_on_references'] = ["dont_remove_just_map", \
                                                                                 "references_for_removal_txt", \
                                                                                 "delimiter-for-iu-remove-ids-from-fastq"]

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        forbidden_params = {}
        forbidden_params['krakenuniq'] = ['--fastq-input', '--paired', '--output']

        self.forbidden_params.update(forbidden_params)

        self.dirs_dict.update({"QC_DIR": "01_QC",
                               "FASTA_DIR": "02_FASTA",
                               "CONTIGS_DIR": "03_CONTIGS",
                               "MAPPING_DIR": "04_MAPPING",
                               "PROFILE_DIR": "05_ANVIO_PROFILE",
                               "MERGE_DIR": "06_MERGED",
                               "TAXONOMY_DIR": "07_TAXONOMY",
                               "SUMMARY_DIR": "08_SUMMARY",
                               "SPLIT_PROFILES_DIR": "09_SPLIT_PROFILES"})

        self.default_config.update({'samples_txt': "samples.txt",
                                    'read_type_suffix': 'auto',
                                    'metaspades': {"additional_params": "--only-assembler", "threads": 7},
                                    'megahit': {"--min-contig-len": min_contig_length_for_assembly, "--memory": 0.4, "threads": 7},
                                    'idba_ud': {"--min_contig": min_contig_length_for_assembly, "threads": 7},
                                    'flye': {"run": False, "threads": 7, "additional_params": ""},
                                    'hifiasm_meta': {"run": False, "threads": 7, "additional_params": ""},
                                    'iu_filter_quality_minoche': {"run": True, "--ignore-deflines": True},
                                    "gzip_fastqs": {"run": True},
                                    "bowtie": {"additional_params": "--no-unal", "threads": 3},
                                    'minimap2_index': {"additional_params": ""},
                                    'minimap2': {"threads": 3, "preset": "map-hifi", "additional_params": "--secondary-seq"},
                                    "samtools_view": {"additional_params": "-F 4"},
                                    "anvi_profile": {"threads": 3, "--overwrite-output-destinations": True},
                                    "anvi_merge": {"--sample-name": "{group}", "--overwrite-output-destinations": True},
                                    "import_percent_of_reads_mapped": {"run": True},
                                    "krakenuniq": {"threads": 3, "--gzip-compressed": True, "additional_params": ""},
                                    "remove_short_reads_based_on_references": {"delimiter-for-iu-remove-ids-from-fastq": " "},
                                    "anvi_cluster_contigs": {"--collection-name": "{driver}"}})


    def init(self):
        super().init()

        # loading the samples.txt file
        samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        self.samples_txt = SamplesTxt(samples_txt_file, expected_format="free")
        self.samples_information = self.samples_txt.as_df()
        # get a list of the sample names
        self.sample_names = self.samples_txt.samples()
        self.run_metaspades = self.get_param_value_from_config(['metaspades', 'run'])

        # read type suffix policy for readsets: 'auto' (default) or 'force'
        read_type_suffix = self.get_param_value_from_config('read_type_suffix') or 'auto'

        # get the readsets, a.k.a. the fundamental unit that will be assembled and/or mapped.
        self.readsets = self.samples_txt.get_readsets(read_type_suffix=read_type_suffix)
        self.readset_ids = [rs['id'] for rs in self.readsets]

        # quick flags about data composition
        self.has_sr = self.samples_txt.has_any_sr()
        self.has_lr = self.samples_txt.has_any_lr()

        # sanity checks: assemblers required when assembling (not in references mode)
        if not self.references_mode:
            # SR present → require exactly one SR assembler enabled
            sr_choices = [bool(self.get_param_value_from_config(['megahit', 'run'])),
                          bool(self.get_param_value_from_config(['metaspades', 'run'])),
                          bool(self.get_param_value_from_config(['idba_ud', 'run'])),]
            if self.has_sr and sum(sr_choices) == 0:
                raise ConfigError("Short-reads detected in samples.txt, but no short-read assembler is enabled "
                                  "(expected one of: megahit, metaspades, idba_ud).")
            if sum(sr_choices) > 1:
                raise ConfigError("Multiple short-read assemblers are enabled; please enable only one of "
                                  "megahit, metaspades, idba_ud.")

            # LR present → require exactly one LR assembler enabled
            lr_choices = [bool(self.get_param_value_from_config(['flye', 'run'])),
                          bool(self.get_param_value_from_config(['hifiasm_meta', 'run'])),]
            if self.has_lr and sum(lr_choices) == 0:
                raise ConfigError("Long-reads detected in samples.txt, but no long-read assembler is enabled "
                                  "(expected one of: flye, hifiasm_meta).")
            if sum(lr_choices) > 1:
                raise ConfigError("Multiple long-read assemblers are enabled; please enable only one of "
                                  "flye, hifiasm_meta.")

        # sanity check for conda env: use either conda_yaml or conda_env, not both
        for tool in ['flye','hifiasm_meta','minimap2','minimap2_index','bowtie_build','bowtie2','megahit','metaspades','idba_ud']:
            y = self.get_param_value_from_config([tool, 'conda_yaml'])
            n = self.get_param_value_from_config([tool, 'conda_env'])
            if (y and y.strip()) and (n and n.strip()):
                raise ConfigError(f"For '{tool}', please set only one of 'conda_yaml' (YAML path) "
                                  f"or 'conda_env' (existing env name), not both.")


        self.use_scaffold_from_metaspades = self.get_param_value_from_config(['metaspades', 'use_scaffolds'])
        self.use_scaffold_from_idba_ud = self.get_param_value_from_config(['idba_ud', 'use_scaffolds'])
        self.run_qc = self.get_param_value_from_config(['iu_filter_quality_minoche', 'run']) == True
        self.run_summary = self.get_param_value_from_config(['anvi_summarize', 'run']) == True
        self.run_split = self.get_param_value_from_config(['anvi_split', 'run']) == True
        self.references_mode = self.get_param_value_from_config('references_mode')
        self.fasta_txt_file = self.get_param_value_from_config('fasta_txt')
        self.profile_databases = {}

        # just some extra checks. TO BE UTLRA-SAFE
        if len(self.sample_names) < 1:
            raise WorkflowError("No samples found in samples.txt")

        self.references_for_removal_txt = self.get_param_value_from_config(['remove_short_reads_based_on_references',\
                                                                            'references_for_removal_txt'])
        if self.references_for_removal_txt:
            self.load_references_for_removal()

        self.collections_txt = self.get_param_value_from_config('collections_txt')
        if self.collections_txt:
            self.load_collections()
        elif self.run_summary:
            raise ConfigError('If you want to run anvi-summarize you must provide a collections_txt file')
        elif self.run_split:
            raise ConfigError('If you want to run anvi-split you must provide a collections_txt file')

        self.init_kraken()
        self.init_references_txt()

        # Set the PROFILE databases paths variable:
        for group in self.group_names:
            if self.group_sizes[group] > 1:
                self.profile_databases[group] = os.path.join(self.dirs_dict["MERGE_DIR"], group, "PROFILE.db")
            else:
                if not self.references_mode:
                    # Use the concrete readset id of the single member (e.g., S1 or S3_SR)
                    member_ids = self.assembly_members.get(group, [])
                    if not member_ids:
                        raise ConfigError(f"Internal error: no members recorded for single-member group '{group}'.")
                    single_member = member_ids[0]
                else:
                    # References mode keeps legacy behavior: raw groups map has the sample name
                    single_member = self.samples_txt.groups()[group][0]
                self.profile_databases[group] = os.path.join(self.dirs_dict["PROFILE_DIR"], group, single_member, "PROFILE.db")


    def get_metagenomics_target_files(self):

        target_files = []

        target_files.extend(list(self.profile_databases.values()))

        # for groups of size 1 we create a message file
        message_file_for_groups_of_size_1 = [os.path.join(self.dirs_dict["MERGE_DIR"], g, "README.txt") \
                            for g in self.group_names if self.group_sizes[g] == 1]
        target_files.extend(message_file_for_groups_of_size_1)


        contigs_annotated = [os.path.join(self.dirs_dict["CONTIGS_DIR"],\
                             g + "-annotate_contigs_database.done") for g in self.group_names]
        target_files.extend(contigs_annotated)

        if self.run_qc:
            qc_report = os.path.join(self.dirs_dict["QC_DIR"], "qc-report.txt")
            target_files.append(qc_report)

        if self.references_for_removal_txt:
            filter_report = os.path.join(self.dirs_dict["QC_DIR"], "short-read-removal-report.txt")
            target_files.append(filter_report)

        if self.collections:
            for group in self.collections.keys():
                target_files.append(self.get_collection_import_flag(group))

        if self.run_summary:
            summary = [os.path.join(self.dirs_dict["SUMMARY_DIR"], g + "-SUMMARY") for g in self.collections.keys()]
            target_files.extend(summary)

        if self.run_split:
            split = [os.path.join(self.dirs_dict["SPLIT_PROFILES_DIR"],\
                                  g + "-split.done")\
                                  for g in self.collections.keys() if not 'default_collection' in self.collections[g]]
            target_files.extend(split)

        targets_files_for_binning = self.get_target_files_for_anvi_cluster_contigs()
        if targets_files_for_binning:
            target_files.extend(targets_files_for_binning)

        return target_files


    def get_collection_import_flag(self, group):
        ''' Return the flag for collection import (either default collection or from file).'''
        if not self.collections[group].get('default_collection'):
            flag = os.path.join(self.dirs_dict["MERGE_DIR"], group, 'collection-import.done')
        else:
            flag = os.path.join(self.dirs_dict["MERGE_DIR"], group, 'default-collection-import.done')

        return flag


    def init_references_txt(self):
        if self.references_mode:
            if not self.fasta_txt_file:
                raise ConfigError("In refrences mode, you need to also fill in the `fasta_txt` value in your config file.")

            if not filesnpaths.is_file_exists(self.fasta_txt_file, dont_raise=True):
                raise ConfigError('You know the path you have for `fasta_txt` in your config file? There is no such file on your disk :(')

        if not self.references_mode:
            # if it is reference mode then the group names have been assigned in the contigs Snakefile
            # if it is not reference mode and no groups are supplied in the samples_txt then group names are sample names
            self.group_names = self.sample_names

        if self.fasta_txt_file and not self.references_mode:
            raise ConfigError("In order to use reference fasta files you must set "
                          "\"'references_mode': true\" in your config file, yet "
                          "you didn't, but at the same time you supplied the following "
                          "fasta_txt: %s. So we don't know what to do with this "
                          "fasta_txt" % self.fasta_txt_file)

        # Collecting information regarding groups.
        if self.samples_txt.has_groups():
            # if groups were specified then members of a groups will be co-assembled.
            self.group_names = self.samples_txt.group_names()
            # creating a dictionary with groups as keys and number of samples in
            # the groups as values
            #self.group_sizes = self.samples_information['group'].value_counts().to_dict()
            self.group_sizes = self.samples_txt.group_sizes()

            if self.references_mode:
                # sanity check to see that groups specified in samples.txt match
                # the names of fasta.
                mismatch = set(self.group_names) - set(self.contigs_information.keys())
                if mismatch:
                    raise ConfigError("Group names specified in the samples.txt "
                                      "file must match the names of fasta "
                                      "in the fasta.txt file. These are the "
                                      "mismatches: %s" % mismatch)
                groups_in_contigs_information_but_not_in_samples_txt = set(self.contigs_information.keys()) - set(self.group_names)
                if groups_in_contigs_information_but_not_in_samples_txt:
                    run.warning('The following group names appear in your fasta_txt '
                                'but do not appear in your samples_txt. Maybe this is '
                                'ok with you, but we thought you should know. This means '
                                'that the metagenomics workflow will simply ignore these '
                                'groups.')

        else:
            if self.references_mode:
                # if the user didn't provide a group column in the samples.txt,
                # in references mode the default is 'all_against_all'.
                run.warning("No groups were provided in your samples_txt,\
                             hence 'all_against_all' mode has been automatically\
                             set to True.")
                self.set_config_param('all_against_all', True)
            else:
                # if no groups were specified then each sample would be assembled
                # separately
                run.warning("No groups were specified in your samples_txt. This is fine. "
                            "But we thought you should know. Any assembly will be performed "
                            "on individual samples (i.e. NO co-assembly).")
                # TODO: maybe SamplesTxt should automatically add a group colunn if none exists?
                self.samples_information['group'] = self.samples_information['sample']
                self.group_names = list(self.sample_names)
                self.group_sizes = dict.fromkeys(self.group_names,1)

        if self.get_param_value_from_config('all_against_all'):
            # in all_against_all, the size of each group is as big as the number
            # of samples.
            self.group_sizes = dict.fromkeys(self.group_names,len(self.sample_names))


        if not self.references_mode and not (self.get_param_value_from_config(['anvi_script_reformat_fasta','run']) == True):
            # in assembly mode (i.e. not in references mode) we always have
            # to run reformat_fasta. The only reason for this is that
            # the megahit output is temporary, and if we dont run
            # reformat_fasta we will delete the output of meghit at the
            # end of the workflow without saving a copy.
            raise ConfigError('You seem to be interested in running the metagenomics workflow in assembly mode. '
                              'In this mode you can\'t skip `anvi_script_reformat_fasta` rule, which means you need '
                              'to add this rule to your config file with `"run": true` directive for this to work.')


    def init_kraken(self):
        '''Making sure the sample names and file paths the provided kraken.txt file are valid'''
        kraken_txt = self.get_param_value_from_config('kraken_txt')
        self.run_krakenuniq = self.get_param_value_from_config(['krakenuniq', 'run']) == True

        if kraken_txt:
            if self.get_param_value_from_config(['krakenuniq', 'run']) == True:
                raise ConfigError("You supplied a kraken_txt file (\"%s\") but you set krakenuniq "
                                  "to run in the config file. anvi'o is confused and "
                                  "is officially going on a strike. Ok, let's clarify, "
                                  "having a kraken_txt file means you already ran krakenuniq "
                                  "and want us to use those results, and yet you set krakenuniq "
                                  "to run again? why? Ok, time to strike. Bye!" % kraken_txt)

            # if a kraken_txt was supplied then let's run kraken by default
            self.run_krakenuniq = True

            kraken_annotation_dict = u.get_TAB_delimited_file_as_dictionary(kraken_txt)
            if next(iter(next(iter(kraken_annotation_dict.values())).keys())) != "path":
                raise ConfigError("Your kraken annotation file, '%s', is not formatted properly "
                                  "anvi'o expected it to have two columns only and the second column "
                                  "should have a header 'path'." % kraken_txt)
            samples_in_kraken_txt = set(kraken_annotation_dict.keys())
            # get a list of the sample names
            sample_names = set(self.samples_information['sample'])

            wrong_samples_in_kraken_txt = samples_in_kraken_txt - sample_names
            if wrong_samples_in_kraken_txt:
                raise ConfigError("Your kraken annotation file, '%s', contains samples that "
                                  "are not in your samples_txt file, '%s'. Here is an example "
                                  "of such a sample: %s." % (kraken_txt, self.get_param_value_from_config('samples_txt'), next(iter(wrong_samples_in_kraken_txt))))

            missing_samples_in_kraken_txt = sample_names - samples_in_kraken_txt
            if missing_samples_in_kraken_txt:
                raise ConfigError("Your kraken annotation file, '%s', is missing samples that "
                                  "are in your samples_txt file, '%s'. This is not allowed. "
                                  "Here is an example of such a sample: %s." % (kraken_txt, self.get_param_value_from_config('samples_txt'), wrong_samples_in_kraken_txt[0]))
            self.kraken_annotation_dict = kraken_annotation_dict

        if self.get_param_value_from_config(['krakenuniq', 'run']):
            if not self.get_param_value_from_config(['krakenuniq', '--db']):
                raise ConfigError('In order to run krakenuniq, you must provide a path to '
                                  'a database using the --db parameter in the config file.')


    def load_collections(self):
        ''' Load the collections_txt file, run some sanity checks, and figure out params for anvi_import_collection'''
        collections = u.get_TAB_delimited_file_as_dictionary(self.collections_txt)
        bad_groups = [g for g in collections if g not in self.group_names]
        if bad_groups:
                raise ConfigError('Some of the names in your collection_txt '
                                  'file ("%s") don\'t match the names of the '
                                  'groups in your samples_txt/fasta_txt. '
                                  'Here are the names that don\'t match: %s. '
                                  'And here are the group names we expect to find: '
                                  '%s' % (self.collections_txt, ', '.join(bad_groups), ', '.join(self.group_names)))
        for group in collections:
            default_collection = collections[group].get('default_collection')

            if default_collection:
                # User can specify either a default collection OR collection from file
                not_allowed_params = {'collection_name', 'collection_file', 'bins_info', 'contigs_mode'}
                if any([collections[group][key] for key in not_allowed_params if key in collections[group].keys()]):
                    raise ConfigError('We encountered the following problem with your '
                                      'collections_txt file ("%s"): you can choose '
                                      'either using a default collection OR importing '
                                      'a collection from a file. Yet, for "%s", you specificy '
                                      'a default collection AND also specify some of the following '
                                      'parameters: %s.' % (self.collections_txt, group, ", ".join(not_allowed_params)))

                collections[group]['collection_name'] = 'DEFAULT'
                collections[group]['contigs_mode'] = ''

            else:
                if not filesnpaths.is_file_exists(collections[group]['collection_file'], dont_raise=True):
                    raise ConfigError('We encountered the following problem with your '
                                      'collections_txt file ("%s"): you did not specify '
                                      'a valid collection file for "%s".' % (self.collections_txt, group))

                if not collections[group]['collection_name']:
                    raise ConfigError('You must specify a name for each collection in your collections_txt')
                u.check_collection_name(collections[group]['collection_name'])
                if collections[group].get('bins_info'):
                    filesnpaths.is_file_exists(collections[group]['bins_info'])
                    collections[group]['bins_info'] = '--bins-info %s' % collections[group]['bins_info']
                else:
                    collections[group]['bins_info'] = ''
                if collections[group].get('contigs_mode'):
                    collections[group]['contigs_mode'] = '--contigs-mode'
                else:
                    collections[group]['contigs_mode'] = ''
        self.collections = collections


    def load_references_for_removal(self):
        """Load and perform some sanity checks on the references for removal"""
        self.references_for_removal = u.get_TAB_delimited_file_as_dictionary(self.references_for_removal_txt)
        # adding the references_for_removal to the fasta_information dict
        self.fasta_information.update(self.references_for_removal)

        for sample in self.references_for_removal.keys():
            try:
                u.check_sample_id(sample)
            except ConfigError as e:
                raise ConfigError("While processing the references for removal txt file ('%s'), anvi'o ran into the following error: "
                                  "%s" % (self.samples_txt_file, e))

        files_that_end_with_gz = []
        for ref_dict in self.references_for_removal.values():
            if 'path' not in ref_dict:
                raise ConfigError('Yor references for removal txt file is not formatted properly. It must have only two columns '
                                  'with the headers "reference" and "path".')
            if ref_dict['path'].endswith('.gz'):
                filesnpaths.is_file_exists(ref_dict['path'])
                files_that_end_with_gz.append(ref_dict['path'])
            else:
                # if the file is not compressed then we can verify that it is a fasta file
                filesnpaths.is_file_fasta_formatted(ref_dict['path'])

        if files_that_end_with_gz:
            run.warning('The following reference for removal files are compressed: %s. '
                        'That\'s fine, but it means that we will skip the '
                        'sanity check to verify that this is actually '
                        'a properly formatted fasta file. Things are '
                        'probably Ok, this is just one of these occasions '
                        'in which anvi\'o is oversharing.' % ', '.join(files_that_end_with_gz))

        if self.references_mode:
            # Make sure that the user didn't give the same name to references and references_for_removal
            ref_name_in_both = [r for r in self.references_for_removal if r in self.contigs_information]
            if ref_name_in_both:
                raise ConfigError('You must have unique names for your fasta files in your fasta txt file '
                                  'and your references for removal txt file. These are the names that appear '
                                  'in both: %s' % ', '.join(ref_name_in_both))
        dont_remove = self.get_param_value_from_config(['remove_short_reads_based_on_references', 'dont_remove_just_map'])
        if not dont_remove:
            self.remove_short_reads_based_on_references = True


    def get_assembly_software_list(self):
        return ['megahit', 'idba_ud', 'metaspades']


    def gen_report_with_references_for_removal_info(self, filtered_id_files, output_file_name):
        ''' If mapping was done to reference for removal then we create a report with the results.'''
        report_dict = {}
        for filename in filtered_id_files:
            sample = os.path.basename(filename).split("-ids-to-remove.txt")[0]
            ids = set(open(filename).read().splitlines())
            report_dict[sample] = {}
            report_dict[sample]['number_of_filtered_reads'] = len(ids)
        u.store_dict_as_TAB_delimited_file(report_dict, output_file_name, headers=["sample", 'number_of_filtered_reads'])


    def get_fasta(self, wildcards):
        if wildcards.group in self.references_for_removal:
            # if it's a reference for removal then we just want to use the
            # raw fasta file, and there is no need to reformat or assemble
            contigs = self.get_raw_fasta(wildcards)
        elif self.get_param_value_from_config(['anvi_script_reformat_fasta','run']):
            contigs = self.dirs_dict["FASTA_DIR"] + "/{group}/{group}-contigs.fa".format(group=wildcards.group)
        else:
            contigs = self.get_raw_fasta(wildcards)
        return contigs


    def get_raw_fasta(self, wildcards, remove_gz_suffix=True):
        if self.references_mode or wildcards.group in self.references_for_removal:
            # in 'reference mode' the input is the reference fasta
            contigs = super(MetagenomicsWorkflow, self).get_raw_fasta(wildcards, remove_gz_suffix=remove_gz_suffix)
        else:
            # by default the input fasta is the assembly output
            contigs = self.dirs_dict["FASTA_DIR"] + "/%s/final.contigs.fa" % wildcards.group
        return contigs


    def get_target_files_for_anvi_cluster_contigs(self):
        import anvio.workflows as w
        w.D(self.get_param_value_from_config(['anvi_cluster_contigs', 'run']))
        if self.get_param_value_from_config(['anvi_cluster_contigs', 'run']) is not True:
            # the user doesn't want to run this
            return
        w.D('hi2')
        requested_drivers = self.get_param_value_from_config(['anvi_cluster_contigs', '--driver'])
        if not requested_drivers:
            raise ConfigError('You must specify a driver for anvi_cluster_contigs. '
                              'You specified \'"run": true\' for anvi_cluster_contigs, \
                               but provided no driver.')

        if type(requested_drivers) != list:
            requested_drivers = [requested_drivers]

        incompatible_drivers = [d for d in requested_drivers if d not in list(driver_modules['binning'].keys())]
        if incompatible_drivers:
            raise ConfigError('The following drivers were listed in the config file for rule anvi_cluster_contigs '
                              'but they are not familiar to anvi\'o: %s' % ', '.join(incompatible_drivers))

        # TODO: we should make sure drivers are installed. Maybe Ozcan or Meren are willing to do this?

        for d in driver_modules['binning'].keys():
            # let's make sure no parameters were set for a driver that was not listed
            additional_parameters = self.get_param_value_from_config(['anvi_cluster_contigs', self.get_param_name_for_binning_driver(d)])
            if additional_parameters:
                if d not in requested_drivers:
                    raise ConfigError('You set the following parameters: "%s" for %s, but you did not '
                                     'specify it as one of the drivers that should be used by anvi_cluster_contigs. '
                                     'In order to reduce room for mistakes, we don\'t allow this.' % (additional_parameters, d))

        collection_name = self.get_param_value_from_config(['anvi_cluster_contigs', '--collection-name'])
        if '{driver}' not in collection_name:
            if len(requested_drivers) > 1:
                raise ConfigError('When using multiple binning algorithms, the --collection-name '
                                  'for rule anvi_cluster_contigs must contain '
                                  'the key word "{driver}" (including those curly brackets). '
                                  'It appears you changed the --collection-name, which by '
                                  'default is simply "{driver}" to: "%s". That\'s fine, \
                                   but only as long as you have "{driver}" appear somewhere in \
                                   the name, because otherwise the collections made by different \
                                   binning algorithms would try to override each other, since they \
                                   would all end up having the same name' % collection_name)

            # if the key word '{driver}' is not in the collection name then it is a static collection name
            example_collection_name = collection_name
        else:
            # if the key word '{driver}' IS in the collection name, then let's take one
            # driver as example when we check if the collection name is valid.
            example_collection_name = collection_name.format(driver=requested_drivers[0])
        try:
            u.check_collection_name(example_collection_name)
        except ConfigError as e:
            raise ConfigError('%s is not an acceptable collection name for anvi_cluster_contigs. '
                  'We tried it for one of the drivers you requested and this is the '
                  'error we got: %s' % (collection_name, e))

        groups_of_size_one = []
        target_files = []
        for d in requested_drivers:
            for g in self.group_names:
                if self.group_sizes[g] > 1:
                    # add flag files of binning to the target files
                    # only if the group is larger than 1 (because otherwise, there is no merged profile)
                    target_files.append(self.get_flag_file_for_binning_driver(g, d))
                else:
                    groups_of_size_one.append(g)
        if groups_of_size_one:
            run.warning('You requested to run anvi_cluster_contigs, but it will not run for '
                        '%s, since there is only one sample in this group, and hence there '
                        'will be no merged profile, which is required for anvi-cluster-contigs.' % ', '.join(groups_of_size_one))

        return target_files


    def get_flag_file_for_binning_driver(self, group, driver):
        return os.path.join(self.dirs_dict['MERGE_DIR'], group + '-' + driver + '.done')


    def get_param_name_for_binning_driver(self, driver):
        return '--additional-params' + '-' + driver
