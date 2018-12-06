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
from anvio.errors import ConfigError, FilesNPathsError
from anvio.workflows import WorkflowSuperClass
from anvio.workflows.contigs import ContigsDBWorkflow


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
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

        self.target_files = [] # TODO: Once we update all other workflows then this will be initiated in WorkflowSuperClass
        self.samples_information = {}
        self.kraken_annotation_dict = {}
        self.run_krakenhll = None
        self.run_metaspades = None
        self.use_scaffold_from_metaspades = None
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
                     'anvi_gen_contigs_database', 'anvi_export_gene_calls', 'centrifuge',\
                     'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam', 'idba_ud',\
                     'anvi_profile', 'annotate_contigs_database', 'anvi_merge', 'import_percent_of_reads_mapped',\
                     'krakenhll', 'krakenhll_mpa_report', 'import_kraken_hll_taxonomy', 'metaspades',\
                     'remove_short_reads_based_on_references', 'anvi_summarize', 'anvi_split'])

        self.general_params.extend(['samples_txt', "references_mode", "all_against_all",\
                                    "kraken_txt", "collections_txt"])

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule
        rule_acceptable_params_dict['iu_gen_configs'] = ["--r1-prefix", "--r2-prefix"]
        rule_acceptable_params_dict['iu_filter_quality_minoche'] = ['run', '--visualize-quality-curves', '--ignore-deflines', '--limit-num-pairs', '--print-qual-scores', '--store-read-fate']
        rule_acceptable_params_dict['gzip_fastqs'] = ["run"]
        rule_acceptable_params_dict['anvi_summarize'] = ["additional_params", "run"]
        rule_acceptable_params_dict['anvi_split'] = ["additional_params", "run"]
        rule_acceptable_params_dict['metaspades'] = ["run", "additional_params", "use_scaffolds"]
        rule_acceptable_params_dict['megahit'] = ["run", "--min-contig-len", "--min-count", "--k-min",
                                                  "--k-max", "--k-step", "--k-list",
                                                  "--no-mercy", "--no-bubble", "--merge-level",
                                                  "--prune-level", "--prune-depth", "--low-local-ratio",
                                                  "--max-tip-len", "--no-local", "--kmin-1pass",
                                                  "--presets", "--memory", "--mem-flag",
                                                  "--use-gpu", "--gpu-mem", "--keep-tmp-files",
                                                  "--tmp-dir", "--continue", "--verbose"]
        rule_acceptable_params_dict['idba_ud'] = ["run", "--mink", "--maxk", "--step", "--inner_mink",
                                                  "--inner_step", "--prefix", "--min_count",
                                                  "--min_support", "--seed_kmer", "--min_contig",
                                                  "--similar", "--max_mismatch", "--min_pairs",
                                                  "--no_bubble", "--no_local", "--no_coverage",
                                                  "--no_correct", "--pre_correction"]
        rule_acceptable_params_dict['bowtie'] = ["additional_params"]
        rule_acceptable_params_dict['samtools_view'] = ["additional_params"]
        rule_acceptable_params_dict['anvi_profile'] = ["--overwrite-output-destinations", "--sample-name", "--report-variability-full",
                                                        "--skip-SNV-profiling", "--profile-SCVs", "--description",
                                                        "--skip-hierarchical-clustering", "--distance", "--linkage", "--min-contig-length",
                                                        "--min-mean-coverage", "--min-coverage-for-variability", "--cluster-contigs",
                                                        "--contigs-of-interest", "--queue-size", "--write-buffer-size", "--max-contig-length"]
        rule_acceptable_params_dict['annotate_contigs_database'] = []
        rule_acceptable_params_dict['merge_fastas_for_co_assembly'] = []
        rule_acceptable_params_dict['merge_fastqs_for_co_assembly'] = []
        rule_acceptable_params_dict['anvi_merge'] = ["--sample-name", "--description", "--skip-hierarchical-clustering",
                                                     "--enforce-hierarchical-clustering", "--distance", "--linkage",
                                                     "--skip-concoct-binning", "--overwrite-output-destinations"]
        rule_acceptable_params_dict['import_percent_of_reads_mapped'] = ["run"]
        rule_acceptable_params_dict['krakenhll'] = ["additional_params", "run", "--db", "--gzip-compressed"]
        rule_acceptable_params_dict['krakenhll_mpa_report'] = ["additional_params"]
        rule_acceptable_params_dict['import_kraken_hll_taxonomy'] = ["--min-abundance"]
        rule_acceptable_params_dict['remove_short_reads_based_on_references'] = ["dont_remove_just_map", \
                                                                                 "references_for_removal_txt", \
                                                                                 "delimiter-for-iu-remove-ids-from-fastq"]

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        forbidden_params = {}
        forbidden_params['krakenhll'] = ['--fastq-input', '--paired', '--output']

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
                                    'metaspades': {"additional_params": "--only-assembler", "threads": 7},
                                    'megahit': {"--min-contig-len": min_contig_length_for_assembly, "--memory": 0.4, "threads": 7},
                                    'idba_ud': {"--min_contig": min_contig_length_for_assembly, "threads": 7},
                                    'iu_filter_quality_minoche': {"run": True, "--ignore-deflines": True},
                                    "gzip_fastqs": {"run": True},
                                    "bowtie": {"additional_params": "--no-unal", "threads": 3},
                                    "samtools_view": {"additional_params": "-F 4"},
                                    "anvi_profile": {"threads": 3, "--sample-name": "{sample}", "--overwrite-output-destinations": True},
                                    "anvi_merge": {"--sample-name": "{group}", "--overwrite-output-destinations": True},
                                    "import_percent_of_reads_mapped": {"run": True},
                                    "krakenhll": {"threads": 3, "--gzip-compressed": True, "additional_params": "--preload"},
                                    "remove_short_reads_based_on_references": {"delimiter-for-iu-remove-ids-from-fastq": " "}})


    def init(self):
        super().init()

        # loading the samples.txt file
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(self.samples_txt_file)
        # getting the samples information (names, [group], path to r1, path to r2) from samples.txt
        self.samples_information = pd.read_csv(self.samples_txt_file, sep='\t', index_col=False)


        # get a list of the sample names
        self.sample_names = list(self.samples_information['sample'])
        self.run_metaspades = self.get_param_value_from_config(['metaspades', 'run'])
        self.use_scaffold_from_metaspades = self.get_param_value_from_config(['metaspades', 'use_scaffolds'])
        self.run_qc = self.get_param_value_from_config(['iu_filter_quality_minoche', 'run']) == True
        self.run_summary = self.get_param_value_from_config(['anvi_summarize', 'run']) == True
        self.run_split = self.get_param_value_from_config(['anvi_split', 'run']) == True
        self.references_mode = self.get_param_value_from_config('references_mode', repress_default=True)
        self.fasta_txt_file = self.get_param_value_from_config('fasta_txt', repress_default=True)

        self.references_for_removal_txt = self.get_param_value_from_config(['remove_short_reads_based_on_references',\
                                                                            'references_for_removal_txt'],\
                                                                           repress_default=True)
        if self.references_for_removal_txt:
            self.load_references_for_removal()

        self.collections_txt = self.get_param_value_from_config('collections_txt')
        if self.collections_txt:
            self.load_collections()
        elif self.run_summary:
            raise ConfigError('If you want to run anvi-summarize you must provide a collections_txt file')
        elif self.run_split:
            raise ConfigError('If you want to run anvi-split you must provide a collections_txt file')

        self.init_target_files()


    def init_target_files(self):
        target_files = []
        merged_profiles = [os.path.join(self.dirs_dict["MERGE_DIR"], g, "PROFILE.db") for g in self.group_names]
        target_files.extend(merged_profiles)

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
            import_collection_done = [os.path.join(self.dirs_dict["MERGE_DIR"],\
                                                   g,\
                                                   'collection-import.done')\
                                                   for g in self.collections.keys()]

            target_files.extend(import_collection_done)

        if self.run_summary:
            summary = [os.path.join(self.dirs_dict["SUMMARY_DIR"], g + "-SUMMARY") for g in self.collections.keys()]
            target_files.extend(summary)

        if self.run_split:
            split = [os.path.join(self.dirs_dict["SPLIT_PROFILES_DIR"],\
                                  g + "-split.done")\
                                  for g in self.collections.keys()]
            target_files.extend(split)

        self.target_files.extend(target_files)


    def sanity_checks(self):
        self.sanity_check_for_samples_txt()
        self.sanity_check_for_kraken()
        self.sanity_check_for_refereces_txt()


    def sanity_check_for_refereces_txt(self):
        if self.references_mode:
            try:
                filesnpaths.is_file_exists(self.fasta_txt_file)
            except FilesNPathsError as e:
                raise ConfigError('In references mode you must supply a fasta_txt file.')


        if not self.references_mode:
            # if it is reference mode then the group names have been assigned in the contigs Snakefile
            # if it is not reference mode and no groups are supplied in the samples_txt then group names are sample names
            self.group_names = self.sample_names

        if self.fasta_txt_file and not self.references_mode:
            raise ConfigError("In order to use reference fasta files you must set\
                           \"'references_mode': true\" in your config file, yet\
                           you didn't, but at the same time you supplied the following\
                           fasta_txt: %s. So we don't know what to do with this\
                           fasta_txt" % self.fasta_txt_file)

        # Collecting information regarding groups.
        if "group" in self.samples_information.columns:
            # if groups were specified then members of a groups will be co-assembled.
            self.group_names = list(self.samples_information['group'].unique())
            # creating a dictionary with groups as keys and number of samples in
            # the groups as values
            self.group_sizes = self.samples_information['group'].value_counts().to_dict()

            if self.references_mode:
                # sanity check to see that groups specified in samples.txt match
                # the names of fasta.
                mismatch = set(self.group_names) - set(self.contigs_information.keys())
                if mismatch:
                    raise ConfigError("Group names specified in the samples.txt \
                                       file must match the names of fasta \
                                       in the fasta.txt file. These are the \
                                       mismatches: %s" % mismatch)
                groups_in_contigs_information_but_not_in_samples_txt = set(self.contigs_information.keys()) - set(self.group_names)
                if groups_in_contigs_information_but_not_in_samples_txt:
                    run.warning('The following group names appear in your fasta_txt\
                                 but do not appear in your samples_txt. Maybe this is\
                                 ok with you, but we thought you should know. This means\
                                 that the metagenomics workflow will simply ignore these\
                                 groups.')

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
                run.warning("No groups were specified in your samples_txt. This is fine.\
                             But we thought you should know. Any assembly will be performed\
                             on individual samples (i.e. NO co-assembly).")
                self.samples_information['group'] = self.samples_information['sample']
                self.group_names = list(self.sample_names)
                self.group_sizes = dict.fromkeys(self.group_names,1)

        if self.get_param_value_from_config('all_against_all', repress_default=True):
            # in all_against_all, the size of each group is as big as the number
            # of samples.
            self.group_sizes = dict.fromkeys(self.group_names,len(self.sample_names))


        if not self.references_mode and not (self.get_param_value_from_config(['anvi_script_reformat_fasta','run']) == True):
            # in assembly mode (i.e. not in references mode) we always have
            # to run reformat_fasta. The only reason for this is that
            # the megahit output is temporary, and if we dont run
            # reformat_fasta we will delete the output of meghit at the
            # end of the workflow without saving a copy.
            raise ConfigError("You can't skip reformat_fasta in assembly mode \
                                please change your config.json file")


    def sanity_check_for_samples_txt(self):
        if 'sample' not in self.samples_information.columns.values:
            raise ConfigError("You know what. This '%s' file does not look anything like\
                               a samples file." % self.samples_txt_file)

        for sample in self.samples_information['sample']:
            try:
                u.check_sample_id(sample)
            except ConfigError as e:
                raise ConfigError("While processing the samples txt file ('%s'), anvi'o ran into the following error: \
                                   %s" % (self.samples_txt_file, e))

        fastq_file_names = list(self.samples_information['r1']) + list(self.samples_information['r2'])
        bad_fastq_names = [s for s in fastq_file_names if (not s.endswith('.fastq') and not s.endswith('.fastq.gz'))]
        if bad_fastq_names:
            raise ConfigError("We require tha all fastq file names end with either '.fastq' \
                               or '.fastq.gz'. Some or all of the file names in %s aren't formatted \
                               accordingly. These are the file names we don't like: %s" % (self.samples_txt_file, ', '.join(bad_fastq_names)))


    def sanity_check_for_kraken(self):
        '''Making sure the sample names and file paths the provided kraken.txt file are valid'''
        kraken_txt = self.get_param_value_from_config('kraken_txt')
        self.run_krakenhll = self.get_param_value_from_config(['krakenhll', 'run']) == True

        if kraken_txt:
            if self.get_param_value_from_config(['krakenhll', 'run']) == True:
                raise ConfigError("You supplied a kraken_txt file (\"%s\") but you set krakenhll \
                                   to run in the config file. anvi'o is confused and \
                                   is officially going on a strike. Ok, let's clarify, \
                                   having a kraken_txt file means you already ran krakenhll \
                                   and want us to use those results, and yet you set krakenhll \
                                   to run again? why? Ok, time to strike. Bye!" % kraken_txt)

            # if a kraken_txt was supplied then let's run kraken by default
            self.run_krakenhll = True

            kraken_annotation_dict = u.get_TAB_delimited_file_as_dictionary(kraken_txt)
            if next(iter(next(iter(kraken_annotation_dict.values())).keys())) != "path":
                raise ConfigError("Your kraken annotation file, '%s', is not formatted properly \
                                   anvi'o expected it to have two columns only and the second column \
                                   should have a header 'path'." % kraken_txt)
            samples_in_kraken_txt = set(kraken_annotation_dict.keys())
            # get a list of the sample names
            sample_names = set(self.samples_information['sample'])

            wrong_samples_in_kraken_txt = samples_in_kraken_txt - sample_names
            if wrong_samples_in_kraken_txt:
                raise ConfigError("Your kraken annotation file, '%s', contains samples that \
                                   are not in your samples_txt file, '%s'. Here is an example \
                                   of such a sample: %s." % (kraken_txt, self.get_param_value_from_config('samples_txt'), next(iter(wrong_samples_in_kraken_txt))))

            missing_samples_in_kraken_txt = sample_names - samples_in_kraken_txt
            if missing_samples_in_kraken_txt:
                raise ConfigError("Your kraken annotation file, '%s', is missing samples that \
                                   are in your samples_txt file, '%s'. This is not allowed. \
                                   Here is an example of such a sample: %s." % (kraken_txt, self.get_param_value_from_config('samples_txt'), wrong_samples_in_kraken_txt[0]))
            self.kraken_annotation_dict = kraken_annotation_dict

        if self.get_param_value_from_config(['krakenhll', 'run']):
            if not self.get_param_value_from_config(['krakenhll', '--db']):
                raise ConfigError('In order to run krakenhll, you must provide a path to \
                                   a database using the --db parameter in the config file.')


    def load_collections(self):
        ''' Load the collections_txt file, run some sanity checks, and figure out params for anvi_import_collection'''
        collections = u.get_TAB_delimited_file_as_dictionary(self.collections_txt)
        bad_groups = [g for g in collections if g not in self.group_names]
        if bad_groups:
                raise ConfigError('Some of the names in your collection_txt \
                                   file ("%s") don\'t match the names of the \
                                   groups in your samples_txt/fasta_txt. \
                                   Here are the names that don\'t match: %s. \
                                   And here are the group names we expect to find: \
                                   %s' % (self.collections_txt, ', '.join(bad_groups), ', '.join(self.group_names)))
        for group in collections:
            filesnpaths.is_file_exists(collections[group]['collection_file'])
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
                raise ConfigError("While processing the references for removal txt file ('%s'), anvi'o ran into the following error: \
                                   %s" % (self.samples_txt_file, e))

        files_that_end_with_gz = []
        for ref_dict in self.references_for_removal.values():
            if 'path' not in ref_dict:
                raise ConfigError('Yor references for removal txt file is not formatted properly. It must have only two columns \
                                   with the headers "reference" and "path".')
            if ref_dict['path'].endswith('.gz'):
                filesnpaths.is_file_exists(ref_dict['path'])
                files_that_end_with_gz.append(ref_dict['path'])
            else:
                # if the file is not compressed then we can verify that it is a fasta file
                filesnpaths.is_file_fasta_formatted(ref_dict['path'])

        if files_that_end_with_gz:
            run.warning('The following reference for removal files are compressed: %s. \
                         That\'s fine, but it means that we will skip the \
                         sanity check to verify that this is actually \
                         a properly formatted fasta file. Things are \
                         probably Ok, this is just one of these occasions \
                         in which anvi\'o is oversharing.' % ', '.join(files_that_end_with_gz))

        if self.references_mode:
            # Make sure that the user didn't give the same name to references and references_for_removal
            ref_name_in_both = [r for r in self.references_for_removal if r in self.contigs_information]
            if ref_name_in_both:
                raise ConfigError('You must have unique names for your fasta files in your fasta txt file \
                                   and your references for removal txt file. These are the names that appear \
                                   in both: %s' % ', '.join(ref_name_in_both))
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
