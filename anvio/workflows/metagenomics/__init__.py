# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import anvio
import pandas as pd
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio import utils as u
from anvio.errors import ConfigError
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
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # know thyself.
        self.name = 'metagenomics'
        self.samples_information = {}
        self.kraken_annotation_dict = {}

        # initialize the base class
        ContigsDBWorkflow.__init__(self)

        self.rules.extend(['iu_gen_configs', 'iu_filter_quality_minoche', 'gen_qc_report', 'gzip_fastqs',\
                     'fq2fa', 'merge_fastas_for_co_assembly', 'megahit',\
                     'anvi_gen_contigs_database', 'anvi_export_gene_calls', 'centrifuge',\
                     'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam', 'idba_ud', \
                     'anvi_profile', 'annotate_contigs_database', 'anvi_merge', 'import_percent_of_reads_mapped', \
                     'krakenhll', 'krakenhll_mpa_report', 'import_kraken_hll_taxonomy'])

        self.general_params.extend(["samples_txt", "references_mode", "all_against_all", \
                                    "kraken_txt"])

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule
        rule_acceptable_params_dict['iu_gen_configs'] = ["--r1-prefix", "--r2-prefix"]
        rule_acceptable_params_dict['iu_filter_quality_minoche'] = ['run', '--visualize-quality-curves', '--ignore-deflines', '--limit-num-pairs', '--print-qual-scores', '--store-read-fate']
        rule_acceptable_params_dict['gzip_fastqs'] = ["run"]
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
        rule_acceptable_params_dict['anvi_merge'] = ["--sample-name", "--description", "--skip-hierarchical-clustering",
                                                     "--enforce-hierarchical-clustering", "--distance", "--linkage",
                                                     "--skip-concoct-binning", "--overwrite-output-destinations"]
        rule_acceptable_params_dict['import_percent_of_reads_mapped'] = ["run"]
        rule_acceptable_params_dict['krakenhll'] = ["additional_params", "run", "--db", "--gzip-compressed"]
        rule_acceptable_params_dict['krakenhll_mpa_report'] = ["additional_params"]
        rule_acceptable_params_dict['import_kraken_hll_taxonomy'] = ["--min-abundance"]

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
                               "TAXONOMY_DIR": "07_TAXONOMY"})

        self.default_config.update({'samples_txt': "samples.txt",
                                    'megahit': {"--min-contig-len": min_contig_length_for_assembly, "--memory": 0.4, "threads": 11},
                                    'idba_ud': {"--min_contig": min_contig_length_for_assembly, "threads": 11},
                                    'iu_filter_quality_minoche': {"run": True, "--ignore-deflines": True, "threads": 2},
                                    "gzip_fastqs": {"run": True},
                                    "bowtie_build": {"threads": 10},
                                    "bowtie": {"additional_params": "--no-unal", "threads": 10},
                                    "samtools_view": {"additional_params": "-F 4", "threads": 4},
                                    "anvi_init_bam": {"threads": 4},
                                    "anvi_profile": {"threads": 5, "--sample-name": "{sample}", "--overwrite-output-destinations": True},
                                    "anvi_merge": {"--sample-name": "{group}", "--overwrite-output-destinations": True},
                                    "import_percent_of_reads_mapped": {"run": True},
                                    "krakenhll": {"threads": 12, "--gzip-compressed": True, "additional_params": "--preload"}})


    def init(self):
        super().init()

        # loading the samples.txt file
        samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        filesnpaths.is_file_exists(samples_txt_file)
        # getting the samples information (names, [group], path to r1, path to r2) from samples.txt
        self.samples_information = pd.read_csv(samples_txt_file, sep='\t', index_col=False)

        if 'sample' not in self.samples_information.columns.values:
            raise ConfigError("You know what. This '%s' file does not look anything like\
                               a samples file." % samples_txt_file)

        for sample in self.samples_information['sample']:
            try:
                u.check_sample_id(sample)
            except ConfigError as e:
                raise ConfigError("While processing the samples txt file ('%s'), anvi'o ran into the following error: \
                                   %s" % (samples_txt_file, e))

        self.sanity_check_for_kraken()


    def sanity_check_for_kraken(self):
        '''Making sure the sample names and file paths the provided kraken.txt file are valid'''
        kraken_txt = self.get_param_value_from_config('kraken_txt')

        if kraken_txt:
            if self.get_param_value_from_config(['krakenhll', 'run']) == False:
                raise ConfigError("You supplied a kraken_txt file, %s, but you set krakenhll \
                                   not to run in the config file. anvi'o is confused and \
                                   is officially going on a strike." % kraken_txt)

            if 'krakenhll' not in self.config:
                raise ConfigError('You provided a kraken_txt, but you didnt set any parameters \
                                   for krakenhll. As a minimum, you must provide the path to \
                                   the krakenhll database using the --db parameter in the config file.')

            # if a kraken_txt was supplied then let's run kraken by default
            self.config['krakenhll']['run'] = True

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
        

    def get_assembly_software_list(self):
        return ['megahit', 'idba_ud']
