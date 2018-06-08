# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o contigs workflows.
"""


import anvio
import anvio.terminal as terminal

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

        # initialize the base class
        ContigsDBWorkflow.__init__(self)

        self.rules.extend(['iu_gen_configs', 'iu_filter_quality_minoche', 'gen_qc_report', 'gzip_fastqs',\
                     'fq2fa', 'merge_fastas_for_co_assembly', 'megahit',\
                     'anvi_gen_contigs_database', 'anvi_export_gene_calls', 'centrifuge',\
                     'anvi_import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
                     'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam', 'idba_ud', \
                     'anvi_profile', 'annotate_contigs_database', 'anvi_merge', 'import_percent_of_reads_mapped'])

        self.general_params.extend(["samples_txt", "references_mode", "all_against_all"])

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
                                                        "--contigs-of-interest", "--queue-size", "--write-buffer-size"]
        rule_acceptable_params_dict['annotate_contigs_database'] = []
        rule_acceptable_params_dict['anvi_merge'] = ["--sample-name", "--description", "--skip-hierarchical-clustering",
                                                     "--enforce-hierarchical-clustering", "--distance", "--linkage",
                                                     "--skip-concoct-binning", "--overwrite-output-destinations"]
        rule_acceptable_params_dict['import_percent_of_reads_mapped'] = ["run"]

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.dirs_dict.update({"QC_DIR": "01_QC",
                               "FASTA_DIR": "02_FASTA",
                               "CONTIGS_DIR": "03_CONTIGS",
                               "MAPPING_DIR": "04_MAPPING",
                               "PROFILE_DIR": "05_ANVIO_PROFILE",
                               "MERGE_DIR": "06_MERGED"})

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
                                    "import_percent_of_reads_mapped": {"run": True}})


    def get_assembly_software_list(self):
        return ['megahit', 'idba_ud']
