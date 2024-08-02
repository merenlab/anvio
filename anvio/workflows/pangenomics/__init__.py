# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o pangenomics workflows.
"""

import os

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.workflows.contigs import ContigsDBWorkflow
from anvio.workflows.phylogenomics import PhylogenomicsWorkflow


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class PangenomicsWorkflow(PhylogenomicsWorkflow, ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='pangenomics')

        self.pan_project_name = None
        self.valid_sequence_sources_for_phylogeny = ['gene_clusters', 'hmm']
        self.sequence_source_for_phylogeny = None
        self.tree_name = None

        # initialize the base class
        PhylogenomicsWorkflow.__init__(self)

        self.rules.extend(['anvi_gen_genomes_storage',
                           'anvi_pan_genome',
                           'anvi_get_sequences_for_gene_clusters',
                           'import_phylogenetic_tree_to_pangenome',
                           'anvi_compute_genome_similarity'])

        self.general_params.extend(["project_name",
                                    "fasta_txt",
                                    "internal_genomes",
                                    "external_genomes",
                                    "sequence_source_for_phylogeny"])

        self.dirs_dict.update({"FASTA_DIR": "01_FASTA",
                               "CONTIGS_DIR": "02_CONTIGS",
                               "PAN_DIR": "03_PAN"})

        self.default_config.update({"fasta_txt": "fasta.txt",
                                    "anvi_pan_genome": {"threads": 7},
                                    "import_phylogenetic_tree_to_pangenome": {'tree_name': 'phylogeny'},
                                    "anvi_compute_genome_similarity": {"run": False}})

        pan_params = ["--project-name", "--genome-names", "--skip-alignments",\
                     "--align-with", "--exclude-partial-gene-calls", "--use-ncbi-blast",\
                     "--minbit", "--mcl-inflation", "--min-occurrence",\
                     "--min-percent-identity", "--description",\
                     "--overwrite-output-destinations", "--skip-hierarchical-clustering",\
                     "--enforce-hierarchical-clustering", "--distance", "--linkage"]
        self.rule_acceptable_params_dict['anvi_pan_genome'] = pan_params

        storage_params = ["--gene-caller"]
        self.rule_acceptable_params_dict['anvi_gen_genomes_storage'] = storage_params

        seq_params = ["--gene-cluster-id", "--gene-cluster-ids-file",
                      "--collection-name", "--bin-id",
                      "--min-num-genomes-gene-cluster-occurs", "--max-num-genomes-gene-cluster-occurs",
                      "--min-num-genes-from-each-genome", "--max-num-genes-from-each-genome",
                      "--max-num-gene-clusters-missing-from-genome", "--min-functional-homogeneity-index",
                      "--max-functional-homogeneity-index", "--min-geometric-homogeneity-index",
                      "--max-geometric-homogeneity-index", "--add-into-items-additional-data-table",
                      "--concatenate-gene-clusters", "--separator", "--align-with"]
        self.rule_acceptable_params_dict['anvi_get_sequences_for_gene_clusters'] = seq_params

        import_params = ['--just-do-it', 'tree_name']
        self.rule_acceptable_params_dict['import_phylogenetic_tree_to_pangenome'] = import_params

        self.rule_acceptable_params_dict['anvi_compute_genome_similarity'] = ['run', 'additional_params']


    def init(self):
        ''' backend stuff (mostly sanity checks) specific for the phylogenomics workflow'''
        super().init()
        self.internal_genomes_file = self.get_param_value_from_config('internal_genomes')
        self.external_genomes_file = self.get_param_value_from_config('external_genomes')
        self.input_for_anvi_gen_genomes_storage = self.get_internal_and_external_genomes_files()
        self.project_name = self.get_param_value_from_config("project_name")
        self.pan_project_name = self.get_param_value_from_config(["anvi_pan_genome", "--project-name"])

        if self.pan_project_name:
            run.warning("You chose to set the '--project-name' parameter for 'anvi_pan_genome'. That is OK. "
                        "But just so you know, if you haven't supplied this, then we would have taken the value "
                        "from 'project_name' in your config file to also be the project name for 'anvi_pan_genome'.")
        else:
            self.pan_project_name = self.project_name

        self.sequence_source_for_phylogeny = self.get_param_value_from_config('sequence_source_for_phylogeny')
        if self.sequence_source_for_phylogeny == 'gene_clusters':
            GC_sequences = os.path.join(self.dirs_dict["PHYLO_DIR"], self.project_name + "-GC-sequences.fa")
            self.use_hmms_for_phylogeny = False
            self.phylogenomics_sequence_file = GC_sequences

        self.tree_name = self.get_param_value_from_config(['import_phylogenetic_tree_to_pangenome', 'tree_name'])
        self.pan_db_path = os.path.join(self.dirs_dict["PAN_DIR"], self.pan_project_name + "-PAN.db")
        self.input_for_anvi_compute_genome_similarity = {"pan_db": self.pan_db_path}
        self.input_for_anvi_compute_genome_similarity.update(self.get_internal_and_external_genomes_files())
        self.anvi_compute_genome_similarity_flag = os.path.join(self.dirs_dict["PAN_DIR"], self.project_name + "anvi_compute_genome_similarity.done")
        self.anvi_compute_genome_similarity_output_dir = os.path.join(self.dirs_dict["PAN_DIR"], self.project_name + "-ANI-OUTPUT")


    def get_pangenomics_target_files(self):
        target_files = []
        target_files.append(os.path.join(self.dirs_dict["PAN_DIR"], self.pan_project_name + "-PAN.db"))

        if self.sequence_source_for_phylogeny:
            target_files.append(self.get_phylogeny_imported_flag())

        if self.get_param_value_from_config(['anvi_compute_genome_similarity', 'run']):
            target_files.append(self.anvi_compute_genome_similarity_flag)

        return target_files


    def sanity_checks(self):
        if (not self.internal_genomes_file) and (not self.external_genomes_file):
            raise ConfigError("You must provide a path to either internal_genomes_file or external_genomes_file "
                                  "or both.")
        if not self.project_name:
            raise ConfigError("You must provide a project name in your config file.")

        if self.sequence_source_for_phylogeny:
            if self.sequence_source_for_phylogeny not in self.valid_sequence_sources_for_phylogeny:
                raise ConfigError('%s is not a valid sequence_source_for_phylogeny. '
                                  'We only know: %s' % (self.sequence_source_for_phylogeny,\
                                   ', '.join(self.valid_sequence_sources_for_phylogeny)))


    def get_phylogeny_imported_flag(self):
        return os.path.join(self.dirs_dict["PAN_DIR"], self.project_name + '-' + self.tree_name + "-phylogeny-imported.done")
