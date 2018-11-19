# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o pangenomics workflows.
"""

import anvio
import anvio.terminal as terminal

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


class PangenomicsWorkflow(ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='pangenomics')

        # initialize the base class
        ContigsDBWorkflow.__init__(self)

        self.rules.extend(['anvi_gen_genomes_storage',
                           'anvi_pan_genome'])

        self.general_params.extend(["project_name",
                                    "fasta_txt",
                                    "internal_genomes",
                                    "external_genomes"])

        self.dirs_dict.update({"FASTA_DIR": "01_FASTA",
                               "CONTIGS_DIR": "02_CONTIGS",
                               "PAN_DIR": "03_PAN"})

        self.default_config.update({"fasta_txt": "fasta.txt",
                                    "anvi_pan_genome": {"threads": 7}})

        pan_params = ["--project-name", "--genome-names", "--skip-alignments",\
                     "--align-with", "--exclude-partial-gene-calls", "--use-ncbi-blast",\
                     "--minbit", "--mcl-inflation", "--min-occurrence",\
                     "--min-percent-identity", "--sensitive", "--description",\
                     "--overwrite-output-destinations", "--skip-hierarchical-clustering",\
                     "--enforce-hierarchical-clustering", "--distance", "--linkage"]
        self.rule_acceptable_params_dict['anvi_pan_genome'] = pan_params

        storage_params = ["--gene-caller"]
        self.rule_acceptable_params_dict['anvi_gen_genomes_storage'] = storage_params


    def init(self):
        ''' backhand stuff (mostly sanity checks) specific for the phylogenomics workflow'''
        super().init()
        self.internal_genomes_file = self.get_param_value_from_config('internal_genomes')
        self.external_genomes_file = self.get_param_value_from_config('external_genomes')
        self.input_for_anvi_gen_genomes_storage = self.get_internal_and_external_genomes_files()
        self.project_name = self.get_param_value_from_config("project_name")
        self.pan_project_name = self.get_param_value_from_config(["anvi_pan_genome", "--project-name"])

        self.sanity_checks()

    def sanity_checks(self):
        if (not self.internal_genomes_file) and (not self.external_genomes_file):
            raise ConfigError("You must provide a path to either internal_genomes_file or external_genomes_file\
                                   or both.")
        if not self.project_name:
            raise ConfigError("You must provide a project name in your config file.")

        if self.pan_project_name:
            run.warning('you chose to set the "--project-name" parameter for "anvi_pan_genome". That is ok\
                         but just so you know, if you haven\'t supplied this, then we would have taken the value\
                         from "project_name" in your config file to also be the project name for "anvi_pan_genome"')
        else:
            self.pan_project_name = self.project_name
