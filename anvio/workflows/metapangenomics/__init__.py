# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o metapangenomics workflows.
"""


import anvio
import anvio.terminal as terminal

from anvio.workflows import WorkflowSuperClass
from anvio.workflows import init_workflow_super_class
from anvio.workflows.contigs import ContigsDBWorkflow
from anvio.workflows.metagenomics import MetagenomicsWorkflow
from anvio.workflows.pangenomics import PangenomicsWorkflow


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class MetaPangenomicsWorkflow(MetagenomicsWorkflow, PangenomicsWorkflow, ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        init_workflow_super_class(self, args, workflow_name='metapangenomics')

        self.run = run
        self.progress = progress

        # know thyself.
        self.name = 'metapangenomics'

        # initialize the base classes
        PangenomicsWorkflow.__init__(self)
        MetagenomicsWorkflow.__init__(self)

        self.rules.extend([])

        self.general_params.extend(['metapangenome_fastas_txt'])

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule
        rule_acceptable_params_dict[''] = []

        # updating the dict from upstream
        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.dirs_dict.update({"QC_DIR": "01_QC",
                               "FASTA_DIR": "02_FASTA",
                               "CONTIGS_DIR": "03_CONTIGS",
                               "MAPPING_DIR": "04_MAPPING",
                               "PROFILE_DIR": "05_ANVIO_PROFILE",
                               "MERGE_DIR": "06_MERGED"})

        self.default_config.update({'metapangenome_fastas_txt': 'metapangenome-fastas.txt'})
