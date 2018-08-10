# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o phylogenomics workflows.
"""


import anvio
import anvio.terminal as terminal

from anvio.workflows import WorkflowSuperClass
from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


class PhylogenomicsWorkflow(WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        # if a regular instance of `ContigsDBWorkflow` is being generated, we
        # expect it to have a parameter `args`. if there is no `args` given, we
        # assume the class is being inherited as a base class from within another
        if args:
            if len(self.__dict__):
                raise ConfigError("Something is wrong. You are ineriting `PhylogenomicsWorkflow` from \
                                   within another class, yet you are providing an `args` parameter.\
                                   This is not alright.")
            self.args = args
            self.name = 'phylogenomics'
        else:
            if not len(self.__dict__):
                raise ConfigError("When you are *not* inheriting `PhylogenomicsWorkflow` from within\
                                   a super class, you must provide an `args` parameter.")

            if 'name' not in self.__dict__:
                raise ConfigError("The super class trying to inherit `PhylogenomicsWorkflow` does not\
                                   have a set `self.name`. Which means there may be other things\
                                   wrong with it, hence anvi'o refuses to continue.")

        self.run = run
        self.progress = progress

        # initialize the base class
        WorkflowSuperClass.__init__(self)

        self.rules.extend([])

        self.general_params.extend([])

        self.dirs_dict.update({"FASTA_DIR": "01_FASTA",
                               "CONTIGS_DIR": "02_CONTIGS"})

        self.default_config.update({})

        self.rule_acceptable_params_dict[''] = []
