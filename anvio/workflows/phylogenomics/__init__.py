# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o phylogenomics workflows.
"""


import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.workflows import WorkflowSuperClass
from anvio.workflows.contigs import ContigsDBWorkflow
from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"

run = terminal.Run()
progress = terminal.Progress()

class PhylogenomicsWorkflow(ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='phylogenomics')

        # initialize the base class
        ContigsDBWorkflow.__init__(self)

        self.input_for_anvi_get_sequences_for_hmm_hits = {}
        self.internal_genomes_file = ''
        self.external_genomes_file = ''

        # initialize the base class

        self.rules.extend(['anvi_get_sequences_for_hmm_hits', 'trimal', 'iqtree'])

        self.general_params.extend(['project_name', 'internal_genomes', 'external_genomes'])

        self.dirs_dict.update({"PHYLO_DIR": "01_PHYLOGENOMICS"})

        self.default_config.update({'anvi_get_sequences_for_hmm_hits': {'--return-best-hit': True,
                                                                        '--align-with': 'famsa',
                                                                        '--concatenate-genes': True,
                                                                        '--get-aa-sequences': True,
                                                                        '--hmm-sources': 'Campbell_et_al'},
                                    'trimal': {'-gt': 0.5},
                                    'iqtree': {'threads': 8, '-m': 'WAG', '-bb': 1000}})

        get_sequences_params = ['--return-best-hit', \
                                '--separator', '--align-with', '--min-num-bins-gene-occurs', \
                                '--max-num-genes-missing-from-bin', '--concatenate-genes', \
                                '--get-aa-sequences', '--gene-names', '--hmm-sources']
        self.rule_acceptable_params_dict['anvi_get_sequences_for_hmm_hits'] = get_sequences_params
        self.rule_acceptable_params_dict['trimal'] = ['-gt', 'additional_params']
        self.rule_acceptable_params_dict['iqtree'] = ['-m', '-bb', 'additional_params']


    def init(self):
        ''' backhand stuff (mostly sanity checks) specific for the phylogenomics workflow'''
        super().init()

        self.internal_genomes_file = self.get_param_value_from_config('internal_genomes')
        self.external_genomes_file = self.get_param_value_from_config('external_genomes')
        self.input_for_anvi_get_sequences_for_hmm_hits = self.get_internal_and_external_genomes_files()

        self.sanity_checks()


    def sanity_checks(self):
        if not self.get_rule_param('anvi_get_sequences_for_hmm_hits', '--return-best-hit'):
            run.warning('You changed the value for --return-best-hit for the rule anvi_get_sequences_for_hmm_hits \
                         to something other than the default value, which is "true", while we allow you to do it \
                         this is likely to break things, we trust that you know what you are doing, but advise you \
                         to proceed with caution.')
