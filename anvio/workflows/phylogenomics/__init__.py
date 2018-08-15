# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o phylogenomics workflows.
"""


import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.workflows import WorkflowSuperClass
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

        self.input_for_anvi_get_sequences_for_hmm_hits = {}
        self.internal_genomes_file = ''
        self.external_genomes_file = ''

        # initialize the base class
        WorkflowSuperClass.__init__(self)

        self.rules.extend(['anvi_get_sequences_for_hmm_hits', 'trimal', 'iqtree'])

        self.general_params.extend(['project_name'])

        self.dirs_dict.update({"PHYLO_DIR": "01_PHYLOGENOMICS"})

        self.default_config.update({'anvi_get_sequences_for_hmm_hits': {'--return-best-hit': True,
                                                                        '--align-with': 'famsa',
                                                                        '--concatenate-genes': True,
                                                                        '--get-aa-sequences': True,
                                                                        '--hmm-sources': 'Campbell_et_al'},
                                    'trimal': {'-gt': 0.5},
                                    'iqtree': {'threads': 8, '-m': 'WAG', '-bb': 1000}})

        get_sequences_params = ['--external-genomes', '--internal-genomes', '--return-best-hit', \
                                '--separator', '--align-with', '--min-num-bins-gene-occurs', \
                                '--max-num-genes-missing-from-bin', '--concatenate-genes', \
                                '--get-aa-sequences', '--gene-names', '--hmm-sources']
        self.rule_acceptable_params_dict['anvi_get_sequences_for_hmm_hits'] = get_sequences_params
        self.rule_acceptable_params_dict['trimal'] = ['-gt', 'additional_params']
        self.rule_acceptable_params_dict['iqtree'] = ['-m', '-bb', 'additional_params']


    def init(self):
        ''' backhand stuff (mostly sanity checks) specific for the phylogenomics workflow'''
        super().init()

        internal_genomes_file = self.get_param_value_from_config(['anvi_get_sequences_for_hmm_hits', '--internal-genomes'])
        external_genomes_file = self.get_param_value_from_config(['anvi_get_sequences_for_hmm_hits', '--external-genomes'])

        if not internal_genomes_file and not external_genomes_file:
            raise ConfigError('You must provide either an external genomes file or internal genomes file \
                               for the rule anvi_get_sequences_for_hmm_hits')

        # here we do a little trick to make sure the rule can expect either one or both
        self.input_for_anvi_get_sequences_for_hmm_hits = {"internal_genomes_file": external_genomes_file,
                                                          "external_genomes_file": internal_genomes_file}

        if internal_genomes_file:
            filesnpaths.is_file_exists(internal_genomes_file)
            self.input_for_anvi_get_sequences_for_hmm_hits['internal_genomes_file'] = internal_genomes_file
            self.internal_genomes_file = internal_genomes_file

        if external_genomes_file:
            filesnpaths.is_file_exists(external_genomes_file)
            self.input_for_anvi_get_sequences_for_hmm_hits['external_genomes_file'] = external_genomes_file
            self.external_genomes_file = external_genomes_file

        if not self.get_rule_param('anvi_get_sequences_for_hmm_hits', '--return-best-hit'):
            run.warning('You changed the value for --return-best-hit for the rule anvi_get_sequences_for_hmm_hits \
                         to something other than the default value, which is "true", while we allow you to do it \
                         this is likely to break things, we trust that you know what you are doing, but advise you \
                         to proceed with caution.')
