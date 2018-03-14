# -*- coding: utf-8
# pylint: disable=line-too-long
"""Helper functions for the anvi'o snakemake workflows"""

import os

import anvio
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Run, Progress
from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"
__status__ = "Development"


# Mock progress object that will not report anything, for general clarity.
progress = Progress()
progress.verbose = False

run = Run()
run.verbose = False


class WorkflowSuperClass:
    def __init__(self, config):
        self.config = config
        self.rules = []
        self.rule_acceptable_params_dict = {}
        self.dirs_dict = {}
        self.general_params = []


    def init(self):
        for rule in self.rules:
            if rule not in self.rule_acceptable_params_dict:
                self.rule_acceptable_params_dict[rule] = []

            params_that_all_rules_must_accept = ['threads']
            for param in params_that_all_rules_must_accept:
                if param not in self.rule_acceptable_params_dict[rule]:
                    self.rule_acceptable_params_dict[rule].append(param)

            general_params_that_all_workflows_must_accept = ['output_dirs']
            for param in general_params_that_all_workflows_must_accept:
                if param not in self.general_params:
                    self.general_params.append(param)


        self.dirs_dict = get_dir_names(self.config)

        # make sure that config file doesn't have garbage
        self.check_config()

        # create log dir if it doesn't exist
        os.makedirs(dirs_dict["LOGS_DIR"], exist_ok=True)


    def check_config(self):
        acceptable_params = set(self.rules + self.general_params)
        wrong_params = [p for p in self.config if p not in acceptable_params]
        if wrong_params:
            raise ConfigError("some of the parameters in your config file are not familiar to us. \
                        Here is a list of the wrong parameters: %s. This workflow only accepts \
                        the following general parameters: %s. And these are the rules in this \
                        workflow: %s." % (wrong_params, self.general_params, self.rules))

        self.check_rule_params()


    def check_rule_params(self):
        for rule in self.rules:
            if rule in self.config:
                wrong_params = [p for p in self.config[rule] if p not in self.rule_acceptable_params_dict[rule]]
                if wrong_params:
                    raise ConfigError("some of the parameters in your config file for rule %s are not familiar to us. \
                                Here is a list of the wrong parameters: %s. The only acceptable \
                                parameters for this rule are %s." % (rule, wrong_params, self.rule_acceptable_params_dict[rule]))


    def save_empty_config_in_json_format(self, filename='empty_config.json'):
        import json
        filesnpaths.is_output_file_writable(filename)

        empty_config = self.get_empty_config()

        open(filename, 'w').write(json.dumps(empty_config, indent=4))


    def get_empty_config(self):
        ''' This returns a dictionary with all the possible configurables for a workflow'''

        empty_config = {}

        for rule in self.rules:
            empty_config[rule] = {}
            for param in self.rule_acceptable_params_dict[rule]:
                empty_config[rule][param] = ''

        for param in self.general_params:
            empty_config[param] = ''

        return empty_config


# The config file contains many essential configurations for the workflow
# Setting the names of all directories
dirs_dict = {"LOGS_DIR"     : "00_LOGS"         ,\
             "QC_DIR"       : "01_QC"           ,\
             "ASSEMBLY_DIR" : "02_ASSEMBLY"     ,\
             "CONTIGS_DIR"  : "03_CONTIGS"      ,\
             "MAPPING_DIR"  : "04_MAPPING"      ,\
             "PROFILE_DIR"  : "05_ANVIO_PROFILE",\
             "MERGE_DIR"    : "06_MERGED"       ,\
             "PAN_DIR"      : "07_PAN"          ,\
             "FASTA_DIR"    : "01_FASTA"        ,\
             "LOCI_DIR"     : "04_LOCI_FASTAS"   \
}


########################################
# Helper functions
########################################

def A(_list, d, default_value = ""):
    '''
        A helper function to make sense of config details.
        string_list is a list of strings (or a single string)
        d is a dictionary

        this function checks if the strings in x are nested values in y.
        For example if x = ['a','b','c'] then this function checkes if the
        value y['a']['b']['c'] exists, if it does then it is returned
    '''
    if type(_list) is not list:
        # converting to list for the cases of only one item
        _list = [_list]
    while _list:
        a = _list.pop(0)
        if a in d:
            d = d[a]
        else:
            return default_value
    return d


def B(config, _rule, _param, default=''):
    # helper function for params
    val = A([_rule, _param], config, default)
    if val:
        if isinstance(val, bool):
            # the param is a flag so no need for a value
            val = ''
        return _param + ' ' + val
    else:
        return ''


# a helper function to get the user defined number of threads for a rule
def T(config, rule_name, N=1): return A([rule_name,'threads'], config, default_value=N)


def get_dir_names(config):
    ########################################
    # Reading some definitions from config files (also some sanity checks)
    ########################################
    DICT = dirs_dict
    for d in A("output_dirs", config):
        # renaming folders according to the config file, if the user specified.
        if d not in DICT:
            # making sure the user is asking to rename an existing folder.
            raise ConfigError("You define a name for the directory '%s' in your "\
                              "config file, but the only available folders are: "\
                              "%s" % (d, DICT))

        DICT[d] = A(d,config["output_dirs"])
    return DICT

def get_path_to_workflows_dir():
    # this returns a path
    base_path = os.path.dirname(__file__)
    return base_path


def warning_for_param(config, rule, param, wildcard, our_default=None):
    value = A([rule, param], config)
    if value:
        warning_message = 'You chose to define %s for the rule %s in the config file as %s.\
                           while this is allowed, know that you are doing so at your own risk.\
                           The reason this is risky is because this rule uses a wildcard/wildcards\
                           and hence is probably running more than once, and this might be cause a problem.\
                           In case you wanted to know, these are the wildcards used by this rule: %s' % (param, rule, value, wildcard)
        if our_default:
            warning_message = warning_message + ' Just so you are aware, if you dont provide a value\
                                                 in the config file, the default value is %s' % wildcard
        run.warning(warning_message)
