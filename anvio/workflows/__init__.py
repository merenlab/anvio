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

available_workflows = ['contigs']

class WorkflowSuperClass:
    def __init__(self, config):
        self.config = config
        self.rules = []
        self.rule_acceptable_params_dict = {}
        self.dirs_dict = {}
        self.general_params = []
        self.default_config = {}


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

        self.dirs_dict.update(
                {
                    "LOGS_DIR": "00_LOGS"
                }
                            )

        self.default_config = self.get_default_config()

        self.dirs_dict.update(self.config.get("output_dirs"))

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

        wrong_dir_names = [d for d in self.config.get("output_dirs") if d not in self.dirs_dict]
        if wrong_dir_names:
            raise ConfigError("some of the directory names in your config file are not familiar to us. \
                        Here is a list of the wrong directories: %s. This workflow only has \
                        the following directories: %s." % (" ".join(wrong_dir_names), " ".join(list(self.dirs_dict.keys()))))

        self.check_rule_params()


    def get_default_config(self):
        c = self.fill_empty_config_params(self.default_config)
        c.update(self.dirs_dict)
        return c


    def check_rule_params(self):
        for rule in self.rules:
            if rule in self.config:
                wrong_params = [p for p in self.config[rule] if p not in self.rule_acceptable_params_dict[rule]]
                if wrong_params:
                    raise ConfigError("some of the parameters in your config file for rule %s are not familiar to us. \
                                Here is a list of the wrong parameters: %s. The only acceptable \
                                parameters for this rule are %s." % (rule, wrong_params, self.rule_acceptable_params_dict[rule]))


    def save_empty_config_in_json_format(self, filename='empty_config.json'):
        self.save_config_in_json_format(filename, self.get_empty_config())


    def save_default_config_in_json_format(self, filename='default_config.json'):
        self.save_config_in_json_format(filename, self.default_config)


    def save_config_in_json_format(self, filename, config):
        import json
        filesnpaths.is_output_file_writable(filename)

        open(filename, 'w').write(json.dumps(config, indent=4))


    def get_empty_config(self):
        ''' This returns a dictionary with all the possible configurables for a workflow'''
        return self.fill_empty_config_params(config={})


    def fill_empty_config_params(self, config={}):
        ''' Takes a config dictionary and assigns an empty string to any parameter that wasnt defined in the config'''
        new_config = config.copy()

        for rule in self.rules:
            if rule not in config:
                new_config[rule] = {}
            for param in self.rule_acceptable_params_dict[rule]:
                new_config[rule][param] = new_config[rule].get(param, '')

        for param in self.general_params:
            if param not in config:
                new_config[param] = ''

        return new_config


    def get_param_value_from_config(self, _list):
        '''
            A helper function to make sense of config details.
            string_list is a list of strings (or a single string)

            this function checks if the strings in x are nested values in self.config.
            For example if x = ['a','b','c'] then this function checkes if the
            value self.config['a']['b']['c'] exists, if it does then it is returned
        '''
        d = self.config
        default_dict = self.default_config
        return_default = False
        if type(_list) is not list:
            # converting to list for the cases of only one item
            _list = [_list]
        while _list:
            a = _list.pop(0)
            print(default_dict)
            print(a)
            default_dict = default_dict[a]
            if a in d:
                d = d[a]
            else:
                return_default = True

        if return_default:
            return default_dict
        else:
            return d


    def get_rule_param(self, _rule, _param):
        '''
            returns the parameter as an input argument

            this function works differently for two kinds of parameters (flags, vs. non-flags)
            For example the parameter --use-ncbi-blast is a flag hence if the config file contains
            "--use-ncbi-blast": true, then this function would return "--use-ncbi-blast"

            for a non-flag the value of the parameter would also be included, for example,
            if the config file contains "--min-occurence: 5" then this function will return
            "--min-occurence 5"
        '''
        val = self.get_param_value_from_config([_rule, _param])
        if val:
            if isinstance(val, bool):
                # the param is a flag so no need for a value
                val = ''
            return _param + ' ' + val
        else:
            return ''


    def T(self, rule_name): return self.get_param_value_from_config([rule_name,'threads']) if self.get_param_value_from_config([rule_name,'threads']) else 1


# The config file contains many essential configurations for the workflow
# Setting the names of all directories
dirs_dict = {"LOGS_DIR"     : "00_LOGS"         ,\
             "QC_DIR"       : "01_QC"           ,\
             "FASTA_DIR"    : "02_FASTA"     ,\
             "CONTIGS_DIR"  : "03_CONTIGS"      ,\
             "MAPPING_DIR"  : "04_MAPPING"      ,\
             "PROFILE_DIR"  : "05_ANVIO_PROFILE",\
             "MERGE_DIR"    : "06_MERGED"       ,\
             "PAN_DIR"      : "07_PAN"          ,\
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
                           and hence is probably running more than once, and this might cause a problem.\
                           In case you wanted to know, these are the wildcards used by this rule: %s' % (param, rule, value, wildcard)
        if our_default:
            warning_message = warning_message + ' Just so you are aware, if you dont provide a value\
                                                 in the config file, the default value is %s' % wildcard
        run.warning(warning_message)
