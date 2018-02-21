# -*- coding: utf-8
# pylint: disable=line-too-long
"""Helper functions for the anvi'o snakemake workflows"""

import os
import sys

import anvio

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

# create log dir if it doesn't exist
os.makedirs(dirs_dict["LOGS_DIR"], exist_ok=True)

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
