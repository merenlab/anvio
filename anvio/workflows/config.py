"""Configuration helper functions for anvi'o workflows."""

import json

import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


run = terminal.Run()


# The config file contains many essential configurations for the workflow
# Setting the names of all directories
dirs_dict = {"LOGS_DIR"     : "00_LOGS"         ,
             "QC_DIR"       : "01_QC"           ,
             "FASTA_DIR"    : "02_FASTA"     ,
             "CONTIGS_DIR"  : "03_CONTIGS"      ,
             "MAPPING_DIR"  : "04_MAPPING"      ,
             "PROFILE_DIR"  : "05_ANVIO_PROFILE",
             "MERGE_DIR"    : "06_MERGED"       ,
             "PAN_DIR"      : "07_PAN"          ,
             "LOCI_DIR"     : "04_LOCI_FASTAS"
}


def A(_list, d, default_value=""):
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
        try:
            d = d[a]
        except:
            return default_value
    return d


def B(config, _rule, _param, default=''):
    """Return a formatted CLI argument for a rule parameter from a config dictionary."""
    val = A([_rule, _param], config, default)
    if val:
        if isinstance(val, bool):
            # the param is a flag so no need for a value
            val = ''
        return _param + ' ' + str(val)
    else:
        return ''


def T(config, rule_name, N=1):
    """Return the configured thread count for a rule from a config dictionary."""
    return A([rule_name,'threads'], config, default_value=N)


def get_workflow_name_and_version_from_config(config_file, dont_raise=False):
    """Return the workflow name and config version from a workflow config file."""
    filesnpaths.is_file_json_formatted(config_file)
    config = json.load(open(config_file))
    workflow_name = config.get('workflow_name')
    # Notice that if there is no config_version then we return "0".
    # This is in order to accomodate early config files that had no such parameter.
    version = config.get('config_version', "0")

    if (not dont_raise) and (not workflow_name):
        raise ConfigError('Config files must contain a workflow_name.')

    return (workflow_name, version)


def get_dir_names(config, dont_raise=False):
    """Return workflow output directory names after applying config overrides."""
    ########################################
    # Reading some definitions from config files (also some sanity checks)
    ########################################
    DICT = dirs_dict
    for d in A("output_dirs", config):
        # renaming folders according to the config file, if the user specified.
        if d not in DICT and not dont_raise:
            # making sure the user is asking to rename an existing folder.
            raise ConfigError("You define a name for the directory '%s' in your "
                              "config file, but the only available folders are: "
                              "%s" % (d, DICT))

        DICT[d] = A(d,config["output_dirs"])
    return DICT


def check_for_risky_param_change(config, rule, param, wildcard, our_default=None):
    """Warn when a wildcard-sensitive rule parameter differs from the default."""
    value = A([rule, param], config)
    if value != our_default:
        warning_message = 'You chose to define %s for the rule %s in the config file as %s.\
                           while this is allowed, know that you are doing so at your own risk.\
                           The reason this is risky is because this rule uses a wildcard/wildcards\
                           and hence is probably running more than once, and this might cause a problem.\
                           In case you wanted to know, these are the wildcards used by this rule: %s' % (param, rule, value, wildcard)
        if our_default:
            warning_message = warning_message + ' Just so you are aware, if you dont provide a value\
                                                 in the config file, the default value is %s' % wildcard
        run.warning(warning_message)


def get_fields_for_fasta_information():
    """ Return a list of legitimate column names for fasta.txt files"""
    # Notice we don't include the name of the first column because
    # utils.get_TAB_delimited_file_as_dictionary doesn't really care about it.
    return ["path", "external_gene_calls", "gene_functional_annotation"]
