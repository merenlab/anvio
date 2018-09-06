# -*- coding: utf-8
# pylint: disable=line-too-long
"""Helper functions for the anvi'o snakemake workflows"""

import os
import sys
import json
import copy
import snakemake

import anvio
import anvio.utils as u
import anvio.errors as errors
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
r = errors.remove_spaces


class WorkflowSuperClass:
    def __init__(self):
        if 'args' not in self.__dict__:
            raise ConfigError("You need to initialize `WorkflowSuperClass` from within a class that\
                               has a member `self.args`.")

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        # FIXME: it is redundant to have both config and config_file
        # Alon and Meren will discuss how to do this right.
        # but for now this is how it is so things would work.
        # basically, when this class is called from a snakefile
        # then a config (dictionary) will be provided.
        # When this class is called from anvi-run-snakemake-workflow
        # for sanity checks  and things like that, then a config file
        # will be provided.
        self.config = A('config')
        self.config_file = A('config_file')
        self.default_config_output_path = A('get_default_config')
        self.save_workflow_graph = A('save_workflow_graph')
        self.list_dependencies = A('list_dependencies')
        self.dry_run_only = A('dry_run')
        self.additional_params = A('additional_params')

        if self.additional_params:
            run.warning("OK, SO THIS IS SERIOUS, AND WHEN THINGS ARE SERIOUS THEN WE USE CAPS. \
                         WE SEE THAT YOU ARE USING --additional-params AND THAT'S GREAT, BUT WE \
                         WANT TO REMIND YOU THAT ANYTHING THAT FOLLOWS --additional-params WILL \
                         BE CONSIDERED AS A snakemake PARAM THAT IS TRANSFERRED TO snakemake DIRECTLY. \
                         So make sure that these don't include anything that you didn't mean to \
                         include as an additional param: %s." % ', '.join(str(i) for i in self.additional_params))

        if self.save_workflow_graph:
            self.dry_run_only = True

        # if this class is being inherited from a snakefile that was 'included' from
        # within another snakefile.
        self.slave_mode = A('slave_mode')

        if self.config_file:
            filesnpaths.is_file_json_formatted(self.config_file)
            self.config = json.load(open(self.config_file))

        self.rules = []
        self.rule_acceptable_params_dict = {}
        self.dirs_dict = {}
        self.general_params = []
        self.default_config = {}
        self.rules_dependencies = {}
        self.forbidden_params = {}


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

        self.dirs_dict.update({"LOGS_DIR": "00_LOGS"})

        self.default_config = self.get_default_config()

        # the user requested to get a default config path for the workflow
        if self.default_config_output_path:
            self.save_default_config_in_json_format(self.default_config_output_path)
            sys.exit()
        elif not self.config:
            raise ConfigError("You need a config file to run this :/ If you need help to start preparing\
                               a config file for the anvi'o %s workflow, you can try the `--get-default-config`\
                               flag." % (self.name))

        self.dirs_dict.update(self.config.get("output_dirs", ''))

        # create log dir if it doesn't exist
        os.makedirs(self.dirs_dict["LOGS_DIR"], exist_ok=True)

        # lets check everything
        if not self.slave_mode:
            self.check_config()
            self.check_rule_params()


    def go(self, skip_dry_run=False):
        """Do the actual running"""

        if self.save_workflow_graph or (not skip_dry_run):
            self.dry_run()

        if self.dry_run_only:
            return

        # snakemake.main() accepts an `argv` parameter, but then the code has mixed responses to
        # that, and at places continues to read from sys.argv in a hardcoded manner. so we have to
        # overwrite our argv here.
        original_sys_argv = copy.deepcopy(sys.argv)

        sys.argv = ['snakemake',
                    '--snakefile',
                    get_workflow_snake_file_path(self.args.workflow),
                    '--configfile',
                    self.args.config_file]

        if self.additional_params:
            sys.argv.extend(self.additional_params)

        if self.list_dependencies:
            sys.argv.extend(['--dryrun', '--printshellcmds'])
            snakemake.main()
            sys.exit(0)
        else:
            sys.argv.extend(['-p'])
            snakemake.main()

        # restore the `sys.argv` to the original for the sake of sakity (totally made up word,
        # but you already know what it measn. you're welcome.)
        sys.argv = original_sys_argv


    def dry_run(self, workflow_graph_output_file_path_prefix='workflow'):
        """Not your regular dry run.

           The purpose of this function is to make sure there is a way to check for
           workflow program dependencies before the workflow is actually run. this way,
           if there is a `check_workflow_program_dependencies` call at the end of the
           snake file `get_workflow_snake_file_path(self.name)`, it can be called with
           a compiled snakemake `workflow` instance."""

        if self.slave_mode:
            return

        self.progress.new('Bleep bloop')
        self.progress.update('Quick dry run for an initial sanity check ...')
        args = ['snakemake', '--snakefile', get_workflow_snake_file_path(self.name), \
                '--configfile', self.config_file, '--dryrun', '--quiet']

        if self.save_workflow_graph:
            args.extend(['--dag'])

        log_file_path = filesnpaths.get_temp_file_path()
        u.run_command(args, log_file_path)
        self.progress.end()

        # here we're getting the graph info from the log file like a dirty hacker
        # we are (it still may be better to do it elsewhere more appropriate .. so
        # we can look more decent or whatever):
        if self.save_workflow_graph:
            lines = open(log_file_path, 'rU').readlines()

            try:
                line_of_interest = [line_no for line_no in range(0, len(lines)) if lines[line_no].startswith('digraph')][0]
            except IndexError:
                raise ConfigError("Oh no. Anvi'o was trying to generate a DAG output for you, but something must have\
                                   gone wrong in a step prior. Something tells anvi'o that if you take a look at the\
                                   log file here, you may be able to figure it out: '%s'. Sorry!" % log_file_path)
            open(workflow_graph_output_file_path_prefix + '.dot', 'w').write(''.join(lines[line_of_interest:]))

            self.run.info('Workflow DOT file', workflow_graph_output_file_path_prefix + '.dot')

            if u.is_program_exists('dot', dont_raise=True):
                dot_log_file = filesnpaths.get_temp_file_path()
                u.run_command(['dot', '-Tpng', workflow_graph_output_file_path_prefix + '.dot', '-o', workflow_graph_output_file_path_prefix + '.png'], dot_log_file)
                os.remove(dot_log_file)
                self.run.info('Workflow PNG file', workflow_graph_output_file_path_prefix + '.png')
            else:
                self.run.warning("Well, anvi'o was going to try to save a nice PNG file for your workflow\
                                  graph, but clearly you don't have `dot` installed on your system. That's OK. You\
                                  have your dot file now, and you can Google 'how to view dot file on [your operating\
                                  system goes here]', and install necessary programs (like .. `dot`).")

        os.remove(log_file_path)


    def check_workflow_program_dependencies(self, snakemake_workflow_object, dont_raise=False):
        """This function gets a snakemake workflow object and checks whether each shell command
           exists in the path.
        """

        if self.slave_mode:
            return

        shell_programs_needed = [r.shellcmd.strip().split()[0] for r in snakemake_workflow_object.rules if r.shellcmd]

        shell_programs_missing = [s for s in shell_programs_needed if not u.is_program_exists(s)]

        run.warning(None, 'Shell programs for the workflow')
        run.info('Needed', ', '.join(shell_programs_needed))
        run.info('Missing', ', '.join(shell_programs_missing) or 'None', nl_after=1)

        if len(shell_programs_missing):
            if dont_raise:
                return
            else:
                raise ConfigError("This workflow will not run without those missing programs are no longer\
                                   missing :(")

    def check_config(self):
        acceptable_params = set(self.rules + self.general_params)
        wrong_params = [p for p in self.config if p not in acceptable_params]
        if wrong_params:
            raise ConfigError("some of the parameters in your config file are not familiar to us. \
                        Here is a list of the wrong parameters: %s. This workflow only accepts \
                        the following general parameters: %s. And these are the rules in this \
                        workflow: %s." % (wrong_params, self.general_params, self.rules))

        wrong_dir_names = [d for d in self.config.get("output_dirs", '') if d not in self.dirs_dict]
        if wrong_dir_names:
            raise ConfigError("some of the directory names in your config file are not familiar to us. \
                        Here is a list of the wrong directories: %s. This workflow only has \
                        the following directories: %s." % (" ".join(wrong_dir_names), " ".join(list(self.dirs_dict.keys()))))


    def get_default_config(self):
        c = self.fill_empty_config_params(self.default_config)
        c["output_dirs"] = self.dirs_dict
        return c


    def check_additional_params(self, rule):
        ''' Check if the user is trying to use additional_params to set a param that is hard coded'''
        params = []
        if 'additional_params' in self.config[rule].keys() and self.forbidden_params.get(rule):
            # if the rule has 'additional_params' we need to make sure
            # that the user didn't include forbidden params there as well
            params = self.config[rule]['additional_params'].split(' ')

            bad_params = [p for p in self.forbidden_params.get(rule) if p in params]
            if bad_params:
                raise ConfigError("You are not allowed to set the following parameter/s: \
                                   %s for rule %s. These parameters are hard-coded. If you \
                                   are confused or upset please refer to an anvi'o developer \
                                   or a friend for support." % (', '.join(bad_params), rule))


    def check_rule_params(self):
        for rule in self.rules:
            if rule in self.config:
                wrong_params = [p for p in self.config[rule] if p not in self.rule_acceptable_params_dict[rule]]
                if wrong_params:
                    raise ConfigError("some of the parameters in your config file for rule %s are not familiar to us. \
                                Here is a list of the wrong parameters: %s. The only acceptable \
                                parameters for this rule are %s." % (rule, wrong_params, self.rule_acceptable_params_dict[rule]))

                self.check_additional_params(rule)


    def save_empty_config_in_json_format(self, file_path='empty_config.json'):
        self.save_config_in_json_format(file_path, self.get_empty_config())
        self.run.info("Empty config file", "Stored for workflow '%s' as '%s'." % (self.name, file_path))


    def save_default_config_in_json_format(self, file_path='default_config.json'):
        self.save_config_in_json_format(file_path, self.default_config)
        self.run.info("Default config file", "Stored for workflow '%s' as '%s'." % (self.name, file_path))


    def save_config_in_json_format(self, file_path, config):
        filesnpaths.is_output_file_writable(file_path)
        open(file_path, 'w').write(json.dumps(config, indent=4))


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


    def set_config_param(self, _list, value):
        d = self.config
        if type(_list) is not list:
            # converting to list for the cases of only one item
            _list = [_list]
        while _list:
            a = _list.pop(0)
            if a not in d:
                d[a] = {}
            f = d
            d = d[a]
            print(d)
        f[a] = value


    def get_param_value_from_config(self, _list, repress_default=False):
        '''
            A helper function to make sense of config details.
            string_list is a list of strings (or a single string)

            this function checks if the strings in x are nested values in self.config.
            For example if x = ['a','b','c'] then this function checkes if the
            value self.config['a']['b']['c'] exists, if it does then it is returned

            repress_default - If there is a default defined for the parameter (it would be defined
            under self.default_config), and the user didn't supply a parameter
            then the default will be returned. If this flad (repress_default) is set to True
            then this behaviour is repressed and instead an empty string would be returned.

        '''
        d = self.config
        default_dict = self.default_config
        if type(_list) is not list:
            # converting to list for the cases of only one item
            _list = [_list]
        while _list:
            a = _list.pop(0)
            default_dict = default_dict[a]
            try:
                d = d.get(a, None)
            except:
                # we continue becuase we want to get the value from the default config
                continue

        if (d is not None) or repress_default:
            return d
        else:
            return default_dict


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
            return _param + ' ' + str(val)
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
        return _param + ' ' + str(val)
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


def get_workflow_snake_file_path(workflow):
    workflow_dir = os.path.join(get_path_to_workflows_dir(), workflow)

    if not os.path.isdir(workflow_dir):
        raise ConfigError("Anvi'o does not know about the workflow '%s' :/")

    snakefile_path = os.path.join(workflow_dir, 'Snakefile')

    if not os.path.exists(snakefile_path):
        raise ConfigError("The snakefile path for the workflow '%s' seems to be missing :/" % workflow)

    return snakefile_path


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
