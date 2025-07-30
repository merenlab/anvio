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
from anvio.version import versions_for_db_types


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
r = errors.remove_spaces

workflow_config_version = versions_for_db_types['config']

class WorkflowSuperClass:
    def __init__(self):
        if 'args' not in self.__dict__:
            raise ConfigError("You need to initialize `WorkflowSuperClass` from within a class that "
                              "has a member `self.args`.")

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
        self.skip_version_check = A('skip_version_check')

        if self.additional_params:
            run.warning("OK, SO THIS IS SERIOUS, AND WHEN THINGS ARE SERIOUS THEN WE USE CAPS. "
                        "WE SEE THAT YOU ARE USING --additional-params AND THAT'S GREAT, BUT WE "
                        "WANT TO REMIND YOU THAT ANYTHING THAT FOLLOWS --additional-params WILL "
                        "BE CONSIDERED AS A snakemake PARAM THAT IS TRANSFERRED TO snakemake DIRECTLY. "
                        "So make sure that these don't include anything that you didn't mean to "
                        "include as an additional param: %s." % ', '.join(str(i) for i in self.additional_params))

        if self.save_workflow_graph:
            self.dry_run_only = True

        # if this class is being inherited from a snakefile that was 'included' from
        # within another snakefile.
        self.this_workflow_is_inherited_by_another = A('this_workflow_is_inherited_by_another')

        if self.config_file:
            filesnpaths.is_file_json_formatted(self.config_file)
            self.config = json.load(open(self.config_file))

        if self.config and ('config_version' not in self.config or 'workflow_name' not in self.config):
            raise ConfigError(f"You config file '{self.config_file}' is probably too old as it is missing some key varaibles "
                              f"in it. Please make sure your config includes entries for both `config_version` and "
                              f"`workflow_name`.")

        # If a config exists, check that its current
        if self.config and not self.skip_version_check:
            if self.config['config_version'] != workflow_config_version:
                raise ConfigError(f"Anvi'o couldn't get things moving because the version of your config file is out "
                                  f"of date (your version: {self.config['config_version']}; up-to-date version: "
                                  f"{workflow_config_version}). Not a problem though, simply run `anvi-migrate {self.config_file} "
                                  f"--migrate-safely` and it will be updated. Then re-run the command producing "
                                  f"this error message.")

        self.rules = []
        self.rule_acceptable_params_dict = {}
        self.dirs_dict = {}
        self.general_params = []
        self.default_config = {}
        self.rules_dependencies = {}
        self.forbidden_params = {}


    def init(self):
        run.warning('We are initiating parameters for the %s workflow' % self.name)

        for rule in self.rules:
            if rule not in self.rule_acceptable_params_dict:
                self.rule_acceptable_params_dict[rule] = []

            params_that_all_rules_must_accept = self.get_params_that_all_rules_accept()
            for param in params_that_all_rules_must_accept:
                if param not in self.rule_acceptable_params_dict[rule]:
                    self.rule_acceptable_params_dict[rule].append(param)

            for param in self.get_global_general_params():
                if param not in self.general_params:
                    self.general_params.append(param)

        self.dirs_dict.update({"LOGS_DIR": "00_LOGS"})

        self.default_config = self.get_default_config()

        # the user requested to get a default config path for the workflow
        if self.default_config_output_path:
            self.save_default_config_in_json_format(self.default_config_output_path)
            sys.exit()
        elif not self.config:
            raise ConfigError("You need a config file to run this :/ If you need help to start preparing "
                              "a config file for the anvi'o %s workflow, you can try the `--get-default-config` "
                              "flag." % (self.name))

        self.dirs_dict.update(self.config.get("output_dirs", ''))

        # create log dir if it doesn't exist
        os.makedirs(self.dirs_dict["LOGS_DIR"], exist_ok=True)

        # lets check everything
        if not self.this_workflow_is_inherited_by_another:
            self.check_config()
            self.check_rule_params()


    def get_params_that_all_rules_accept(self):
        return ['threads']


    def get_global_general_params(self):
        ''' Return a list of the general parameters that are always acceptable.'''
        return ['output_dirs', 'max_threads', 'config_version', 'workflow_name']


    def sanity_checks(self):
        ''' each workflow has its own sanity checks but we only run these when we go'''
        # each workflow will overide this function with specific things to check
        pass


    def warn_user_regarding_param_with_wildcard_default_value(self, rule_name, param, wildcard_name):
        try:
            default_value = self.default_config[rule_name][param]
        except KeyError:
            raise ConfigError('Someone is trying to read default values for parameters that '
                              'dont have default values. These are the offending rule names and '
                              'parameter: %s, %s' % (rule_name, param))
        check_for_risky_param_change(self.config, rule_name, param, wildcard_name, default_value)


    def get_max_num_cpus_requested_by_the_workflow(self):
        """A semi-smart heuristic to find out what we shold set for `--cores`.

        If the user has defined an additional parameter '--cores', or '-c' in their `anvi-run-workflow`
        command, then we are all good: that is our max threads value. The only thing this function in that
        case would be to make sure that value exceeds the maximum number of threads asked for by any rule.

        If there are no `--cores` or `-c` parameter among additional params, then we have to figure
        something out and send it to snakemake. In that case this function will check if the user set a
        global `max_threads` pameter (and that it exceeds the number of threads asked for by any rule),
        and send that value to snakemake.

        If there are no `--cores` or `-c` parameter among additional params, and there is no value set
        for the global `max_threads` in the config file, then this function will find the max num threads
        requested by any rule, and will send it to snakeamke as the max num CPUs.
        """

        # figure out whether the user sent `--cores` or `-c` as additional params
        if self.additional_params and '--cores' in self.additional_params:
            num_cores_additional_param = self.additional_params[self.additional_params.index('--cores') + 1]
        elif self.additional_params and '-c' in self.additional_params:
            num_cores_additional_param = self.additional_params[self.additional_params.index('-c') + 1]
        else:
            num_cores_additional_param = None

        if num_cores_additional_param:
            try:
                num_cores_additional_param = int(num_cores_additional_param)
            except:
                raise ConfigError(f"So anvi'o was trying to make sense of the `--cores` value you set via the "
                                  f"`--additional-params` mechanism in your command, but this is what it found "
                                  f"there instead of a proper integer: '{num_cores_additional_param}' :( Anvi'o "
                                  f"is not sure if it is you who is screwing something up, or it is anvi'o :/")

        # figure out the maximum number of threads requested by any single job
        max_num_threads_set_by_any_rule = 0
        rule_name_for_max_num_threads_value = None
        for rule in self.rules:
            threads = A([rule, 'threads'], self.config)

            if not threads:
                continue
            else:
                threads = int(threads)

            if threads > max_num_threads_set_by_any_rule:
                max_num_threads_set_by_any_rule = threads
                rule_name_for_max_num_threads_value = rule

        # if the user set `--cores`, make sure that value is larger than the max_num_threads_set_by_any_rule
        if num_cores_additional_param:
            if num_cores_additional_param < max_num_threads_set_by_any_rule:
                raise ConfigError(f"Funny story here. You've set the max number of CPUs to be used for your workflow to "
                                  f"{num_cores_additional_param} via the `--cores` additional parameter. But then, the "
                                  f"rule {rule_name_for_max_num_threads_value} is asking for {max_num_threads_set_by_any_rule}. "
                                  f"Das ist nicht gut. Please either decrease the number of threads for your rule, or "
                                  f"increase the number of CPUs you are asking from your system to make sure your job is "
                                  f"not going to get killed by a manager process later.")
            else:
                # the user set `--cores`, all good here. We will return nothing, and therefore we will not
                # add any additional parameter to our call.
                return None

        # if we are here, it means the user did not set a `--cores` parameter through the additional params
        # mechanism. we basically need to figure out the maximum number of threads set as a global parameter
        global_max_threads_value = None
        if 'max_threads' in self.config and self.config.get('max_threads'):
            global_max_threads_value = int(self.config.get('max_threads'))
        else:
            global_max_threads_value = None

        # now it is time to make a suggestion as our max threads value for snakemake so it can be set as `--cores` argument.`
        if not global_max_threads_value:
            if max_num_threads_set_by_any_rule == 0:
                raise ConfigError("You have not set a global `max_threads` value in your config and none of "
                                  "your rules includes a value for their `threads`. In this case it is impossible "
                                  "to determine the maximum number of CPU cores we should ask snakemake to use to run your "
                                  "workflow :/ Please either define a `max_threads` value in your config file, OR include "
                                  "a `--cores XX` parameter via the `--additional-params` mechanism of `anvi-run-workflow`, "
                                  "where XX is the maximum number of CPUs you wish to assign to this job.")
            else:
                self.run.warning(f"You haven't set a global `max_threads` value in your config file, and you haven't declared the "
                                 f"maximum number of CPUs you wish to assign to this job via `--cores`.  But anvi'o found out "
                                 f"that for the rule `{rule_name_for_max_num_threads_value}` you asked for {max_num_threads_set_by_any_rule} "
                                 f"threads , which is the largest number of threads among all rules. So, anvi'o will set up your job by "
                                 f"passing {max_num_threads_set_by_any_rule} as the maximum number of CPU cores to be used by your workflow. "
                                 f"If you'd like to assign more resources for your workflow, please either use the global `max_threads` "
                                 f"parameter in your config file, or include a `--cores XX` parameter via the `--additional-params` "
                                 f"mechanism, where XX is the maximum number of CPUs you wish to assign to this job.",
                                 header="MAX NUMBER OF CPU CORES HEURISTIC", lc="green")
                return max_num_threads_set_by_any_rule
        else:
            return global_max_threads_value


    def go(self, skip_dry_run=False):
        """Do the actual running"""

        max_num_cpus_requested_by_the_workflow = self.get_max_num_cpus_requested_by_the_workflow()

        self.sanity_checks()

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
            if max_num_cpus_requested_by_the_workflow:
                sys.argv.extend(['--dryrun', '--printshellcmds', '--cores', f'{max_num_cpus_requested_by_the_workflow}'])
            else:
                sys.argv.extend(['--dryrun', '--printshellcmds'])
            snakemake.main()
            sys.exit(0)
        else:
            if max_num_cpus_requested_by_the_workflow:
                sys.argv.extend(['-p', '--cores', f'{max_num_cpus_requested_by_the_workflow}'])
            else:
                sys.argv.extend(['-p'])
            snakemake.main()

        # restore the `sys.argv` to the original for the sake of sakity (totally made up word,
        # but you already know what it means. you're welcome.)
        sys.argv = original_sys_argv


    def dry_run(self, workflow_graph_output_file_path_prefix='workflow'):
        """Not your regular dry run.

           The purpose of this function is to make sure there is a way to check for
           workflow program dependencies before the workflow is actually run. this way,
           if there is a `check_workflow_program_dependencies` call at the end of the
           snake file `get_workflow_snake_file_path(self.name)`, it can be called with
           a compiled snakemake `workflow` instance."""

        if self.this_workflow_is_inherited_by_another:
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
            lines = open(log_file_path, 'r').readlines()

            try:
                line_of_interest = [line_no for line_no in range(0, len(lines)) if lines[line_no].startswith('digraph')][0]
            except IndexError:
                raise ConfigError("Oh no. Anvi'o was trying to generate a DAG output for you, but something must have "
                                  "gone wrong in a step prior. Something tells anvi'o that if you take a look at the "
                                  "log file here, you may be able to figure it out: '%s'. Sorry!" % log_file_path)
            open(workflow_graph_output_file_path_prefix + '.dot', 'w').write(''.join(lines[line_of_interest:]))

            self.run.info('Workflow DOT file', workflow_graph_output_file_path_prefix + '.dot')

            if u.is_program_exists('dot', dont_raise=True):
                dot_log_file = filesnpaths.get_temp_file_path()
                u.run_command(['dot', '-Tpdf', workflow_graph_output_file_path_prefix + '.dot', '-o', workflow_graph_output_file_path_prefix + '.pdf'], dot_log_file)
                os.remove(dot_log_file)
                self.run.info('Workflow PDF file', workflow_graph_output_file_path_prefix + '.pdf')
            else:
                self.run.warning("Well, anvi'o was going to try to save a nice PDF file for your workflow "
                                 "graph, but clearly you don't have `dot` installed on your system. If you "
                                 "install `graphviz` (i.e., by running \"mamba install -c conda-forge graphviz\" "
                                 "or \"conda install -c conda-forge graphviz\") you will no longer see this "
                                 "message, and instead get your workflow as a PDF file (alternatively you can "
                                 "always vuew the file `workflow.dot` using any dot viewer).")

        os.remove(log_file_path)


    def check_workflow_program_dependencies(self, snakemake_workflow_object, dont_raise=True):
        """Check whether each shell command in a snakemake_workflow_object exists in PATH

        Parameters
        ==========
        snakemake_workflow_object: snakemake.workflow
            Source code of this object found at
            https://snakemake.readthedocs.io/en/stable/_modules/snakemake/workflow.html

        Notes
        =====
        - FIXME Not all of the programs identified here will _actually_ be used in the workflow.
          Finding out which commands will actually be used requires building the DAG and then
          finding the appropriate place in the Snakemake API where we can expose this information.
          See https://github.com/merenlab/anvio/issues/1316 for discussion.
        """

        if self.this_workflow_is_inherited_by_another:
            return

        shell_programs_needed = [r.shellcmd.strip().split()[0] for r in snakemake_workflow_object.rules if r.shellcmd]

        shell_programs_missing = [s for s in shell_programs_needed if not u.is_program_exists(s, dont_raise=dont_raise)]

        run.warning(None, 'Shell programs for the workflow')
        run.info('Needed', ', '.join(shell_programs_needed))
        run.info('Missing', ', '.join(shell_programs_missing) or 'None', nl_after=1)

        if len(shell_programs_missing):
            if dont_raise:
                return
            else:
                raise ConfigError("This workflow will not run without those missing programs are no longer "
                                  "missing :(")

    def check_config(self):
        if not self.config.get('config_version'):
            raise ConfigError("Config files must include a config_version. If this is news to you, and/or you don't know what "
                              "version your config should be, please run in your terminal the command `anvi-migrate %s` to "
                              "upgrade your config file." % (self.config_file))

        if not self.config.get('workflow_name'):
            raise ConfigError('Config files must contain a workflow_name. You can simply add a line that goes like '
                              '"workflow_name": "metagenomics" wherever appropriate.')

        acceptable_params = set(self.rules + self.general_params)
        wrong_params = [p for p in self.config if p not in acceptable_params]
        if wrong_params:
            raise ConfigError("Some of the parameters in your config file are not familiar to us. "
                              "Here is a list of the wrong parameters: %s. This workflow only accepts "
                              "the following general parameters: %s. And these are the rules in this "
                              "workflow: %s." % (wrong_params, self.general_params, self.rules))

        wrong_dir_names = [d for d in self.config.get("output_dirs", '') if d not in self.dirs_dict]
        if wrong_dir_names:
            raise ConfigError("Some of the directory names in your config file are not familiar to us. "
                              "Here is a list of the wrong directories: %s. This workflow only has "
                              "the following directories: %s." % (" ".join(wrong_dir_names), " ".join(list(self.dirs_dict.keys()))))

        ## make sure max_threads is an integer number
        max_threads = self.get_param_value_from_config('max_threads')
        if max_threads:
            try:
                int(max_threads)
            except:
                raise ConfigError('"max_threads" must be an integer value. The value you provided in the config file is: "%s"' % max_threads)


    def get_default_config(self):
        c = self.fill_empty_config_params(self.default_config)
        c["output_dirs"] = self.dirs_dict
        c["config_version"] = workflow_config_version
        c["workflow_name"] = self.name
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
                raise ConfigError("You are not allowed to set the following parameter/s: "
                                  "%s for rule %s. These parameters are hard-coded. If you "
                                  "are confused or upset please refer to an anvi'o developer "
                                  "or a friend for support." % (', '.join(bad_params), rule))


    def check_rule_params(self):
        for rule in self.rules:
            if rule in self.config:
                wrong_params = [p for p in self.config[rule] if p not in self.rule_acceptable_params_dict[rule]]
                if wrong_params:
                    raise ConfigError("some of the parameters in your config file for rule %s are not familiar to us. "
                               "Here is a list of the wrong parameters: %s. The only acceptable "
                               "parameters for this rule are %s." % (rule, wrong_params, self.rule_acceptable_params_dict[rule]))

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
        f[a] = value


    def get_param_value_from_config(self, _list):
        '''
            A helper function to make sense of config details.
            string_list is a list of strings (or a single string)

            this function checks if the strings in x are nested values in self.config.
            For example if x = ['a','b','c'] then this function checkes if the
            value self.config['a']['b']['c'] exists, if it does then it is returned.
            If it does not exist then None is returned.
        '''
        d = self.config
        if type(_list) is not list:
            # converting to list for the cases of only one item
            _list = [_list]
        while _list:
            a = _list.pop(0)
            try:
                d = d.get(a, "")
            except:
                return ""

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
        if val is not None and val != '':
            if isinstance(val, bool):
                # the param is a flag so no need for a value
                if val:
                    return _param
            else:
                return _param + ' ' + str(val)
        return ''


    def T(self, rule_name):
        max_threads = self.get_param_value_from_config("max_threads")
        if not max_threads:
            max_threads = float("Inf")
        threads = self.get_param_value_from_config([rule_name,'threads'])
        if threads:
            try:
                if int(threads) > float(max_threads):
                    return int(max_threads)
                else:
                    return int(threads)
            except:
                raise ConfigError('"threads" must be an integer number. In your config file you provided "%s" for '
                                  'the number of threads for rule "%s"' % (threads, rule_name))
        else:
            return 1


    def init_workflow_super_class(self, args, workflow_name):
        '''
            if a regular instance of workflow object is being generated, we
            expect it to have a parameter `args`. if there is no `args` given, we
            assume the class is being inherited as a base class from within another.

            For a regular instance of a workflow this function will set the args
            and init the WorkflowSuperClass.
        '''
        if args:
            if len(self.__dict__):
                raise ConfigError("Something is wrong. You are inheriting %s from \
                                   within another class, yet you are providing an `args` parameter.\
                                   This is not alright." % type(self))
            self.args = args
            self.name = workflow_name
            WorkflowSuperClass.__init__(self)
            self.run = run
            self.progress = progress
        else:
            if not len(self.__dict__):
                raise ConfigError("When you are *not* inheriting %s from within "
                                  "a super class, you must provide an `args` parameter." % type(self))
            if 'name' not in self.__dict__:
                raise ConfigError("The super class trying to inherit %s does not "
                                  "have a set `self.name`. Which means there may be other things "
                                  "wrong with it, hence anvi'o refuses to continue." % type(self))


    def get_internal_and_external_genomes_files(self):
            internal_genomes_file = self.get_param_value_from_config('internal_genomes')
            external_genomes_file = self.get_param_value_from_config('external_genomes')

            if not internal_genomes_file and not external_genomes_file:
                raise ConfigError('You must provide either an external genomes file or internal genomes file')

            fasta_txt_file = self.get_param_value_from_config('fasta_txt')
            if fasta_txt_file:
                if not external_genomes_file:
                    raise ConfigError("You provided a fasta_txt, but didn't specify a path for an external-genomes file. "
                                      "If you wish to use external genomes, you must specify a name for the external-genomes "
                                      "file, using the `external_genomes` parameter in your config file. Just to clarify: "
                                      "the external genomes file DOESN'T HAVE TO EXIST. Anvi'o can create it for you by "
                                      "using the information you supplied in the `fasta_txt` file, but you still must specify "
                                      "a name for the external-genomes file. For example, you could use \"external_genomes\": \"external-genomes.txt\" "
                                      "(but feel free to be creative with the naming of your external-genomes file).")

                filesnpaths.is_file_tab_delimited(fasta_txt_file)

            # here we do a little trick to make sure the rule can expect either one or both
            d = {"internal_genomes_file": external_genomes_file,
                 "external_genomes_file": internal_genomes_file}

            if internal_genomes_file:
                filesnpaths.is_file_tab_delimited(internal_genomes_file)
                d['internal_genomes_file'] = internal_genomes_file

            if external_genomes_file:
                if not filesnpaths.is_file_exists(external_genomes_file, dont_raise=True):
                    run.warning('There is no file %s. No worries, one will be created for you.' % external_genomes_file)
                else:
                    filesnpaths.is_file_tab_delimited

                d['external_genomes_file'] = external_genomes_file

            return d


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
        try:
            d = d[a]
        except:
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


def D(debug_message, debug_log_file_path=".SNAKEMAKEDEBUG"):
    with open(debug_log_file_path, 'a') as output:
            output.write(terminal.get_date() + '\n')
            output.write(str(debug_message) + '\n\n')


# a helper function to get the user defined number of threads for a rule
def T(config, rule_name, N=1): return A([rule_name,'threads'], config, default_value=N)


def get_dir_names(config, dont_raise=False):
    ########################################
    # Reading some definitions from config files (also some sanity checks)
    ########################################
    DICT = dirs_dict
    for d in A("output_dirs", config):
        # renaming folders according to the config file, if the user specified.
        if d not in DICT and not dont_raise:
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
        raise ConfigError("Anvi'o does not know about the workflow '%s' :/" % workflow)

    snakefile_path = os.path.join(workflow_dir, 'Snakefile')

    if not os.path.exists(snakefile_path):
        raise ConfigError("The snakefile path for the workflow '%s' seems to be missing :/" % workflow)

    return snakefile_path


def check_for_risky_param_change(config, rule, param, wildcard, our_default=None):
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


def get_workflow_module_dict():
    from anvio.workflows.contigs import ContigsDBWorkflow
    from anvio.workflows.metagenomics import MetagenomicsWorkflow
    from anvio.workflows.pangenomics import PangenomicsWorkflow
    from anvio.workflows.phylogenomics import PhylogenomicsWorkflow
    from anvio.workflows.trnaseq import TRNASeqWorkflow
    from anvio.workflows.ecophylo import EcoPhyloWorkflow
    from anvio.workflows.sra_download import SRADownloadWorkflow

    workflows_dict = {'contigs': ContigsDBWorkflow,
                      'metagenomics': MetagenomicsWorkflow,
                      'pangenomics': PangenomicsWorkflow,
                      'phylogenomics': PhylogenomicsWorkflow,
                      'trnaseq': TRNASeqWorkflow,
                      'ecophylo': EcoPhyloWorkflow,
                      'sra_download': SRADownloadWorkflow}

    return workflows_dict


def get_workflow_name_and_version_from_config(config_file, dont_raise=False):
    filesnpaths.is_file_json_formatted(config_file)
    config = json.load(open(config_file))
    workflow_name = config.get('workflow_name')
    # Notice that if there is no config_version then we return "0".
    # This is in order to accomodate early config files that had no such parameter.
    version = config.get('config_version', "0")

    if (not dont_raise) and (not workflow_name):
        raise ConfigError('Config files must contain a workflow_name.')

    return (workflow_name, version)
