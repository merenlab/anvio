import os
import sys
import subprocess

import anvio
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.terminal import Run, Progress, get_date


def serialize_args(args, single_dash=False, use_underscore=False, skip_keys=None, translate=None):
    cmdline = []
    for param, value in args.__dict__.items():
        if isinstance(skip_keys, list):
            if param in skip_keys:
                continue

        if translate and param in translate:
            param = translate[param]

        dash = '-' if single_dash else '--'

        if not use_underscore:
            param = param.replace('_', '-')

        if value is True:
            cmdline.append('%s%s' % (dash, param))
        elif value is not False and value is not None:
            cmdline.append('%s%s' % (dash, param))
            cmdline.append(str(value))

    return cmdline



def format_cmdline(cmdline):
    """Takes a cmdline for `run_command` or `run_command_STDIN`, and makes it beautiful."""
    if not cmdline or (not isinstance(cmdline, str) and not isinstance(cmdline, list)):
        raise ConfigError("You made utils/commandline::format_cmdline upset. The parameter you sent to run kinda sucks. It should be string "
                           "or list type. Note that the parameter `shell` for subprocess.call in this `run_command` function "
                           "is always False, therefore if you send a string type, it will be split into a list prior to being "
                           "sent to subprocess.")

    if isinstance(cmdline, str):
        cmdline = [str(x) for x in cmdline.split(' ')]
    else:
        cmdline = [str(x) for x in cmdline]

    return cmdline



def run_command(cmdline, log_file_path, first_line_of_log_is_cmdline=True, remove_log_file_if_exists=True):
    """ Uses subprocess.call to run your `cmdline`

    Parameters
    ==========
    cmdline : str or list
        The command to be run, e.g. "echo hello" or ["echo", "hello"]
    log_file_path : str or Path-like
        All stdout from the command is sent to this filepath

    Raises ConfigError if ret_val < 0, or on OSError.  Does NOT raise if program terminated with exit code > 0.
    """
    cmdline = format_cmdline(cmdline)

    if anvio.DEBUG:
        Progress().reset()
        Run().info("[DEBUG] `run_command` is running", \
                   ' '.join(['%s' % (('"%s"' % str(x)) if ' ' in str(x) else ('%s' % str(x))) for x in cmdline]), \
                   nl_before=1, nl_after=1, mc='red', lc='yellow')

    filesnpaths.is_output_file_writable(log_file_path)

    if remove_log_file_if_exists and os.path.exists(log_file_path):
        os.remove(log_file_path)

    try:
        if first_line_of_log_is_cmdline:
            with open(log_file_path, "a") as log_file: log_file.write('# DATE: %s\n# CMD LINE: %s\n' % (get_date(), ' '.join(cmdline)))

        log_file = open(log_file_path, 'a')
        ret_val = subprocess.call(cmdline, shell=False, stdout=log_file, stderr=subprocess.STDOUT)
        log_file.close()

        # This can happen in POSIX due to signal termination (e.g., SIGKILL).
        if ret_val < 0:
            raise ConfigError("Command failed to run. What command, you say? This: '%s'" % ' '.join(cmdline))
        else:
            return ret_val
    except OSError as e:
        raise ConfigError("command was failed for the following reason: '%s' ('%s')" % (e, cmdline))



def start_command(cmdline, log_file_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, first_line_of_log_is_cmdline=True, remove_log_file_if_exists=True):
    """Start a command using subprocess.Popen, returning an object that can be monitored."""
    cmdline = format_cmdline(cmdline)

    if anvio.DEBUG:
        Progress().reset()
        Run().info("[DEBUG] `start_command`",
                   ' '.join(['%s' % (('"%s"' % str(x)) if ' ' in str(x) else ('%s' % str(x))) for x in cmdline]),
                   nl_before=1, nl_after=1, mc='red', lc='yellow')

    filesnpaths.is_output_file_writable(log_file_path)

    if remove_log_file_if_exists and os.path.exists(log_file_path):
        os.remove(log_file_path)

    try:
        if first_line_of_log_is_cmdline:
            with open(log_file_path, 'a') as log_file:
                log_file.write(f"# DATE: {get_date()}\n# CMD LINE: {' '.join(cmdline)}\n")

        p = subprocess.Popen(cmdline, stdout=stdout, stderr=stderr)
        return p
    except OSError as e:
        raise ConfigError("The command failed for the following reason: '%s' ('%s')" % (e, cmdline))



def run_command_STDIN(cmdline, log_file_path, input_data, first_line_of_log_is_cmdline=True, remove_log_file_if_exists=True):
    """Uses subprocess.Popen and sends data to your `cmdline` through STDIN"""
    cmdline = format_cmdline(cmdline)

    filesnpaths.is_output_file_writable(log_file_path)

    if remove_log_file_if_exists and os.path.exists(log_file_path):
        os.remove(log_file_path)

    try:
        if first_line_of_log_is_cmdline:
            with open(log_file_path, "a") as log_file: log_file.write('# DATE: %s\n# CMD LINE: %s\n' % (get_date(), ' '.join(cmdline)))

        p = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        ret_val = p.communicate(input=input_data.encode('utf-8'))[0]
        return ret_val.decode()
    except OSError as e:
        raise ConfigError("command was failed for the following reason: '%s' ('%s')" % (e, cmdline))



def get_command_output_from_shell(cmd_line):
    ret_code = 0

    try:
        out_bytes = subprocess.check_output(cmd_line.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        out_bytes = e.output.decode("utf-8")
        ret_code = e.returncode

    return out_bytes, ret_code



def get_cmd_line():
    c_argv = []
    for i in sys.argv:
        if ' ' in i:
            c_argv.append('"%s"' % i)
        else:
            c_argv.append(i)
    return ' '.join(c_argv)



def get_gene_caller_ids_from_args(gene_caller_ids, delimiter=','):
    gene_caller_ids_set = set([])
    if gene_caller_ids:
        if os.path.exists(gene_caller_ids):
            gene_caller_ids_set = set([g.strip() for g in open(gene_caller_ids, 'r').readlines()])
        else:
            gene_caller_ids_set = set([g.strip() for g in gene_caller_ids.split(delimiter)])

    try:
        gene_caller_ids_set = set([int(float(g)) for g in gene_caller_ids_set])
    except:
        g = gene_caller_ids_set.pop()
        raise ConfigError("The gene calls you provided do not look like gene callers anvi'o is used to working with :/ Here is "
                          "one of them: '%s' (%s)." % (g, type(g)))
    return gene_caller_ids_set


class RunInDirectory(object):
    """ Run any block of code in a specified directory. Return to original directory

    Parameters
    ==========
    run_dir : str or Path-like
        The directory the block of code should be run in
    """

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.cur_dir = os.getcwd()
        if not os.path.isdir(self.run_dir):
            raise ConfigError("RunInDirectory :: %s is not a directory." % str(self.run_dir))


    def __enter__(self):
        os.chdir(self.run_dir)


    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.cur_dir)
