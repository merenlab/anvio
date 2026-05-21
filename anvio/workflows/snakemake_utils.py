"""Small utility functions used by anvi'o Snakefiles."""

import os
import re

import anvio.terminal as terminal


def D(debug_message, debug_log_file_path=".SNAKEMAKEDEBUG"):
    """Append a timestamped debug message to a Snakemake debug log."""
    with open(debug_log_file_path, 'a') as output:
            output.write(terminal.get_date() + '\n')
            output.write(str(debug_message) + '\n\n')


def regex_from_ids(ids):
    """Return a wildcard constraint regex from a list of IDs."""
    return r"(?:{})".format("|".join(re.escape(str(i)) for i in ids)) if ids else r"DO_NOT_MATCH"


def get_conda_yaml_path(workflow, tool):
    """Return the absolute conda YAML path configured for a workflow tool."""
    path = workflow.get_param_value_from_config([tool, 'conda_yaml'])
    if not path:
        return None

    return os.path.abspath(os.path.expanduser(path))


def get_conda_env_prefix(workflow, tool):
    """Return the command prefix for a configured existing conda environment."""
    name = workflow.get_param_value_from_config([tool, 'conda_env'])
    return f"conda run -n {name}" if name else ""


def gunzip_file(file_path, log, shell, output_path=None):
    """Uncompress a gzipped file and return the uncompressed path."""
    uncompressed_file_path = output_path or os.path.splitext(file_path)[0]
    shell("gunzip < %s > %s 2>> %s" % (file_path, uncompressed_file_path, log))
    return uncompressed_file_path
