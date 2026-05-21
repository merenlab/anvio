"""Path helper functions for anvi'o workflows."""

import os

from anvio.errors import ConfigError


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


def get_workflow_rule_file_path(workflow, filename='main.smk'):
    workflow_dir = os.path.join(get_path_to_workflows_dir(), workflow)

    if not os.path.isdir(workflow_dir):
        raise ConfigError("Anvi'o does not know about the workflow '%s' :/" % workflow)

    rule_file_path = os.path.join(workflow_dir, 'rules', filename)

    if not os.path.exists(rule_file_path):
        raise ConfigError("The rules file '%s' for the workflow '%s' seems to be missing :/" % (filename, workflow))

    return rule_file_path
