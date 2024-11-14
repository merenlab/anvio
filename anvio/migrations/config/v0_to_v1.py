#!/usr/bin/env python
# -*- coding: utf-8

import sys
import json
import argparse

import anvio
import anvio.workflows as w
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

run = terminal.Run()
progress = terminal.Progress()


def migrate(config_path):
    if config_path is None:
        raise ConfigError("No config path is given.")

    workflow_name, version = w.get_workflow_name_and_version_from_config(
        config_path, dont_raise=True
    )

    if not workflow_name:
        raise ConfigError(
            "Your config must include a workflow_name. For example "
            "if this config file is used for the metagenomics workflow "
            'then add \'"workflow_name": "metagenomics"\' to your config.'
        )

    if version != current_version:
        raise ConfigError(
            "Version of this config file is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    progress.new("Upgrading your config")
    progress.update("...")

    # You must skip version check, otherwise anvi'o complains that the config is out of date. DUH,
    # that's why we are migrating.
    args = argparse.Namespace(
        workflow=workflow_name, config_file=config_path, skip_version_check=True
    )

    workflow_module_dict = w.get_workflow_module_dict()
    workflow_object = workflow_module_dict[workflow_name](args)
    config = workflow_object.config

    default_config = workflow_object.default_config
    new_config = config.copy()
    new_config["config_version"] = "1"

    ## Deal with special cases
    special_params = ["fasta_txt", "references_for_removal", "references_mode"]

    for param in default_config:
        if param in special_params:
            # if the param belongs to special params then we skip it
            continue
        elif type(default_config[param]) == dict:
            # otherwise update config rule parameters
            new_config[param] = default_config[param]
            new_config[param].update(config.get(param, ""))
        else:
            # if it's not a dict then it's a general parameter
            # update the general parameter
            new_config[param] = config.get(param, default_config[param])

    filesnpaths.is_output_file_writable(config_path)
    open(config_path, "w").write(json.dumps(new_config, indent=4))

    progress.end()
    run.info_single(
        "The config file version is now %s. This upgrade brought back any default value that was "
        "previously removed from your config file. It will not change anything about the "
        "configuration of the resulting workflow and you can just carry on your work. "
        "" % (next_version),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade a workflow config from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "config", metavar="CONFIG", help="Config at version %s" % current_version
    )
    parser.add_argument(*anvio.A("workflow"), **anvio.K("workflow"))
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.config)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
