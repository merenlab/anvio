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

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

def migrate(config_path):
    if config_path is None:
        raise ConfigError("No config path is given.")

    anvio.QUIET = True
    workflow_name, version = w.get_workflow_name_and_version_from_config(config_path, dont_raise=True)

    if not workflow_name:
        raise ConfigError('Your config must include a workflow_name. For example '
                          'if this config file is used for the metagenomics workflow '
                          'then add \'"workflow_name": "metagenomics"\' to your config.')

    if version != current_version:
        raise ConfigError("Version of this config file is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Upgrading your config")
    progress.update("...")

    # You must skip version check, otherwise anvi'o complains that the config is out of date. DUH,
    # that's why we are migrating.
    args = argparse.Namespace(workflow=workflow_name, config_file=config_path, skip_version_check=True)

    workflow_object = w.get_workflow_module_dict()[workflow_name](args)
    config = workflow_object.config
    config['config_version'] = '2'

    if 'anvi_gen_contigs_database' in config:
        if '--external-gene-calls' in config['anvi_gen_contigs_database']:
            del config['anvi_gen_contigs_database']['--external-gene-calls']

    filesnpaths.is_output_file_writable(config_path)
    open(config_path, 'w').write(json.dumps(config, indent=4))

    progress.end()

    anvio.QUIET = False

    run.info_single("The config file version is now %s. This upgrade removed --external-gene-calls from "
                    "the rule anvi_gen_contigs_database, if it existed at all. This was a redundant parameter, "
                    "since external gene calls are supplied to anvi-run-workflow as an external_gene_calls "
                    "columns provided in the fasta_txt" % (next_version), nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade a workflow config from version %s to version %s' % (current_version, next_version))
    parser.add_argument('config', metavar = 'CONFIG', help = 'Config at version %s' % current_version)
    parser.add_argument(*anvio.A("workflow"), **anvio.K("workflow"))
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.config)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
