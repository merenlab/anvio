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
        raise ConfigError("No config path is given :/")

    # do we have write access?
    filesnpaths.is_output_file_writable(config_path)

    # learn the workflow name and the version
    workflow_name, version = w.get_workflow_name_and_version_from_config(config_path, dont_raise=True)

    # is this the right version?
    if version != current_version:
        raise ConfigError(f"This config file is not v{current_version} (hence, this script cannot really do anything).")

    progress.new("Migrating the config file")
    progress.update("...")

    # You must skip version check, otherwise anvi'o complains that the config is out of date. DUH,
    # that's why we are migrating.
    anvio.QUIET = True
    args = argparse.Namespace(workflow=workflow_name, config_file=config_path, skip_version_check=True)
    workflow_object = w.get_workflow_module_dict()[workflow_name](args)
    config = workflow_object.config
    anvio.QUIET = False

    if 'gunzip_fasta' in config:
        del config['gunzip_fasta']

    # set it to the new version
    config['config_version'] = next_version

    with open(config_path, 'w') as output:
        output.write(json.dumps(config, indent=4))

    progress.end()
    run.info_single(f"The config file version is now {next_version}. This upgrade removes the `gunzip_fasta` "
                    f"rule from your config (if it was there). The contigs workflow no longer needs to decompress "
                    f"gzipped FASTA files before processing them, since all downstream tools (anvi-script-reformat-fasta, "
                    f"anvi-gen-contigs-database, bowtie2-build, minimap2) natively support gzipped input. In addition, "
                    f"the contigs workflow had additional changes that will make your *existing* `anvi-run-workflow` outputs "
                    f"not recognized with the new version of the workflow :/ BUT we have made available a migration script "
                    f"at https://github.com/merenlab/anvio/pull/2553 to help you with this backwards compatibility issue. "
                    f"Please check out the link for more details, and reach out to us on anvi'o Discord if you have any "
                    f"questions regarding how to complete the migration process for your data.",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"The migration script to migrate an anvi'o workflow config from "
            f"v{current_version} to v{next_version}")
    parser.add_argument('config', metavar = 'CONFIG', help = 'A v{current_version} config JSON.')
    parser.add_argument(*anvio.A("workflow"), **anvio.K("workflow"))
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.config)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
