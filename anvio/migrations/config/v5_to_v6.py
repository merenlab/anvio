#!/usr/bin/env python

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

# The optional long-read QC rules added to the metagenomics workflow in v6, with their
# default parameter blocks. These are frozen snapshots of the defaults as they were when
# this migration was written (do not update them if the schema changes later — that is what
# future migrations are for). They are added only to metagenomics configs, and only when a
# rule is not already present, so any values the user already set are preserved untouched.
NEW_METAGENOMICS_QC_RULES = {
    'fastqc_sr': {"threads": 1, "run": False, "run_on_raw": False, "run_on_filtered": False,
                  "conda_yaml": "", "use_anvio_conda_yaml": True, "conda_env": "",
                  "additional_params": ""},
    'filtlong':  {"threads": 1, "run": False, "conda_yaml": "", "use_anvio_conda_yaml": True,
                  "conda_env": "", "--min-length": None, "--max-length": None,
                  "--target-bases": None, "additional_params": ""},
    'nanoplot':  {"threads": 2, "run": False, "run_on_raw": False, "run_on_filtered": False,
                  "conda_yaml": "", "use_anvio_conda_yaml": True, "conda_env": "",
                  "additional_params": ""},
    'multiqc':   {"threads": 1, "run": False, "conda_yaml": "", "use_anvio_conda_yaml": True,
                  "conda_env": "", "additional_params": ""},
}

# v6 also introduced the `use_anvio_conda_yaml` option (default true → use the conda env anvi'o
# ships) on these pre-existing metagenomics tools. New v6 configs get it as `true`, but on a
# MIGRATED config we set it to `false` on purpose (see the loop below): a working v5 config ran
# these tools from $PATH, and flipping them to `true` would make Snakemake build conda envs
# instead — breaking offline/no-conda runs. This frozen list covers the tools whose blocks
# already existed in v5. Like the blocks above, it is a snapshot of v6 and must not be updated
# for later schema changes.
PRE_EXISTING_METAGENOMICS_CONDA_TOOLS = ['megahit', 'metaspades', 'idba_ud',
                                         'flye', 'bowtie', 'minimap2']


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

    # v6 adds optional QC rules (Filtlong, NanoPlot, FastQC, MultiQC) to the metagenomics
    # workflow. Add their default blocks to metagenomics configs if they are not already there.
    # Every rule defaults to `run: False`, so this changes no existing behavior.
    added_rules = []
    path_pinned_tools = []
    if workflow_name == 'metagenomics':
        for rule, default_block in NEW_METAGENOMICS_QC_RULES.items():
            if rule not in config:
                config[rule] = dict(default_block)
                added_rules.append(rule)

        # v6 added the `use_anvio_conda_yaml` option to the pre-existing assembly/mapping tools.
        # On a migrated config we set it to `false` (NOT the v6 default of true): a v5 config that
        # lacked `conda_yaml`/`conda_env` ran the tool from $PATH, and setting `true` would make
        # Snakemake build a conda env for it instead (go()/dry_run() add --use-conda, and the $PATH
        # check is skipped) — turning a working offline/no-conda run into a failing one. `false`
        # keeps the exact v5 behavior; users can opt into the anvi'o-shipped env by editing it to
        # true. We only add the key when absent, and skip tools that already declare their own
        # `conda_yaml`/`conda_env` (those three are mutually exclusive, see MetagenomicsWorkflow.init).
        for tool in PRE_EXISTING_METAGENOMICS_CONDA_TOOLS:
            block = config.get(tool)
            if not isinstance(block, dict) or 'use_anvio_conda_yaml' in block:
                continue
            if block.get('conda_yaml') or block.get('conda_env'):
                continue
            block['use_anvio_conda_yaml'] = False
            path_pinned_tools.append(tool)

    # set it to the new version
    config['config_version'] = next_version

    with open(config_path, 'w') as output:
        output.write(json.dumps(config, indent=4))

    progress.end()

    if added_rules or path_pinned_tools:
        parts = [f"The config file version is now {next_version}."]
        if added_rules:
            parts.append(f"This upgrade added the following optional long-read QC rules to your "
                         f"metagenomics config (all disabled by default, so nothing changes unless you "
                         f"turn them on): {', '.join(added_rules)}.")
        if path_pinned_tools:
            parts.append(f"It also set 'use_anvio_conda_yaml: false' on these tools that did not "
                         f"already declare their own 'conda_yaml'/'conda_env', so they keep running "
                         f"from your $PATH exactly as before: {', '.join(path_pinned_tools)}. If you "
                         f"would rather have anvi'o build a dedicated conda environment for any of "
                         f"them, set 'use_anvio_conda_yaml: true' for that tool.")
        parts.append("One related note for long-read users: anvi'o no longer applies a silent default "
                     "minimap2 preset or Flye read-type — provide them either via the new `lr_technology` "
                     "column in your samples-txt (recommended) or explicitly in this config; your existing "
                     "values (if any) were left untouched. See the metagenomics workflow documentation for details.")
        message = ' '.join(parts)
    else:
        message = (f"The config file version is now {next_version}. No rule changes were necessary for the "
                   f"'{workflow_name}' workflow — this was just a version stamp update.")

    run.info_single(message, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"The migration script to migrate an anvi'o workflow config from "
            f"v{current_version} to v{next_version}")
    parser.add_argument('config', metavar = 'CONFIG', help = f'A v{current_version} config JSON.')
    parser.add_argument(*anvio.A("workflow"), **anvio.K("workflow"))
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.config)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
