#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.workflows as w
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ShaiberAlon', 'semiller10', 'mschecht']
__resources__ = [("Tutorial", "http://merenlab.org/2018/07/09/anvio-snakemake-workflows/")]
__tags__ = ["metagenomics", "phylogenomics", "contigs", "pangenomics"]
__description__ = ("Execute, manage, parallelize, and troubleshoot entire 'omics workflows and "
                   "chain together anvi'o and third party programs")
__requires__ = ["workflow-config"]
__provides__ = ["workflow"]



def main():
    args = get_args()
    run = terminal.Run()

    run.warning('If you publish results from this workflow, please do not forget to cite Snakemake '
                '(doi:10.1093/bioinformatics/bts480)', lc = 'yellow')

    try:
        workflows_dict = w.get_workflow_module_dict()

        if args.list_workflows:
            run.info("Available workflows", ", ".join(list(workflows_dict.keys())))
            sys.exit(0)

        # FIXME: Meren and Alon should discuss these next lines
        # we can't call the snake_file_path from the class so I think
        # we have to do it this way
        if (not args.workflow) or (not args.config_file) and (not args.get_default_config):
            raise ConfigError("You must provide a workflow name AND a config file. You can use --list-workflow "
                              "to learn what workflows are available, and you can use --get-default-config "
                              "if you need help writing your config file.")

        if args.workflow not in workflows_dict:
            raise ConfigError("anvi'o is not familiar with a workflow %s, you can use --list-workflow to "
                              "learn what workflows are available." % args.workflow)

        M = workflows_dict[args.workflow](args)
        M.init()
        M.go(skip_dry_run=args.skip_dry_run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('ESSENTIAL INPUTS', "Things you must provide or this won't work")
    groupA.add_argument(*anvio.A("workflow"), **anvio.K("workflow"))

    groupB = parser.add_argument_group('ADDITIONAL STUFF', "additional stuff")
    groupB.add_argument(*anvio.A("get-default-config"), **anvio.K("get-default-config"))
    groupB.add_argument(*anvio.A("list-workflows"), **anvio.K("list-workflows"))
    groupB.add_argument(*anvio.A("list-dependencies"), **anvio.K("list-dependencies"))
    groupB.add_argument(*anvio.A("config-file"), **anvio.K("config-file"))
    groupB.add_argument(*anvio.A("dry-run"), **anvio.K("dry-run"))
    groupB.add_argument(*anvio.A("skip-dry-run"), **anvio.K("skip-dry-run"))
    groupB.add_argument(*anvio.A("save-workflow-graph"), **anvio.K("save-workflow-graph"))
    groupB.add_argument(*anvio.A("additional-params"), **anvio.K("additional-params"))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
