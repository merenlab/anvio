#!/usr/bin/env python
"""Entry point to the local anvi'o workflow builder interface."""

import sys

import anvio
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError
from anvio.workflowbuilder import WorkflowBuilder

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['bastien-camillo']
__requires__ = []
__provides__ = []
__can_provide__ = []
__description__ = "Launch a local web interface to design anvi'o workflows"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)
        WorkflowBuilder(args, run=run).serve()
    except ConfigError as error:
        print(error)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SERVER CONFIGURATION', "Options for the local workflow builder server.")
    groupA.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupA.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupA.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupA.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
