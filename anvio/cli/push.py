#!/usr/bin/env python
# -*- coding: utf-8
"""A program to 'push' analyses to anvi'server.

   Takes multiple files in, sends them to the server."""

import sys
import getpass
from anvio.argparse import ArgumentParser

import anvio

from anvio.serverAPI import AnviServerAPI
from anvio.errors import AnviServerError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ozcan']
__description__ = "Push stuff to an anvi'server"


def main():
    args = get_args()

    try:
        server = AnviServerAPI(args)
        server.login()
        server.push()
    except AnviServerError as e:
        print(e)
        sys.exit(1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SERVER DETAILS', "Details of how to access to an anvi'server instance.")
    groupA.add_argument(*anvio.A('user'), **anvio.K('user', {'required': True}))
    groupA.add_argument(*anvio.A('api-url'), **anvio.K('api-url'))


    groupB = parser.add_argument_group("PROJECT DETAILS", "What to send to the server")
    groupB.add_argument(*anvio.A('project-name'), **anvio.K('project-name', {'required': True}))
    groupB.add_argument(*anvio.A('tree'), **anvio.K('tree'))
    groupB.add_argument(*anvio.A('items-order'), **anvio.K('items-order'))
    groupB.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file'))
    groupB.add_argument(*anvio.A('view-data'), **anvio.K('view-data'))
    groupB.add_argument(*anvio.A('additional-layers'), **anvio.K('additional-layers'))
    groupB.add_argument(*anvio.A('state'), **anvio.K('state'))
    groupB.add_argument(*anvio.A('description'), **anvio.K('description'))
    groupB.add_argument(*anvio.A('bins'), **anvio.K('bins'))
    groupB.add_argument(*anvio.A('bins-info'), **anvio.K('bins-info'))
    # FIX ME: bring these back
    # groupB.add_argument(*anvio.A('layers-information-file'), **anvio.K('layers-information-file'))
    # groupB.add_argument(*anvio.A('layers-order-file'), **anvio.K('layers-order-file'))

    groupC = parser.add_argument_group("RISKY CLICKS", "As the name suggests!")
    groupC.add_argument(*anvio.A('delete-if-exists'), **anvio.K('delete-if-exists'))

    args = parser.get_args(parser)

    args.password = getpass.getpass('Enter your password for the server: ')

    return args


if __name__ == '__main__':
    main()
