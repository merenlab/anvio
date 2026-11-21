#!/usr/bin/env python
DESCRIPTION = """This program downloads and sets up the ModelSEED Biochemistry database."""

from sys import exit
from argparse import Namespace

import anvio.reactionnetwork as reactionnetwork

from anvio import A, K
from anvio.errors import ConfigError
from anvio import __version__ as VERSION
from anvio.argparse import ArgumentParser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__authors__ = ['semiller10']
__requires__ = ['functions']
__provides__ = ["reaction-ref-data"]
__description__ = DESCRIPTION


def main() -> None:
    args = get_args()

    try:
        reactionnetwork.ModelSEEDDatabase.set_up(dir=args.dir, reset=args.reset)
    except ConfigError as e:
        print(e)
        exit(-1)

def get_args() -> Namespace:
    parser = ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--dir', default=None, type=str,
        help="Directory in which a new subdirectory, 'MODELSEED', is created to hold the database "
             "files. Without this flag, the database is set up in a default directory within the "
             "anvi'o installation (anvio/data/misc/MODELSEED)."
    )
    parser.add_argument(*A('reset'), **K('reset'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
