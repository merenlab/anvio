#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.states import TablesForStates


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db', 'profile-db', 'state-json']
__provides__ = ['state']
__description__ = "Import an anvi'o state into a profile database"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_pan_or_profile_db(args.pan_or_profile_db, genes_db_is_also_accepted=True)
        filesnpaths.is_file_json_formatted(args.state)
        utils.is_this_name_OK_for_database('--name parameter', args.name)

        TablesForStates(args.pan_or_profile_db).store_state(args.name, open(args.state).read())

        run.info('Done', 'State "%s" is added to the database' % args.name)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))
    parser.add_argument('-s', '--state', metavar="STATE_FILE", default = None, required = True,
                        help = "JSON serializable anvi'o state file.")
    parser.add_argument('-n', '--name', metavar = 'STATE_NAME', default = None, required = True,
                        help = 'State name.')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
