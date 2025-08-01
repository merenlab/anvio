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
__requires__ = ['pan-db', 'profile-db', 'state']
__provides__ = ['state-json',]
__description__ = "Export an anvi'o state into a profile database"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_pan_or_profile_db(args.pan_or_profile_db, genes_db_is_also_accepted=True)

        states_access = TablesForStates(args.pan_or_profile_db)
        states = states_access.states

        if not len(states):
            raise ConfigError("But there are no states in this %s database :/" % utils.get_db_type(args.pan_or_profile_db))

        if args.list_states:
            states_access.list_states()
            sys.exit()

        if not args.output_file:
            raise ConfigError("You should provide an output file name (and in fact a state name, too, in case "
                               "you haven't. You can see all available state names using '--list-states' flag).")

        filesnpaths.is_output_file_writable(args.output_file)

        if not args.state:
            raise ConfigError("Sorry! you need to provide state name. Yes. You don't know what states ara available? "
                               "did you try '--list-states' flag? You can press `1` for yes, `2` for no, and `3` for "
                               "maybe.")

        if not args.state in states:
            raise ConfigError("The state name '%s' does not seem to appear in this database. But we have %s instead: "
                               "%s." % ((args.state, 'this one' if len(states) == 1 else 'these ones', ', '.join(list(states.keys())))))

        open(args.output_file, 'w').write(states[args.state]['content'])
        run.info('Output', args.output_file)
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
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    parser.add_argument('-s', '--state', metavar="STATE_NAME", default = None,
                        help = "The state name to export.")
    parser.add_argument(*anvio.A('list-states'), **anvio.K('list-states'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
