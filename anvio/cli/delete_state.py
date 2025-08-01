#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.states import TablesForStates


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["pan-db", "profile-db", "state"]
__description__ = "Delete an anvi'o state from a pan or profile database"


run = terminal.Run()


@terminal.time_program
def main():
    args = get_args()

    try:
        utils.is_pan_or_profile_db(args.pan_or_profile_db, genes_db_is_also_accepted=True)

        states_access = TablesForStates(args.pan_or_profile_db)
        states = states_access.states

        if not len(states):
            raise ConfigError("But there are no states in this %s database :/" % utils.get_db_type(args.pan_or_profile_db))

        if args.list_states:
            states_access.list_states()
            sys.exit()

        if not args.state:
            raise ConfigError("Well, you need to provide the state name you want to delete for this program to work. "
                              "You don't know which state name you want to get rid of? We heard that the '--list-states' "
                              "flag helps for that.")

        if not args.state in states:
            raise ConfigError("The state name '%s' does not seem to appear in this database. But we have %s instead: "
                               "%s." % ((args.state, 'this one' if len(states) == 1 else 'these ones', ', '.join(list(states.keys())))))

        states_access.remove_state(args.state)

        run.info_single('The state "%s" is no more' % args.state, nl_after=1)
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
    parser.add_argument('-s', '--state', metavar="STATE_NAME", default = None,
                        help = "The state name to ... delete :(")
    parser.add_argument(*anvio.A('list-states'), **anvio.K('list-states'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
