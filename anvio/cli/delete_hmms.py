#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.hmmops as hmmops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.hmmhits import TablesForHMMHits
from anvio.dbops import ContigsDatabase


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db", "hmm-source", "hmm-hits"]
__description__ = "Remove HMM hits from an anvi'o contigs database"


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()

    try:
        info_table = hmmops.SequencesForHMMHits(args.contigs_db).hmm_hits_info

        if not len(info_table):
            run.warning("The HMM tables in your contigs databse is empty. Now anvi'o will quit and go back to sleep.")
            sys.exit()

        if args.list_hmm_sources:
            ContigsDatabase(args.contigs_db).list_available_hmm_sources()
            sys.exit()

        if not args.hmm_source and not args.just_do_it:
            raise ConfigError("But you should give this program an HMM souce to remove from the database :/ If you "
                              "want everything to be removed, you can simply add the flag `--just-do-it`. In that "
                              "case anvi'o will remove every HMM source from your database without asking anything.")

        if args.hmm_source and args.just_do_it:
            raise ConfigError("You are both asking anvi'o to remove every HMM source (with `--just-do-it`), and "
                              "specifying a paticular HMM source to be removed (with `--hmm-source`). Since anvi'o "
                              "has this suspicion that you may have no idea what you are doing, she chooses to not "
                              "do anything at all (anvi'o the pacifist). Choose either, and let's pretend this never "
                              "happened.")

        hmm_tables = TablesForHMMHits(args.contigs_db, initializing_for_deletion=True)
        if args.hmm_source:
            if args.hmm_source not in info_table:
                run.warning("Nice try, and anvi'o is proud of you. But it doesn't change the fact that the HMM source "
                            "'%s' was not in your contigs database. That's OK. Now anvi'o will pretend as if nothing "
                            "happened, and quit without a fuss." % args.hmm_source)
                sys.exit(0)

            hmm_tables.remove_source(args.hmm_source)
        elif args.just_do_it:
            for source in info_table:
                hmm_tables.remove_source(source)
        else:
            raise ConfigError("Anvi'o does not know how the hell you ended up in this part of its code. Please find a "
                              "developer, explain how you ended up here, and collect your prize.")
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('hmm-source'), **anvio.K('hmm-source'))
    parser.add_argument(*anvio.A('list-hmm-sources'), **anvio.K('list-hmm-sources'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
