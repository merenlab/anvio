#!/usr/bin/env python
# -*- coding: utf-8

import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_trnaseq_db


current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    is_trnaseq_db(db_path)

    trnaseq_db = db.DB(db_path, None, ignore_version=True)
    if str(trnaseq_db.get_version()) != current_version:
        raise ConfigError("Version of this tRNA-seq database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Migrating things")

    progress.update("Dealing with the self table")
    trnaseq_db.set_meta_value('INDELs_profiled', True)

    progress.update("Updating version")
    trnaseq_db.remove_meta_key_value_pair('version')
    trnaseq_db.set_version(next_version)

    progress.update("Committing changes")
    trnaseq_db.disconnect()

    progress.end()

    message = ("Done! Your tRNA-seq db is now version 2. "
               "Here is a little information on why this update was necessary (in case you have nothing else to read): "
               "1. The prediction of INDELs in a tRNA sequence is time-consuming and may not add much value to your analysis. "
               "You now have the option of skipping INDEL profiling, which always occurred for v1 databases. ")

    run.info_single(message, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade TRNASEQ.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('trnaseq_db', metavar = 'TRNASEQ_DB', help = 'tRNA-seq database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.trnaseq_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
