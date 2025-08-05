#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_contigs_db

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(f"The version of the provided contigs database is {contigs_db.get_version()}, "
                          f"not the required version, {current_version}, so this script cannot upgrade the database.")

    db_altered = False
    progress.new("Migrating")

    # usually we never need to do this, but in this case it is best practice to test if
    # we need an actual removal of the SCGs from a contigs-db since we released new SCGs
    # in the master repository, and some people likely already updated their contigs-dbs.
    # if they have already gone through that, we can save them the trouble.
    scg_taxonomy_was_run = contigs_db.get_meta_value('scg_taxonomy_was_run')
    scg_taxonomy_db_version = contigs_db.get_meta_value('scg_taxonomy_database_version')
    if scg_taxonomy_was_run and scg_taxonomy_db_version == "GTDB: v214.1; Anvi'o: v1":
        # does not need an update.
        pass
    else:
        # needs an update
        contigs_db._exec('''DELETE FROM scg_taxonomy''')
        contigs_db.set_meta_value('scg_taxonomy_was_run', 0)
        contigs_db.set_meta_value('scg_taxonomy_database_version', None)

        db_altered = True

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()

    if db_altered:
        message = (f"Your contigs database is now version {next_version}. Sadly this update removed all SCG taxonomy "
                   f"data in this contigs-db due to a change in the set of SCGs anvi'o now uses for taxonomy estimation. "
                   f"As a result, you will need to re-run anvi-run-scg-taxonomy command on this contigs-db :/ If you "
                   f"would like to learn why this was necessary, please visit https://github.com/merenlab/anvio/issues/2211. "
                   f"We thank you for your patience!")
    else:
        message = ("Since you have already updated your contigs-db with new SCGs, anvi'o simply bumped the version of your "
                   "database rather than removing or editing any data :) Moving on.")

    run.info_single(message, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"A simple script to upgrade an anvi'o contigs database from version {current_version} to version {next_version}")
    parser.add_argument("contigs_db", metavar="CONTIGS_DB", help=f"An anvi'o contigs database of version {current_version}")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
