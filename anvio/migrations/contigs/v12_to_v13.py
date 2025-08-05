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
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Dropping the HMMs ")
    progress.update("...")
    for table_name in ['hmm_hits_info', 'hmm_hits', 'hmm_hits_in_splits']:
        contigs_db.remove_some_rows_from_table(table_name, 'source IN ("Rinke_et_al", "Campbell_et_al", "BUSCO_83_Protista")')

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()
    run.info_single("The contigs database is now %s. Unfortunately this update removed ALL SINGLE-COPY CORE GENE "
                    "HMMs FROM YOUR CONTIGS DATABASE :( We are very sorry about this, but we only did it to be "
                    "able to offer you nicer things. It is best if you re-run `anvi-run-hmms` program from scratch. "
                    "Doing that will not remove any 'non-default' HMM profiles you may have added in this contigs "
                    "database, so you have nothing to worry about." % (next_version), nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
