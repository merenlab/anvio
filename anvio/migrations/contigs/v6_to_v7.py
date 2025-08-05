#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_contigs_db


run = terminal.Run()
progress = terminal.Progress()

current_version = '6'
next_version    = '7'

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_contigs_db(db_path)

    # make sure the version is 2
    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Trying to upgrade the contigs database")
    progress.update('...')

    # drop the old tables:
    try:
        contigs_db._exec('''DROP TABLE %s''' % (t.hmm_hits_info_table_name))
        contigs_db._exec('''DROP TABLE %s''' % (t.hmm_hits_table_name))
        contigs_db._exec('''DROP TABLE %s''' % (t.hmm_hits_splits_table_name))
    except:
        pass
    contigs_db.commit()

    # create new empty ones
    contigs_db.create_table(t.hmm_hits_info_table_name, t.hmm_hits_info_table_structure, t.hmm_hits_info_table_types)
    contigs_db.create_table(t.hmm_hits_table_name, t.hmm_hits_table_structure, t.hmm_hits_table_types)
    contigs_db.create_table(t.hmm_hits_splits_table_name, t.hmm_hits_splits_table_structure, t.hmm_hits_splits_table_types)

    # set the version
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    # bye
    contigs_db.disconnect()

    # bye
    progress.end()
    run.info_single("The contigs database is now %s! It no longer has any HMM hits :/ Don't "
                    "forget to run `anvi-run-hmms` on it!" % (next_version), nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade contigs database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database')
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
