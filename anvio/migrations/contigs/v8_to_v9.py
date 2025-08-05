#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.tables as t
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_contigs_db


run = terminal.Run()
progress = terminal.Progress()

current_version = '8'
next_version    = '9'

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_contigs_db(db_path)

    # make sure the current version matches
    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Trying to upgrade the contigs database")
    progress.update('...')

    contigs_db.remove_meta_key_value_pair('project_name')
    contigs_db.set_meta_value('project_name', "NO_NAME")
    contigs_db.commit()

    # set the version
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    # gene name changes:
    contigs_db._exec('''UPDATE %s SET genes = replace(genes, '%s', '%s') WHERE source LIKE 'Rinke_et_al';''' % (t.hmm_hits_info_table_name, 'Ribosomal_S12', 'Ribosom_S12_S23'))
    contigs_db._exec('''UPDATE %s SET genes = replace(genes, '%s', '%s') Where source LIKE 'Rinke_et_al';''' % (t.hmm_hits_info_table_name, 'UPF0027', 'RtcB'))
    contigs_db._exec('''UPDATE %s SET genes = replace(genes, '%s', '%s') Where source LIKE 'Rinke_et_al';''' % (t.hmm_hits_info_table_name, '‐', '-')) # first - is not ASCII

    # bye
    contigs_db.disconnect()

    # bye to you too
    progress.end()

    dbops.update_description_in_db(db_path, 'No description is given')

    run.info_single("The contigs database is now %s! All this upgrade did was to associate your contigs db with a "
                    "project name (which happened to be 'NO_NAME', because anvi'o likes you very much)" \
                                            % (next_version), nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade contigs database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database')
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
