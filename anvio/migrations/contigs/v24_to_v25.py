#!/usr/bin/env python

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

contig_classification_table_name      = 'contig_classification'
contig_classification_table_structure = ['contig', 'class', 'source', 'tool_classification', 'confidence']
contig_classification_table_types     = [ 'text' ,'numeric', 'text',        'text'          ,   'text'   ]

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(f"The version of the provided contigs database is {contigs_db.get_version()}, not the "
                          f"required version, {current_version}, so this script cannot upgrade the database.")

    progress.new("Adding contig classification table")
    progress.update("...")

    try:
        contigs_db.drop_table(contig_classification_table_name)
    except:
        pass

    contigs_db.create_table(contig_classification_table_name,
                             contig_classification_table_structure,
                             contig_classification_table_types)
    progress.end()

    progress.new("Updating version")
    progress.update("...")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)
    contigs_db.disconnect()
    progress.end()

    run.info_single(f"Congratulations! Your contigs database is now version {next_version}. An empty "
                    f"contig_classification table has been added to support contig-level classification "
                    f"via anvi-import-contig-classification.",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar='CONTIGS_DB', help='Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
