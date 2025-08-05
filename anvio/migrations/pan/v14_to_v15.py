#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

def migrate(db_path):
    pan_db = db.DB(db_path, None, ignore_version = True, skip_rowid_prepend=True)

    progress.new("Bleep bloop")
    progress.update('...')

    for table_name in ['gene_cluster_frequencies', 'gene_cluster_presence_absence']:
        progress.update(f"Working on '{table_name}' ...")

        table_data = pan_db.get_table_as_dict(table_name)

        # know your genomes
        genome_names = sorted(list(list(table_data.values())[0].keys()))

        # drop the old view table
        pan_db._exec(f'DROP TABLE {table_name}')

        # create a new view table!
        pan_db._exec(f'''CREATE TABLE {table_name} (item text, layer text, value numeric)''')

        # fill in the new view data from the old format
        view_data = []
        for gc in table_data:
            for genome_name in genome_names:
                view_data.append((gc, genome_name, table_data[gc][genome_name]), )

        # populate new view table
        pan_db._exec_many('''INSERT INTO %s VALUES (?,?,?)''' % (table_name), view_data)

    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)
    pan_db.disconnect()

    progress.end()

    run.info_single(f"""The pan database is now {next_version}. This upgrade fixed an annoying design decision """
                    f"""that was made in 2014. Nobody: "", Anvi'o: "IT IS NEVER TOO LATE".""",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade pan database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('pan_db', metavar = 'PAN_DB', help = "An anvi'o pan database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
