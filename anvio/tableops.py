# -*- coding: utf-8

""" Table schemas for databases."""

import os

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()



####################################################################################################
#
#     TABLE SUPERCLASS
#
####################################################################################################


class Table(object):
    """Superclass for rudimentary needs and operations for contigs db tables"""
    def __init__(self, db_path, version, run=run, progress=progress, quiet = False):
        if not os.path.exists(db_path):
            raise ConfigError, "Database ('%s') does not exist. You must create one first." % db_path

        self.quiet = quiet
        self.db_type = None
        self.db_path = db_path
        self.version = version
        self.next_available_id = {}

        self.splits_info = None
        self.contigs_info = None
        self.split_length = None

        self.run = run
        self.progress = progress

        database = db.DB(self.db_path, version)
        self.db_type = database.get_meta_value('db_type')

        if self.db_type == 'contigs': 
            # FIXME: a better design is required. the salient point is, "Table" must serve for both profile db
            # and contigs db calls.
            self.split_length = database.get_meta_value('split_length')
            self.contigs_info = database.get_table_as_dict(t.contigs_info_table_name, string_the_key = True)
            self.splits_info  = database.get_table_as_dict(t.splits_info_table_name)
            self.contig_name_to_splits = utils.get_contig_name_to_splits_dict(self.splits_info, self.contigs_info)

        database.disconnect()


    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError, "If you need unique ids, you must call 'set_next_available_id' first"

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1


    def set_next_available_id(self, table):
        database = db.DB(self.db_path, self.version)
        table_content = database.get_table_as_dict(table)
        if table_content:
            self.next_available_id[table] = max(table_content.keys()) + 1
        else:
            self.next_available_id[table] = 0

        database.disconnect()


    def export_contigs_in_db_into_FASTA_file(self):
        if self.db_type != 'contigs':
            return None

        database = db.DB(self.db_path, self.version)
        contig_sequences_table = database.get_table_as_dict(t.contig_sequences_table_name)
        database.disconnect()

        self.progress.new('Exporting contigs into a FASTA file')
        self.progress.update('...')
        contigs_fasta_path = os.path.join(filesnpaths.get_temp_directory_path(), 'contigs.fa')
        contigs_fasta = u.FastaOutput(contigs_fasta_path)
        for contig in contig_sequences_table:
            contigs_fasta.write_id(contig)
            contigs_fasta.write_seq(contig_sequences_table[contig]['sequence'], split=False)

        self.progress.end()
        self.run.info('FASTA for contigs', contigs_fasta_path)

        return contigs_fasta_path


    def delete_entries_for_key(self, table_column, key, tables_to_clear = []):
        # FIXME: this shoudl be in db.py
        # removes rows from each table in 'tables_to_remove' where 'table_column' equals 'value'
        database = db.DB(self.db_path, self.version)

        table_content = database.get_table_as_dict(tables_to_clear[0])
        if key in table_content:
            self.run.warning('Previous entries for "%s" is being removed from "%s"' % (key, ', '.join(tables_to_clear)))
            for table_name in tables_to_clear:
                database._exec('''DELETE FROM %s WHERE %s = "%s"''' % (table_name, table_column, key))

        database.disconnect()


