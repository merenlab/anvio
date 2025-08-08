# -*- coding: utf-8
# pylint: disable=line-too-long

""" Table schemas for databases."""

import os

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"



####################################################################################################
#
#     TABLE SUPERCLASS
#
####################################################################################################


class Table(object):
    """Superclass for rudimentary needs and operations for db tables

    Notes
    =====
    - This class needs to be redesigned so it better serves the purpose for both profile db and
      contigs db calls
    """

    def __init__(self, db_path, version, run=terminal.Run(), progress=terminal.Progress(), quiet=False, simple=False):
        if not db_path:
            raise ConfigError("Table superclass is being initiated without a db path, and it is very "
                               "very concerning :( Anvi'o needs an adult.")

        if not os.path.exists(db_path):
            raise ConfigError("Database ('%s') does not exist. You must create one first." % db_path)

        self.quiet = quiet
        self.db_type = None
        self.db_path = db_path
        self.version = version
        self.next_available_id = {}

        self.split_length = None
        self.genes_are_called = None

        self.run = run
        self.progress = progress

        database = db.DB(self.db_path, version)
        self.db_type = database.get_meta_value('db_type')
        if not simple and self.db_type == 'contigs':
            # FIXME: a better design is required. the salient point is, "Table" must serve for both profile db
            # and contigs db calls.
            self.split_length = database.get_meta_value('split_length')
            self.contigs_info = database.get_table_as_dict(t.contigs_info_table_name, string_the_key=True)
            self.splits_info = database.get_table_as_dict(t.splits_info_table_name)
            self.genes_are_called = database.get_meta_value('genes_are_called')
            self.gene_calls_dict = None
        database.disconnect()

        if not simple and self.db_type == 'contigs':
            self.contig_name_to_splits = utils.get_contig_name_to_splits_dict(self.db_path)


    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError("If you need unique ids, you must call 'set_next_available_id' first")

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1


    def set_next_available_id(self, table):
        # FIXME: This could be a lot faster if entry IDs are ordered. Then we just need the last
        #        entry of the 'entry_id' column.
        database = db.DB(self.db_path, self.version)
        self.next_available_id[table] = database.get_max_value_in_column(table, 'entry_id', value_if_empty=-1) + 1
        database.disconnect()


    def reset_next_available_id_for_table(self, table):
        self.next_available_id[table] = 0


    def delete_entries_for_key(self, table_column, key, tables_to_clear=[]):
        # FIXME: this should be in db.py
        # removes rows from each table in 'tables_to_remove' where 'table_column' equals 'value'
        database = db.DB(self.db_path, self.version)

        table_content = database.get_table_as_dict(tables_to_clear[0])
        if key in table_content:
            self.run.warning('Previous entries for "%s" is being removed from "%s"' % (key, ', '.join(tables_to_clear)))
            for table_name in tables_to_clear:
                database._exec('''DELETE FROM %s WHERE %s = "%s"''' % (table_name, table_column, key))

        database.disconnect()


    def delete_contents_of_table(self, table_name, warning=True):
        database = db.DB(self.db_path, self.version)

        if warning:
            self.run.warning('Contents of the table "%s" is being removed' % (table_name))

        database._exec('''DELETE FROM %s''' % (table_name))

        database.disconnect()


    def init_gene_calls_dict(self):
        if self.db_type != 'contigs':
            return None

        self.progress.new('Initializing the dictionary for gene calls')
        self.progress.update('...')

        database = db.DB(self.db_path, self.version)
        self.gene_calls_dict = database.get_table_as_dict(t.genes_in_contigs_table_name)
        database.disconnect()

        self.progress.end()


