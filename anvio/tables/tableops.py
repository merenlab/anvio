# -*- coding: utf-8
# pylint: disable=line-too-long

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


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
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
    def __init__(self, db_path, version, run=run, progress=progress, quiet=False, simple=False):
        if not db_path:
            raise ConfigError("Table superclass is being initiated without a db path, and it is very\
                                very concerning :( Anvi'o needs an adult.")

        if not os.path.exists(db_path):
            raise ConfigError("Database ('%s') does not exist. You must create one first." % db_path)

        self.quiet = quiet
        self.db_type = None
        self.db_path = db_path
        self.version = version
        self.next_available_id = {}

        self.splits_info = None
        self.contigs_info = None
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
            self.genes_are_called = database.get_meta_value('genes_are_called')
            self.contigs_info = database.get_table_as_dict(t.contigs_info_table_name, string_the_key=True)
            self.splits_info = database.get_table_as_dict(t.splits_info_table_name)
            self.contig_name_to_splits = utils.get_contig_name_to_splits_dict(self.splits_info, self.contigs_info)
            self.gene_calls_dict = None

        database.disconnect()


    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError("If you need unique ids, you must call 'set_next_available_id' first")

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1


    def set_next_available_id(self, table):
        # FIXME: This is one of the least efficient funcitons ever.
        #        we need to get the max int value found in the first
        #        column of the table .. there is no need to get the
        #        entire table as a dict! we even have a function for
        #        that
        #
        #            db.get_max_value_in_column(table_name, 'entry_id')
        #
        database = db.DB(self.db_path, self.version)
        table_content = database.get_table_as_dict(table)
        if table_content:
            self.next_available_id[table] = max(table_content.keys()) + 1
        else:
            self.next_available_id[table] = 0

        database.disconnect()


    def reset_next_available_id_for_table(self, table):
        self.next_available_id[table] = 0


    def export_sequences_table_in_db_into_FASTA_file(self, table=t.contig_sequences_table_name, output_file_path=None, item_names=set([])):
        '''Exports a sequence table from the contigs database.

            - t.contig_sequences_table_name: contig sequences (where item_names are contig names)
            - t.gene_amino_acid_sequences_table_name: amino acid sequences for gene calls (item_names are gene caller ids)


          If `item_names` are specified, only those sequences with matching ids to something in this set will be reported.
          '''

        if self.db_type != 'contigs':
            return None

        if not isinstance(item_names, set):
            raise ConfigError("`item_names` must be of type `set`")

        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)
        else:
            output_file_path = os.path.join(filesnpaths.get_temp_directory_path(), 'aa_sequences.fa')

        database = db.DB(self.db_path, self.version)

        if table not in database.get_table_names():
            raise ConfigError('Trying to export sequences into a FASTA file, but the table\
                                "%s" does not seem to be in this database :/' % (table))

        if 'sequence' not in database.get_table_structure(table):
            raise ConfigError("You requested to store sequences in table '%s' into a FASTA\
                                file, however this table does not seem to be a table that\
                                stores sequence information :(" % table)

        sequences_table = database.get_table_as_dict(table)
        database.disconnect()

        if len(item_names):
            total_num_items_in_db = len(sequences_table)
            item_names_to_remove = set(list(sequences_table.keys())).difference(item_names)

            for item_name in item_names_to_remove:
                if item_name in item_names_to_remove:
                    sequences_table.pop(item_name)

            # who does this for their users:
            num_items_to_be_reported = len(sequences_table)
            optional_info = ("It turned out %d of the item ids you requested was actually in the database." \
                                    % len(sequences_table)) if num_items_to_be_reported != len(item_names) else ''

            if not self.quiet:
                self.run.info_single("You asked anvi'o to report only %d items from a database that contained %d. %s" \
                                            % (len(item_names), total_num_items_in_db, optional_info))

        if not len([sequences_table]):
            raise ConfigError("There are no sequences to report in table '%s'." % (table))

        self.progress.new('Exporting %d sequences into a FASTA file' % len(sequences_table))
        self.progress.update('...')

        sequences_fasta = u.FastaOutput(output_file_path)

        blank_seq_ids_not_reported = set([])

        for seq_id in sequences_table:
            if len(sequences_table[seq_id]['sequence']):
                sequences_fasta.write_id(seq_id)
                sequences_fasta.write_seq(sequences_table[seq_id]['sequence'], split=False)
            else:
                blank_seq_ids_not_reported.add(seq_id)

        self.progress.end()

        if len(blank_seq_ids_not_reported):
            self.run.warning("%d entries in the sequences table had blank sequences :/ This is related to the issue\
                             at https://github.com/merenlab/anvio/issues/565. If this is like mid-2020 and you still\
                             are getting this warning, please find an anvi'o developer and make them feel embarrassed.\
                             If it is earlier than that, then take this as a simple warning to remember that some gene\
                             calls in your downstream analyses may have no amino acid sequences, and that's actuall OK.\
                             This is a very minor issue due to on-the-fly addition of Ribosomal RNA gene calls to the\
                             contigs database, and will unlikely affect anything major. This warning will go away when\
                             anvi'o can seamlessly work with multiple gene callers (which we are looking forward to\
                             implement in the future)." % len(blank_seq_ids_not_reported))

        self.run.info('Sequences', '%d sequences reported.' % (len(sequences_table) - len(blank_seq_ids_not_reported)))
        self.run.info('FASTA', output_file_path)

        return output_file_path


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

