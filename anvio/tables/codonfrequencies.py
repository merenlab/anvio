# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table


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
pp = terminal.pretty_print


class TableForCodonFrequencies(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, utils.get_required_version_for_db(db_path), run=self.run, progress=self.progress)

        self.num_entries = 0
        self.db_entries = []
        self.set_next_available_id(t.variable_codons_table_name)

        self.max_num_entries_in_storage_buffer = 15000


    def append_entry(self, entry):
        self.db_entries.append(entry)

        if len(self.db_entries) > self.max_num_entries_in_storage_buffer:
            self.store()


    def append(self, entry):
        """Append a single entry based on a dictionary

        Parameters
        ==========
        entry : sequence
            values in order they are in the table, entry_id excluded (it will be appended in the
            body of this function)
        """

        db_entry = (self.next_id(t.variable_codons_table_name), *entry)
        self.db_entries.append(db_entry)
        self.num_entries += 1

        if len(self.db_entries) >= self.max_num_entries_in_storage_buffer:
            # everytime we are here, the contenst of self.db_entries will be stored in the
            # database
            self.store()


    def store(self):
        if not len(self.db_entries):
            return

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.variable_codons_table_name, ','.join(['?'] * len(t.variable_codons_table_structure))), self.db_entries)
        database.disconnect()

        if anvio.DEBUG:
            run.info_single("SCVs: %d entries added to the nt variability table." % len(self.db_entries), mc="blue")

        self.db_entries = []
