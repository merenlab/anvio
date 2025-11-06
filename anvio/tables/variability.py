# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TableForVariability(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, utils.get_required_version_for_db(db_path), run=self.run, progress=self.progress)

        self.num_entries = self.get_num_entries()
        self.db_entries = []

        # after getting an instance, we don't want things to keep accumulating
        # in memory. the purpose of the following variable is to ensure whenever
        # the number of entries in `self.db_entries` variable exceeds a certain
        # value, it will be written to the database and the global variable
        # `self.db_entries` will be emptied, saving significant memory space:
        self.max_num_entries_in_storage_buffer = 50000


    def get_num_entries(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        num_entries = database.get_row_counts_from_table(t.variable_nts_table_name)
        database.disconnect()

        return num_entries


    def append_entry(self, entry):
        """FIXME This needs documentation to explain difference between append and append_entry"""

        self.db_entries.append(entry)

        if len(self.db_entries) > self.max_num_entries_in_storage_buffer:
            # everytime we are here, the contents of self.db_entries will be stored in the
            # database
            self.store()


    def append(self, entry):
        """Append a single entry based on a sequence

        Parameters
        ==========
        entry : sequence
            values in order they are in the table, entry_id excluded (it will be appended in the
            body of this function)
        """

        self.db_entries.append(entry)
        self.num_entries += 1

        if len(self.db_entries) >= self.max_num_entries_in_storage_buffer:
            # everytime we are here, the contenst of self.db_entries will be stored in the
            # database
            self.store()


    def store(self):
        if not len(self.db_entries):
            return

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''' % t.variable_nts_table_name, self.db_entries)
        database.disconnect()

        if anvio.DEBUG:
            run.info_single("SNVs: %d entries added to the nt variability table." % len(self.db_entries), mc="green")

        self.db_entries = []
