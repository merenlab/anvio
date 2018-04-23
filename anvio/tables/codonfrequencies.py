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


    def append(self, profile):
        db_entry = tuple([self.next_id(t.variable_codons_table_name)] + [profile[h] for h in t.variable_codons_table_structure[1:]])
        self.db_entries.append(db_entry)
        self.num_entries += 1
        if self.num_entries % 100 == 0:
            self.progress.update('Information for %d codons have been added ...' % self.num_entries)


    def store(self):
        utils.is_profile_db(self.db_path)
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.variable_codons_table_name, ','.join(['?'] * len(t.variable_codons_table_structure))), self.db_entries)
        database.disconnect()



