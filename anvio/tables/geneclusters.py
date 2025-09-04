# -*- coding: utf-8
# pylint: disable=line-too-long

"""TablesForCollections"""

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal

from anvio.tables.tableops import Table
from anvio.dbinfo import is_pan_db
from anvio.utils.database import get_required_version_for_db


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


class TableForGeneClusters(Table):
    """A class to populate gene clusters table in a given pan db.

      Here is an example:

        >>> table = TableForGeneClusters(db_path)
        >>> for ...:
        >>>     table.add({'gene_caller_id': gene_caller_id, 'gene_cluster_id': gene_cluster_id, 'genome_name': genome_name})
        >>> table.store()
    """

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        is_pan_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pan__version__, run, progress)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in t.pan_gene_clusters_table_structure])


    def store(self):
        self.delete_contents_of_table(t.pan_gene_clusters_table_name, warning=False)

        db_entries = [tuple(entry) for entry in self.entries]

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.pan_gene_clusters_table_name, db_entries)
        database.disconnect()


