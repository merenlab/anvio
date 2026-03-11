# -*- coding: utf-8
# pylint: disable=line-too-long

"""TablesForCollections"""

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


class TableForGeneClusters(Table):
    """A class to populate gene clusters table in a given pan db.

    Here is an example:

        >>> table = TableForGeneClusters(db_path)
        >>> for ...:
        >>>     table.add({'gene_caller_id': gene_caller_id, 'gene_cluster_id': gene_cluster_id, 'genome_name': genome_name})
        >>> table.store()

    Please note the very important `gc_tracker_table` parameter. When this paramter is set to True,
    the `store` function will no longer use the default table name for gene cluster storage, but an alternative
    table to keep track of sequence-based gene clusters in structure mode.
    """

    def __init__(self, db_path, gc_tracker_table=False, run=run, progress=progress):
        self.db_path = db_path

        if gc_tracker_table:
            self.table_name, self.table_structure = t.pan_gc_tracker_table_name, t.pan_gc_tracker_table_structure
        else:
            self.table_name, self.table_structure = t.pan_gene_clusters_table_name, t.pan_gene_clusters_table_structure

        utils.is_pan_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pan__version__, run, progress)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in self.table_structure])


    def store(self):
        self.delete_contents_of_table(t.pan_gene_clusters_table_name, warning=False)

        db_entries = [tuple(entry) for entry in self.entries]

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % self.table_name, db_entries)
        database.disconnect()



class TableForPSGCGCAssociations(Table):
    """A class to populate protein structure-informed gene cluster <-> de novo gene cluster associations in a pan db.

      Here is an example:

        >>> table = TableForPSGCGCAssociations(db_path)
        >>> for ...:
        >>>     table.add({'gene_cluster_id': gene_cluster_id, 'protein_structure_informed_gene_cluster_id': psgc_id})
        >>> table.store()
    """

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        utils.is_pan_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pan__version__, run, progress)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in t.pan_gc_psgc_associations_table_structure])


    def store(self):
        self.delete_contents_of_table(t.pan_gc_psgc_associations_table_name, warning=False)

        db_entries = [tuple(entry) for entry in self.entries]

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.pan_gc_psgc_associations_table_name, db_entries)
        database.disconnect()


