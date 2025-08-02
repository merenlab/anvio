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


class TableForNodes(Table):
    """A class to populate pangenome graph nodes in a given pan-graph-db.

      Here is an example:

        >>> table = TableForNodes(db_path)
        >>> for ...:
        >>>     table.add({'syn_gene_cluster_id': syn_gene_cluster_id, 'syn_gene_cluster_type': syn_gene_cluster_type, 'gene_cluster_id': gene_cluster_id, 'gene_calls_json': gene_calls_json})
        >>> table.store()
    """

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        utils.is_pan_graph_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pangraph__version__, run, progress)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in t.pan_graph_nodes_table_structure])


    def store(self):
        self.delete_contents_of_table(t.pan_graph_nodes_table_name, warning=False)

        db_entries = [tuple(entry) for entry in self.entries]

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.pan_graph_nodes_table_name, db_entries)
        database.disconnect()


class TableForEdges(Table):
    """A class to populate pangenome graph edges in a given pan-graph-db"""

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        utils.is_pan_graph_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pangraph__version__, run, progress)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in t.pan_graph_edges_table_structure])


    def store(self):
        self.delete_contents_of_table(t.pan_graph_edges_table_name, warning=False)

        db_entries = [tuple(entry) for entry in self.entries]

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.pan_graph_edges_table_name, db_entries)
        database.disconnect()


