# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.dbops import DBClassFactory
from anvio.errors import ConfigError
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


class TablesForViews(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, utils.get_required_version_for_db(db_path), run, progress)


    def create_new_view(self, data_dict, table_name, table_structure, table_types, view_name=None, append_mode=False):
        """Creates a new view table, and adds an entry for it into the 'views' table.

        Entries in 'views' table appear in various places in the interface. However, we also generate
        view tables to store the type of data we do not wish to display on interfaces, but be able
        access from various other modules. A good example to this is the item_order recipes. When we
        profile a sample, we treat every stplit as their own entity with respect to their mean coverage.
        Although it is great for visualization purposes, it is not useful for item_order purposes since in
        most cases we wish splits to stay together in item_order output. Hence, we create a mean_coverage_splits
        table, where each split holds their own coverage, and we create a mean_coverage_contigs table where each
        split has the coverage of their parent. Clearly the second table is not useful to display. When a table
        is not added as an entry to the 'views' table, then it only exists in the database for other purposes
        than displaying it.

        If a new view does not have a 'view_id', it is not added the 'views' table to provide that flexibility.
        """

        anvio_db = DBClassFactory().get_db_object(self.db_path)

        views_in_db = anvio_db.db.get_table_as_dict(t.views_table_name)

        if not append_mode:
            if view_name and view_name in views_in_db:
                raise ConfigError("TablesForViews speaking: Yo yo yo. You already have a view in the db called '%s'.\
                                    You can't create another one before you get rid of the existing one, because rules."\
                                                                            % view_name)

            # first create the data table:
            anvio_db.db.drop_table(table_name)

        try:
            anvio_db.db.create_table(table_name, table_structure, table_types)
        except:
            if not append_mode:
                raise ConfigError("Table already exists")

        db_entries = [tuple([item] + [data_dict[item][h] for h in table_structure[1:]]) for item in data_dict]
        anvio_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (table_name, ','.join(['?'] * len(table_structure))), db_entries)

        if view_name and view_name not in views_in_db:
            anvio_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % t.views_table_name, (view_name, table_name))

        anvio_db.disconnect()


    def remove(self, view_name, table_names_to_blank=[]):
        anvio_db = DBClassFactory().get_db_object(self.db_path)
        anvio_db.db._exec('''DELETE FROM %s WHERE view_id = "%s"''' % (t.views_table_name, view_name))
        for table_name in table_names_to_blank:
            if table_name in anvio_db.db.get_table_names():
                anvio_db.db._exec('''DELETE FROM %s''' % table_name)
        anvio_db.disconnect()
