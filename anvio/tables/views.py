# -*- coding: utf-8
# pylint: disable=line-too-long

import os

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
        profile a sample, we treat every split as its own entity with respect to its mean coverage.
        Although it is great for visualization purposes, it is not useful for item_order purposes since in
        most cases we wish splits to stay together in item_order output. Hence, we create a mean_coverage_splits
        table, where each split holds its own coverage, and we create a mean_coverage_contigs table where each
        split has the coverage of its parent. Clearly the second table is not useful to display. When a table
        is not added as an entry to the 'views' table, then it only exists in the database for other purposes
        than displaying it.

        If a new view does not have a 'view_id', it is not added the 'views' table to provide that flexibility.
        """

        anvio_db = DBClassFactory().get_db_object(self.db_path)

        views_in_db = anvio_db.db.get_table_as_dict(t.views_table_name)

        if not append_mode:
            if view_name and view_name in views_in_db:
                raise ConfigError("TablesForViews speaking: Yo yo yo. You already have a view in the db '%s' called '%s'. "
                                   "You can't create another one before you get rid of the existing one, because rules."\
                                                                            % (self.db_path, view_name))

            # first create the data table:
            anvio_db.db.drop_table(table_name)

        try:
            anvio_db.db.create_table(table_name, table_structure, table_types)
        except Exception as e:
            # FIXME: the following if statement will omit errors and quietly continue despite the
            # table creation failed. I think we should remove it, and add `create_table` function
            # a new flag, such as `ok_if_exists` and call it in this context as
            # `ok_if_exists=append_mode`.
            if not append_mode:
                raise ConfigError("Something bad happened when anvi'o was trying to create table `%s` in database "
                                  "'%s'. Here is how the part of the code that was about this described the "
                                  "problem: '%s'." % (table_name, self.db_path, str(e)))

        db_entries = [tuple([item] + [data_dict[item][h] for h in table_structure[1:]]) for item in data_dict]

        try:
            anvio_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (table_name, ','.join(['?'] * len(table_structure))), db_entries)
        except Exception as e:
            num_columns = set([len(x) for x in db_entries])
            if len(num_columns) == 1:
                columns_text = "%d columns" % num_columns.pop()
            else:
                columns_text = "%d to %d columns (which is utterly weird)" % (min(num_columns), max(num_columns))

            temp_file_output_path = os.path.abspath(self.db_path) + '-DB_ENTRIES_FOR_SAD_ERROR.txt'
            with open(temp_file_output_path, 'w') as temp_file:
                for entry in db_entries:
                    temp_file.write('\t'.join([str(e) for e in entry]) + '\n')

            raise ConfigError("Something bad happened while anvi'o was trying to insert %d entries with %s into the "
                              "table '%s' which contained a table structure with %d columns in '%s' :( This "
                              "is the error we got back from the database module: \"%s\". Anvi'o created a temporary "
                              "file for you so you can see the contents of the db_entries it tried to add to the "
                              "database, which is here: '%s'." % \
                                    (len(db_entries), columns_text, table_name, len(table_structure), self.db_path, e, temp_file_output_path))

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
