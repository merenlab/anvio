# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.tables as t
import anvio.terminal as terminal

from anvio.dbops import DBClassFactory
from anvio.errors import ConfigError
from anvio.tables.tableops import Table
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


class TablesForViews(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.db_path = db_path

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), self.run, self.progress)


    def sanity_check(self, view_data):
        if not isinstance(view_data, list):
            self.progress.reset()
            raise ConfigError(f"View data must be a list of tuples or list of lists :( Yours is a {type(view_data)}.")

        if len(view_data[0]) != 3:
            self.progress.reset()
            raise ConfigError(f"Each item in the view data list must be a list or tuple with three items (item name, "
                              f"layer name, and data value). The length of yours is {len(view_data[0])}. Not OK.")

        data_tuple_lengths = set([len(d) for d in view_data])
        if len(data_tuple_lengths) > 1:
            self.progress.reset()
            raise ConfigError("Each list or tuple item in the view data list must be of length 3. Your view data is "
                              "composed of mixed lengths :/")

        try:
            [float(d[2]) for d in view_data]
        except:
            self.progress.reset()
            raise ConfigError("The third item in every list or tuple in your view data list must be a numerical value. "
                              "At lesat that is what anvi'o asks for currently. We can get rid of this requirement, but "
                              "that is another discussion. For now, you need to ensure your view data complies with what "
                              "anvi'o expects :/")


    def matrix_to_long_form(self, view_data):
        d = []
        for item in view_data:
            for layer in view_data[item]:
                d.append((item, layer, view_data[item][layer]), )

        return(d)


    def create_new_view(self, view_data, table_name, view_name=None, append_mode=False, skip_sanity_check=False, from_matrix_form=False):
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

        if from_matrix_form:
            view_data = self.matrix_to_long_form(view_data)

        if not skip_sanity_check:
            self.sanity_check(view_data)

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
            anvio_db.db.create_table(table_name, t.view_table_structure, t.view_table_types)
        except Exception as e:
            # FIXME: the following if statement will omit errors and quietly continue despite the
            # table creation failed. I think we should remove it, and add `create_table` function
            # a new flag, such as `ok_if_exists` and call it in this context as
            # `ok_if_exists=append_mode`.
            if not append_mode:
                raise ConfigError("Something bad happened when anvi'o was trying to create table `%s` in database "
                                  "'%s'. Here is how the part of the code that was about this described the "
                                  "problem: '%s'." % (table_name, self.db_path, str(e)))

        anvio_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (table_name, ','.join(['?'] * len(t.view_table_structure))), view_data)

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
