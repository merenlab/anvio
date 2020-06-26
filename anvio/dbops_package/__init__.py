# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

import os
import time

import anvio
import anvio.constants_package as constants_package
import anvio.db as db
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.errors import ConfigError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class Database:
    def __init__(self, db_path, args=None, run=terminal.Run(), progress=terminal.Progress(), quiet=True):
        self.db = None
        self.db_path = db_path
        self.db_type = utils.get_db_type(self.db_path)
        self.db_version = utils.get_required_version_for_db(self.db_path)

        if args:
            A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        else:
            A = anvio.EmptyArgs()
        self.meta_int_keys = A('meta_int_keys') or []
        self.meta_float_keys = A('meta_float_keys') or []
        self.table_info = A('table_info') or []

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):

        if os.path.exists(self.db_path):
            self.db = db.DB(self.db_path, self.db_version)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            for key in self.meta_int_keys:
                try:
                    self.meta[key] = int(self.meta[key])
                except:
                    pass

            for key in self.meta_float_keys:
                try:
                    self.meta[key] = float(self.meta[key])
                except:
                    pass

            self.run.info("%s database', 'An existing database, %s, has been initiated." % (self.db_type, self.db_path), quiet=self.quiet)
        else:
            self.db = None


    def touch(self):
        is_db_ok_to_create(self.db_path, self.db_type)

        self.db = db.DB(self.db_path, self.db_version, new_database=True)

        for table_name, column_names, column_types in self.table_info:
            self.db.create_table(table_name, column_names, column_types)

        return self.db


    def create(self, meta_values={}):
        self.touch()

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

        # know thyself
        self.db.set_meta_value('db_type', self.db_type)

        self.disconnect()

        self.run.info('%s database', 'A new database, %s, has been created.' % (self.db_type, self.db_path), quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


####################################################################################################
#
#     HELPER FUNCTIONS
#
####################################################################################################

def is_db_ok_to_create(db_path, db_type):
    if os.path.exists(db_path):
        raise ConfigError("Anvi'o will not overwrite an existing %s database. "
                          "Please choose a different name or remove the existing database ('%s') first." % (db_type, db_path))

    if not db_path.lower().endswith('.db'):
        raise ConfigError("Please make sure the file name for your new %s db has a '.db' extension. "
                          "Anvi'o developers apologize for imposing their views on how anvi'o databases should be named, "
                          "and are humbled by your cooperation." % db_type)