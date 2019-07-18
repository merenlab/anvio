# -*- coding: utf-8
# pylint: disable=line-too-long

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table
from anvio.errors import ConfigError
from anvio.tables.genecalls import GenesInSplits


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



class TaxonNamesTable(object):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress


    def populate_taxon_names_table(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        db_entries = [tuple([t_name_id] + [self.taxon_names_dict[t_name_id][t_level] for t_level in t.taxon_names_table_structure[1:]]) for t_name_id in self.taxon_names_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.taxon_names_table_name, db_entries)

        database.disconnect()
        self.run.info('Taxon names table', 'Updated with %d unique taxon names' % len(db_entries))
