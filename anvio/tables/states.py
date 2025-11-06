# -*- coding: utf-8
# pylint: disable=line-too-long

import datetime

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table

from anvio.errors import ConfigError

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


class TablesForStates(Table):
    def __init__(self, db_path):
        self.db_path = db_path
        self.states = {}

        if utils.get_db_type(self.db_path) not in ['profile', 'pan', 'structure', 'genes']:
            raise ConfigError("Your database '%s' does not seem to have states table, which anvi'o tries to access.")

        Table.__init__(self, self.db_path, utils.get_required_version_for_db(db_path), run, progress)

        self.init()


    def init(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        self.states = database.get_table_as_dict(t.states_table_name)
        database.disconnect()


    def list_states(self):
        state_names = sorted(list(self.states.keys()))

        self.run.warning('', 'AVAILABLE STATES (%d FOUND)' % (len(self.states)), lc='yellow')
        for state_name in state_names:
            self.run.info_single('%s (last modified on %s)' % (state_name, self.states[state_name]['last_modified']),
                                 nl_after = 1 if state_name == state_names[-1] else 0)


    def get_state(self, state_id):
        if state_id not in self.states:
            return None

        return self.states[state_id]


    def store_state(self, state_id, content, last_modified=None):
        self.remove_state(state_id)

        last_modified = datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S") if not last_modified else last_modified

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec('''INSERT INTO %s VALUES (?,?,?)''' % t.states_table_name, (state_id, content, last_modified))
        self.states = database.get_table_as_dict(t.states_table_name)

        database.disconnect()


    def remove_state(self, state_id):
        self.delete_entries_for_key('name', state_id, [t.states_table_name])
