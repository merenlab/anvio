import anvio.tables as t

class TablesForStates(object):
    def __init__(self, db):
        self.db = db
        self.states = {}

        self.load_states()


    def load_states(self):
        self.states = self.db.get_table_as_dict(t.states_table_name)


    def list_states(self):
        state_names = sorted(list(self.states.keys()))

        self.run.warning('', 'AVAILABLE STATES (%d FOUND)' % (len(self.states)), lc='yellow')
        for state_name in state_names:
            self.run.info_single('%s (last modified on %s)' % (state_name, self.states[state_name]['last_modified']),
                                 nl_after = 1 if state_name == state_names[-1] else 0)


    def get_state(self, state_name):
        if state_name not in self.states:
            return None

        return self.states[state_id]


    def store_state(self, state_id, content, last_modified=None):
        self.remove_state(state_id)

        if not last_modified:
            last_modified = datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S")

        self.db.insert(t.states_table_name,(state_id, content, last_modified, ))
        self.load_states()



    def remove_state(self, state_id):
        self.db.delete_entries_for_key('name', state_id, [t.states_table_name])
        self.load_states()