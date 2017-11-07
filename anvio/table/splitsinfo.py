import anvio.utils as utils
import anvio.tables as t

class InfoTableForSplits(object):
    def __init__(self, db):
        self.db = db
        self.db_entries = []
        self.total_splits = 0


    def append(self, seq_id, sequence, order, start, end, parent_gc_content, parent):
        self.total_splits += 1
        sequence_length = len(sequence)
        db_entry = tuple([seq_id, order, start, end, sequence_length, utils.get_GC_content_for_sequence(sequence), parent_gc_content, parent])
        self.db_entries.append(db_entry)


    def store(self):
        self.db.insert_many(t.splits_info_table_name, self.db_entries)

