# # -*- coding: utf-8
# # pylint: disable=line-too-long


# import anvio
# import anvio.tables as t
# import anvio.utils as utils
# import anvio.terminal as terminal


# __author__ = "Developers of anvi'o (see AUTHORS.txt)"
# __copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
# __credits__ = []
# __license__ = "GPL 3.0"
# __version__ = anvio.__version__
# __maintainer__ = "A. Murat Eren"
# __email__ = "a.murat.eren@gmail.com"
# __status__ = "Development"


# run = terminal.Run()
# progress = terminal.Progress()
# pp = terminal.pretty_print


# class TableForTRNASeedsInfo:
#     def __init__(self):
#         self.db_entries = []
#         self.total_nts = 0
#         self.total_trnaseeds = 0


#     def append(self, seq_id, sequence, gene_start_stops=None):
#         sequence_length = len(sequence)

#         self.total_nts += sequence_length
#         self.total_trnaseeds += 1
#         db_entry = tuple([seq_id, sequence_length])
#         self.db_entries.append(db_entry)

#         return (sequence_length, )


#     def store(self, db):
#         if len(self.db_entries):
#             db._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.trnaseeds_info_table_name, (','.join(['?'] * len(self.db_entries[0])))), self.db_entries)
