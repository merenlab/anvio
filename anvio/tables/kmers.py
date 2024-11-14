# -*- coding: utf-8
# pylint: disable=line-too-long

"""TablesForCollections"""


import anvio
import anvio.kmers as kmers
import anvio.terminal as terminal


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


class KMerTablesForContigsAndSplits:
    def __init__(self, table_name, k=4):
        self.table_name = table_name
        self.kmers_class = kmers.KMers(k)
        self.kmers = sorted(list(self.kmers_class.kmers[k]))

        self.kmer_dict = {}
        self.db_entries = []

        self.kmers_table_structure = ["contig"] + self.kmers
        self.kmers_table_types = ["text"] + ["numeric"] * len(self.kmers)

    def get_kmer_freq(self, sequence):
        return self.kmers_class.get_kmer_frequency(sequence, dist_metric_safe=True)

    def append(self, seq_id, sequence, kmer_freq=None):
        if not kmer_freq:
            kmer_freq = self.kmers_class.get_kmer_frequency(
                sequence, dist_metric_safe=True
            )

        db_entry = tuple([seq_id] + [kmer_freq[kmer] for kmer in self.kmers])
        self.db_entries.append(db_entry)

    def store(self, db):
        db.create_table(
            self.table_name, self.kmers_table_structure, self.kmers_table_types
        )
        db._exec_many(
            """INSERT INTO %s VALUES (%s)"""
            % (self.table_name, (",".join(["?"] * len(self.kmers_table_structure)))),
            self.db_entries,
        )
