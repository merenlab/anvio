# -*- coding: utf-8
'''Storing and retrieving metadata regarding contigs and splits'''

import anvio.tables as t

import anvio
from anvio.terminal import Progress
from anvio.terminal import pretty_print as pp


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


progress = Progress()
progress.verbose = False


class Metadata:
    def __init__(self, p=progress):
        self.metadata_contigs = {}
        self.metadata_splits = {}
        self.progress = p


    def store_metadata_for_contigs_and_splits(self, sample_id, contigs, db):
        self.progress.new('Storing metadata')

        num_contigs = pp(len(contigs))
        cur_contig = 1

        # this loop will get metadata information from Contig instanes and store them into the db
        # at once. this was broken down into about 10 functions, but this structure seems to be the most efficient
        # although it looks crappy:
        for contig_name in contigs:
            self.progress.update("Processing contig %s of %s" % (pp(cur_contig), num_contigs))
            contig = contigs[contig_name]
            contig_metadata = contig.get_metadata_dict()

            self.metadata_contigs[contig.name] = {'contig': contig.name}
            for metadata_field in t.metadata_table_structure[1:]:
                self.metadata_contigs[contig.name][metadata_field] = contig_metadata[metadata_field]

            # contig is done, deal with splits in it:
            for split in contig.splits:
                split_metadata = split.get_metadata_dict()
                self.metadata_splits[split.name] = {'contig': split.name}
                for metadata_field in t.metadata_table_structure[1:]:
                    self.metadata_splits[split.name][metadata_field] = split_metadata[metadata_field]


        self.progress.update("Generating tables ...")
        gen_metadata_tables(self.metadata_splits, self.metadata_contigs, db)
        self.progress.end()



def gen_metadata_tables(metadata_splits, metadata_contigs, db):
    # all objects are ready, creating tables next.
    db.create_table('metadata_splits', t.metadata_table_structure, t.metadata_table_types)
    db_entries = [tuple([split] + [metadata_splits[split][h] for h in t.metadata_table_structure[1:]]) for split in metadata_splits]
    db._exec_many('''INSERT INTO metadata_splits VALUES (?,?,?,?,?,?,?,?,?,?)''', db_entries)

    db.create_table('metadata_contigs', t.metadata_table_structure, t.metadata_table_types)
    db_entries = [tuple([split] + [metadata_contigs[metadata_splits[split]['__parent__']][h] for h in t.metadata_table_structure[1:]]) for split in metadata_splits]
    db._exec_many('''INSERT INTO metadata_contigs VALUES (?,?,?,?,?,?,?,?,?,?)''', db_entries)

    db.commit()

