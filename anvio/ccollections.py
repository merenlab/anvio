# -*- coding: utf-8

"""Implements the collections class (the file name has an extra 'c' to avoid
masking the standard collections library).

If the user have analyzed their metagenome using a metagenome binning software
and identified draft genomes in their data (or by any other means grouped their
contigs based on any criterion), this information can be stored in the
annotation database's collections_* tables. The class implemented here collects
this information from the database, and presents it as an intuitive data structure
for the client.
"""

from itertools import chain

import PaPi.db as db
import PaPi.utils as utils
import PaPi.terminal as terminal

from PaPi.tables import *
from PaPi.utils import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The PaPi Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


def create_blank_collections_tables(db):
    db.create_table(collections_info_table_name, collections_info_table_structure, collections_info_table_types)
    db.create_table(collections_colors_table_name, collections_colors_table_structure, collections_colors_table_types)
    db.create_table(collections_contigs_table_name, collections_contigs_table_structure, collections_contigs_table_types)
    db.create_table(collections_splits_table_name, collections_splits_table_structure, collections_splits_table_types)


class Collections:
    def __init__(self, run = run, progress = progress):
        self.sources_dict = {}


    def populate_sources_dict(self, db_path, version):
        database = db.DB(db_path, version)
        db_type = database.get_meta_value('db_type')
        collections_info_table = database.get_table_as_dict(collections_info_table_name)
        database.disconnect()

        # collections info must be read only if its coming from the annotation database.
        if db_type == 'annotation':
            read_only = True
        elif db_type == 'profile':
            read_only = False
        else:
            raise ConfigError, 'Collections class does not know about this "%s" database type :/' % db_type

        for source in collections_info_table:
            self.sources_dict[source] = collections_info_table[source]
            self.sources_dict[source]['read_only'] = read_only
            self.sources_dict[source]['source_db_path'] = db_path
            self.sources_dict[source]['source_db_version'] = version


    def sanity_check(self, source):
        if source not in self.sources_dict:
            raise ConfigError, 'There is no "%s" I know of. Maybe the populate_sources_dict was not called\
                                for whatever database you are trying to get collections from? (PaPi asks this\
                                rhetorical question to the programmer).'


    def get_collection_dict(self, source):
        self.sanity_check(source)

        c = self.sources_dict[source]

        database = db.DB(c['source_db_path'], c['source_db_version'])
        collections_splits_table = database.get_table_as_dict(collections_splits_table_name)
        database.disconnect()

        # FIXME: this could be resolved with a WHERE clause in the SQL query:
        collection = utils.get_filtered_dict(collections_splits_table, 'source', set([source]))

        collection_dict = {}

        for entry in collection.values():
            source = entry['source']
            cluster_id = entry['cluster_id']
            split = entry['split']

            if collection_dict.has_key(cluster_id):
                collection_dict[cluster_id].append(split)
            else:
                collection_dict[cluster_id] = [split]

        return collection_dict


    def get_collection_colors(self, source):
        self.sanity_check(source)

        c = self.sources_dict[source]

        database = db.DB(c['source_db_path'], c['source_db_version'])
        collections_colors = database.get_table_as_dict(collections_colors_table_name)
        database.disconnect()

        # FIXME: this could be resolved with a WHERE clause in the SQL query:
        collection = utils.get_filtered_dict(collections_colors, 'source', set([source]))

        collection_color_dict = {}

        for entry in collection.values():
            collection_color_dict[entry['cluster_id']] = entry['htmlcolor']

        return collection_color_dict


class TablesForCollections(Table):
    """Populates the collections_* tables, where clusters of contigs and splits are kept"""
    def __init__(self, db_path, version, run=run, progress=progress):
        self.db_path = db_path
        self.version = version

        Table.__init__(self, self.db_path, version, run, progress)

        # set these dudes so we have access to unique IDs:
        self.set_next_available_id(collections_colors_table_name)
        self.set_next_available_id(collections_contigs_table_name)
        self.set_next_available_id(collections_splits_table_name)


    def append(self, source, clusters_dict, cluster_colors = None):
        # remove any pre-existing information for 'source'
        self.delete_entries_for_key('source', source, [collections_info_table_name, collections_contigs_table_name, collections_splits_table_name, collections_colors_table_name])

        num_splits_in_clusters_dict = sum([len(splits) for splits in clusters_dict.values()])
        splits_in_clusters_dict = set(list(chain.from_iterable(clusters_dict.values())))
        if len(splits_in_clusters_dict) != num_splits_in_clusters_dict:
            raise ConfigError, "TablesForCollections::append: %d of the split or contig IDs appear more than once in\
                                your collections input. It is unclear to PaPi how did you manage to do this, but we\
                                cannot go anywhere with this :/" % (num_splits_in_clusters_dict - len(splits_in_clusters_dict))

        database = db.DB(self.db_path, self.version)

        # how many clusters are defined in 'clusters_dict'?
        cluster_ids = clusters_dict.keys()

        # push information about this search result into serach_info table.
        db_entries = tuple([source, num_splits_in_clusters_dict, len(cluster_ids)])
        database._exec('''INSERT INTO %s VALUES (?,?,?)''' % collections_info_table_name, db_entries)

        # populate colors table.
        if not cluster_colors:
            cluster_colors = utils.get_random_colors_dict(cluster_ids)
        db_entries = [(self.next_id(collections_colors_table_name), source, cid, cluster_colors[cid]) for cid in cluster_ids]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % collections_colors_table_name, db_entries)

        # populate splits table
        db_entries = []
        for cluster_id in clusters_dict:
            for split_name in clusters_dict[cluster_id]:
                db_entries.append(tuple([self.next_id(collections_splits_table_name), source, split_name, cluster_id]))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % collections_splits_table_name, db_entries)


        # FIXME: This function can be called to populate the annotation database (via papi-populate-collections), or
        # the profile database. when it is annotation database, the superclass Table has the self.splits_info variable
        # set when it is initialized. however, the Table instance is missing self.splis when it is initialized with
        # the profile database. hence some special controls for annotation db (note that collections_contigs_table is
        # only populated in the annotations database):
        if self.db_type == 'annotation':
            splits_only_in_clusters_dict = [c for c in splits_in_clusters_dict if c not in self.splits_info]
            splits_only_in_db = [c for c in self.splits_info if c not in splits_in_clusters_dict]

            if len(splits_only_in_clusters_dict):
                self.run.warning('%d of %d splits found in "%s" results are not in the database. This may be OK,\
                                          but you must be the judge of it. If this is somewhat surprising, please use caution\
                                          and make sure all is fine before going forward with you analysis.'\
                                                % (len(splits_only_in_clusters_dict), len(splits_in_clusters_dict), source))

            if len(splits_only_in_db):
                self.run.warning('%d of %d splits found in the database were missing from the "%s" results. If this\
                                          does not make any sense, please make sure you know why before going any further.'\
                                                % (len(splits_only_in_db), len(self.splits_info), source))

            # then populate contigs table.
            db_entries = self.process_contigs(source, clusters_dict)
            database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % collections_contigs_table_name, db_entries)

        database.disconnect()

        self.run.info('Collections', '%s annotations for %d splits have been successfully added to the database at "%s".'\
                                        % (source, len(db_entries), self.db_path), mc='green')


    def process_contigs(self, source, clusters_dict):
        db_entries_for_contigs = []

        split_to_cluster_id = {}
        for cluster_id in clusters_dict:
            for split_name in clusters_dict[cluster_id]:
                split_to_cluster_id[split_name] = cluster_id

        contigs_processed = set([])
        for split_name in split_to_cluster_id:
            if split_name not in self.splits_info:
                # which means this split only appears in the input file, but not in the database.
                continue

            contig_name = self.splits_info[split_name]['parent']

            if contig_name in contigs_processed:
                continue
            else:
                contigs_processed.add(contig_name)

            db_entry = tuple([self.next_id(collections_contigs_table_name), source, contig_name, split_to_cluster_id[split_name]])
            db_entries_for_contigs.append(db_entry)

        return db_entries_for_contigs
