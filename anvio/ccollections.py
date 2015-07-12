# -*- coding: utf-8

"""Implements the collections class (the file name has an extra 'c' to avoid
masking the standard collections library).

If the user have analyzed their metagenome using a metagenome binning software
and identified draft genomes in their data (or by any other means binned their
contigs based on any criterion), this information can be stored in the
annotation database's collections_* tables. The class implemented here collects
this information from the database, and presents it as an intuitive data structure
for the client.
"""


import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class Collections:
    def __init__(self, r = run, p = progress):
        self.sources_dict = {}
        self.run = r
        self.progress = p


    def populate_sources_dict(self, db_path, version):
        database = db.DB(db_path, version)
        db_type = database.get_meta_value('db_type')
        collections_info_table = database.get_table_as_dict(t.collections_info_table_name)
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
                                for whatever database you are trying to get collections from? (anvio asks this\
                                rhetorical question to the programmer).'


    def get_collection_dict(self, source):
        self.sanity_check(source)

        c = self.sources_dict[source]

        database = db.DB(c['source_db_path'], c['source_db_version'])
        collections_splits_table = database.get_table_as_dict(t.collections_splits_table_name)
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
        collections_colors = database.get_table_as_dict(t.collections_colors_table_name)
        database.disconnect()

        # FIXME: this could be resolved with a WHERE clause in the SQL query:
        collection = utils.get_filtered_dict(collections_colors, 'source', set([source]))

        collection_color_dict = {}

        for entry in collection.values():
            collection_color_dict[entry['cluster_id']] = entry['htmlcolor']

        return collection_color_dict


    def list_collections(self):
        for collection_id in self.sources_dict:
            c = self.sources_dict[collection_id]
            output = '%s (%d clusters, representing %d splits).' % (collection_id, c['num_clusters'], c['num_splits'])
            self.run.info_single(output)

