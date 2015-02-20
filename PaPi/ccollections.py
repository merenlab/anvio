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

from collections import Counter

import PaPi.fastalib as u
import PaPi.utils as utils
import PaPi.terminal as terminal
import PaPi.annotation as annotation

from PaPi.utils import ConfigError
from PaPi.filesnpaths import FilesNPathsError


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


class Collections:
    def __init__(self, annotation_db_path, source = None, run = run, progress = progress):
        self.collections_dict = {}

        # hi db
        annotation_db = annotation.AnnotationDatabase(annotation_db_path)

        # read info table to get what is available in the db
        collections_info_table = annotation_db.db.get_table_as_dict(annotation.collections_info_table_name)
        self.sources = collections_info_table.keys()

        # read search table (which holds hmmscan hits for splits).
        collections_splits_table = annotation_db.db.get_table_as_dict(annotation.collections_splits_table_name)

        # we're done with the db
        annotation_db.disconnect()

        if source:
            if source not in self.sources:
                raise ConfigError, 'Source "%s" is not one of the binning sources found in the database.' % source

            # filter out sources that are not requested
            self.sources = [source]
            collections_splits_table = utils.get_filtered_dict(collections_splits_table, 'source', set([source]))

        for source in self.sources:
            self.collections_dict[source] = {}

        for entry in collections_splits_table.values():
            source = entry['source']
            cluster_id = entry['cluster_id']
            split = entry['split']

            if self.collections_dict[source].has_key(cluster_id):
                self.collections_dict[source][cluster_id].append(split)
            else:
                self.collections_dict[source][cluster_id] = [split]


    def get_collections_dict(self):
        return self.collections_dict
