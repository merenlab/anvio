# -*- coding: utf-8
# pylint: disable=line-too-long

"""Implements the collections class (the file name has an extra 'c' to avoid
masking the standard collections library).

If the user have analyzed their metagenome using a metagenome binning software
and identified draft genomes in their data (or by any other means binned their
contigs based on any criterion), this information can be stored in the
contigs database's collections_* tables. The class implemented here collects
this information from the database, and presents it as an intuitive data structure
for the client.
"""

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

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
    def __init__(self, r=run, p=progress):
        self.collections_dict = {}
        self.run = r
        self.progress = p


    def populate_collections_dict(self, db_path, version):
        database = db.DB(db_path, version)
        db_type = database.get_meta_value('db_type')
        collections_info_table = database.get_table_as_dict(t.collections_info_table_name)
        database.disconnect()

        # collections info must be read only if its coming from the contigs database.
        if db_type == 'contigs':
            read_only = True
        elif db_type == 'profile':
            read_only = False
        else:
            raise ConfigError, 'Collections class does not know about this "%s" database type :/' % db_type

        for collection_name in collections_info_table:
            self.collections_dict[collection_name] = collections_info_table[collection_name]
            self.collections_dict[collection_name]['read_only'] = read_only
            self.collections_dict[collection_name]['source_db_path'] = db_path
            self.collections_dict[collection_name]['source_db_version'] = version


    def sanity_check(self, collection_name):
        if collection_name not in self.collections_dict:
            raise ConfigError, 'There is no "%s" I know of. Probably something is spelled wrong somewhere? In case you are\
                                a programmer and accessing to the collections from your program, here is a reminder for you:\
                                are you sure `populate_collections_dict` was called for whatever database you are trying to\
                                get collections from?' % collection_name


    def get_collection_dict(self, collection_name):
        self.sanity_check(collection_name)

        c = self.collections_dict[collection_name]

        database = db.DB(c['source_db_path'], c['source_db_version'])
        collections_splits_table = database.get_table_as_dict(t.collections_splits_table_name)
        database.disconnect()

        # FIXME: this could be resolved with a WHERE clause in the SQL query:
        collection_dict_from_db = utils.get_filtered_dict(collections_splits_table, 'collection_name', set([collection_name]))

        collection_dict_to_return = {}

        for entry in collection_dict_from_db.values():
            collection_name = entry['collection_name']
            bin_name = entry['bin_name']
            split = entry['split']

            if bin_name in collection_dict_to_return:
                collection_dict_to_return[bin_name].append(split)
            else:
                collection_dict_to_return[bin_name] = [split]

        return collection_dict_to_return


    def get_bins_info_dict(self, collection_name):
        self.sanity_check(collection_name)

        c = self.collections_dict[collection_name]

        database = db.DB(c['source_db_path'], c['source_db_version'])
        collections_bins_info_table = database.get_table_as_dict(t.collections_bins_info_table_name)
        database.disconnect()

        # FIXME: this could be resolved with a WHERE clause in the SQL query:
        collections_bins_info_table_filtered = utils.get_filtered_dict(collections_bins_info_table, 'collection_name', set([collection_name]))

        bins_info_dict = {}
        for v in collections_bins_info_table_filtered.values():
            bins_info_dict[v['bin_name']] = {'html_color': v['html_color'], 'source': v['source']}

        return bins_info_dict


    def list_collections(self):
        self.run.warning('', 'COLLECTIONS FOUND', lc='yellow')
        for collection_name in self.collections_dict:
            c = self.collections_dict[collection_name]
            output = '%s (%d bins, representing %d splits).' % (collection_name, c['num_bins'], c['num_splits'])
            self.run.info_single(output)


    def export_collection(self, collection_name, output_file_prefix=None):
        self.sanity_check(collection_name)

        if not output_file_prefix:
            output_file_prefix = 'collection-%s' % (collection_name.strip().replace(' ', '-'))

        info_file_path = output_file_prefix + '-info.txt'
        items_file_path = output_file_prefix + '.txt'

        self.run.info('Items file path', items_file_path)
        filesnpaths.is_output_file_writable(items_file_path)

        bins_info = self.get_bins_info_dict(collection_name)
        collection = self.get_collection_dict(collection_name)

        if len(bins_info):
            self.run.info('Info file path', info_file_path)
            info_file = open(info_file_path, 'w')
            for bin_name in bins_info:
                info_file.write('%s\t%s\t%s\n' % (bin_name, bins_info[bin_name]['source'], bins_info[bin_name]['html_color']))
            info_file.close()

        items_file = open(items_file_path, 'w')
        for bin_name in collection:
            for item_name in collection[bin_name]:
                items_file.write('%s\t%s\n' % (item_name, bin_name))


class GetSplitNamesInBins:
    def __init__(self, args):
        # we will fill this in and return it
        self.split_names_of_interest = set([])

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bin_ids_file_path = A('bin_ids_file')
        self.bin_id = A('bin_id')
        self.collection_name = A('collection_name')
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.debug = A('debug')

        if self.bin_ids_file_path and self.bin_id:
            raise ConfigError, 'Either use a file to list all the bin ids (-B), or declare a single bin (-b)\
                                you would like to focus. Not both :/'
        if (not self.bin_ids_file_path) and (not self.bin_id):
            raise ConfigError, "You must either use a file to list all the bin ids (-B) you would like to\
                                focus on, or declare a single bin id (-b) from your collection. You have\
                                not really given anvi'o anything to work with."

        if not self.collection_name:
            raise ConfigError, 'This will not work without a collection ID for your bins :/'

        if self.bin_ids_file_path:
            filesnpaths.is_file_exists(self.bin_ids_file_path)
            self.bins = set([b.strip() for b in open(self.bin_ids_file_path).readlines()])
        if self.bin_id:
            self.bins = set([self.bin_id])

        if not len(self.bins):
            raise ConfigError, 'There is no bin to work with :/'

        self.collections = Collections()
        self.collections.populate_collections_dict(self.profile_db_path, anvio.__profile__version__)

        if self.collection_name not in self.collections.collections_dict:
            raise ConfigError, 'The collection id "%s" does not seem to be in the profile database. These are the\
                                collections that are available through this profile database: %s.'\
                                                    % (self.collection_name, ', '.join(self.collections.collections_dict))

        self.collection_dict = self.collections.get_collection_dict(self.collection_name)

        bins_in_collection = self.collection_dict.keys()

        bins_that_does_not_exist_in_collection = [b for b in self.bins if b not in bins_in_collection]
        if len(bins_that_does_not_exist_in_collection):
            raise ConfigError, 'Some of the bins you requested does not appear to have been described in the collection\
                                "%s". Here is a list of bins that are missing: %s'\
                                        % (self.collection_name, ', '.join(bins_that_does_not_exist_in_collection))



    def get_split_names_only(self):
        split_names_of_interest = []
        for bin_id in self.bins:
            split_names_of_interest.extend(self.collection_dict[bin_id])

        self.split_names_of_interest = set(split_names_of_interest)

        return self.split_names_of_interest


    def get_dict(self):
        d = {}

        for bin_id in self.bins:
            d[bin_id] = set(self.collection_dict[bin_id])

        return d


class GetSequentialBlocksOfSplits:
    """A simple class to identify longest stretches in a list of integers.

       >>> sequentials = SequentialBlocksOfSplits([1, 2, 3, 5, 6, 9])
       >>> print sequentials.blocks
           [[1, 2, 3], [5, 6], [9]]
       >>>
    """

    def __init__(self, l):
        self.l = sorted(list(set(l)))
        self.blocks = []
        self.current_block = []


    def finalize_block(self):
        self.blocks.append(self.current_block)
        self.current_block = []


    def process(self):
        while True:
            if not self.l:
                break

            current = self.l.pop(0)

            if not len(self.current_block) or current == self.current_block[-1] + 1:
                self.current_block.append(current)
            else:
                self.finalize_block()
                self.current_block.append(current)

        self.finalize_block()

        return self.blocks

