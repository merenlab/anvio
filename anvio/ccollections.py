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

import copy

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tables.collections import TablesForCollections
from anvio.utils.database import get_all_item_names_from_the_database, get_required_version_for_db
from anvio.utils.misc import get_filtered_dict


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


class Collections:
    def __init__(self, r=run, p=progress):
        self.collections_dict = {}
        self.run = r
        self.progress = p

        self.db_type = None
        self.db_path = None


    def populate_collections_dict(self, db_path):
        filesnpaths.is_file_exists(db_path)
        self.db_path = db_path

        database = db.DB(db_path, get_required_version_for_db(db_path))
        self.db_type = database.get_meta_value('db_type')
        collections_info_table = database.get_table_as_dict(t.collections_info_table_name)
        database.disconnect()

        # collections info must be read only if its coming from the contigs database.
        if self.db_type == 'contigs':
            read_only = True
        elif self.db_type == 'profile':
            read_only = False
        elif self.db_type:
            read_only = False
        elif self.db_type == 'pan':
            read_only = False
        else:
            raise ConfigError('Collections class does not know about this "%s" database type :/' % self.db_type)

        for collection_name in collections_info_table:
            self.collections_dict[collection_name] = collections_info_table[collection_name]
            self.collections_dict[collection_name]['read_only'] = read_only
            self.collections_dict[collection_name]['source_db_path'] = db_path
            self.collections_dict[collection_name]['source_db_version'] = get_required_version_for_db(db_path)


    def sanity_check(self, collection_name):
        if collection_name not in self.collections_dict:
            raise ConfigError('There is no "%s" I know of. Probably something is spelled wrong somewhere? In case you are '
                              'a programmer and accessing to the collections from your program, here is a reminder for you: '
                              'are you sure `populate_collections_dict` was called for whatever database you are trying to '
                              'get collections from? If you are a user, you can always try to use the `--list-collections` '
                              'flag and hope for the best.' % collection_name)


    def get_trimmed_dicts(self, collection_name, split_names = set([])):
        """Returns collection_dict, bins_info_dict for splits matching split_names, and
        split names that are in the db, but not binned in the collection..

        Any bin that does not have any splits left after removal simply is removed
        from the dictionary"""

        self.progress.new('Recovering collection information for "%s" ...' % collection_name)

        collection_dict = self.get_collection_dict(collection_name)
        bins_info_dict = self.get_bins_info_dict(collection_name)

        self.progress.update('Identifying split names that are in the profile db but do not appear in any bin ...')
        split_names_in_db_but_missing_in_collection = copy.deepcopy(split_names)
        for bin_id in collection_dict:
            split_names_in_db_but_missing_in_collection -= set(collection_dict[bin_id])

        self.progress.update('Identifying bin names that do not have any splits that appear in the profile database ...')
        bins_with_zero_splits_in_profile_db = []
        bin_ids_in_collection = list(collection_dict.keys())
        for bin_id in bin_ids_in_collection:
            # good split names are the ones that appear in `split_names` user sent. so here we will replace
            # the content of each bin with only split names that are 'good' in that sense. in practice, this
            # will ensure that the collection dict will not contain any split name that does not appear in
            # the profile database (i.e., a relevant need can be seen in the `load_collection_mode` function
            # in the interactive.py)
            good_split_names = set([split_name for split_name in collection_dict[bin_id] if split_name in split_names])
            if not len(good_split_names):
                bins_with_zero_splits_in_profile_db.append(bin_id)
                collection_dict.pop(bin_id)
                bins_info_dict.pop(bin_id)
            else:
                collection_dict[bin_id] = good_split_names

        self.progress.end()

        if len(bins_with_zero_splits_in_profile_db):
            self.run.warning('Some of the bins in this collection (precisely %d of %d total) did not contain any '
                             'that appeared in the profile database. There are multiple reasons for why this can '
                             'happen. But one of the common scenario could be this: You imported an external '
                             'collection, and some of the bins you have in that collection contain a small number '
                             'of contigs that were too short to make it into the merged profile. Well, if you would '
                             'like to figure out what might be the scenario for your experiment, here is the list of '
                             'bin names that did not go through: %s.' \
                                % (len(bins_with_zero_splits_in_profile_db), len(collection_dict), ", ".join(bins_with_zero_splits_in_profile_db)))

        return (collection_dict, bins_info_dict, split_names_in_db_but_missing_in_collection)


    def get_collection_dict(self, collection_name):
        self.sanity_check(collection_name)

        c = self.collections_dict[collection_name]

        database = db.DB(c['source_db_path'], c['source_db_version'])
        collection_dict_from_db = database.get_some_rows_from_table_as_dict(t.collections_splits_table_name, 'collection_name="%s"' % collection_name)
        database.disconnect()

        collection_dict_to_return = {}

        for entry in list(collection_dict_from_db.values()):
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
        collections_bins_info_table_filtered = get_filtered_dict(collections_bins_info_table, 'collection_name', set([collection_name]))

        bins_info_dict = {}
        for v in list(collections_bins_info_table_filtered.values()):
            bins_info_dict[v['bin_name']] = {'html_color': v['html_color'], 'source': v['source']}

        return bins_info_dict


    def is_bin_in_collection(self, collection_name, bin_name):
        self.sanity_check(collection_name)

        bins_info_dict = self.get_bins_info_dict(collection_name)

        if bin_name not in bins_info_dict:
            raise ConfigError("The bin '%s' does not seem to be a member of the collection '%s'. If you want to see all bins in "
                              "this collection you can try to add `--list-bins` to your arguments." % (bin_name, collection_name))

        return True


    def list_collections(self):
        self.run.warning('', 'COLLECTIONS FOUND', lc='yellow')
        for collection_name in self.collections_dict:
            c = self.collections_dict[collection_name]
            output = '%s (%d bins, representing %d items).' % (collection_name, c['num_bins'], c['num_splits'])
            self.run.info_single(output)


    def list_bins_in_collection(self, collection_name):
        if collection_name not in self.collections_dict:
            raise ConfigError("The collection name '%s' is not know to anyone here :/ You have to go back, Kate." % collection_name)

        self.run.warning('', 'BINS IN COLLECTION "%s"' % collection_name, lc='yellow')
        bins_info = self.get_bins_info_dict(collection_name)
        for bin_name in sorted(bins_info.keys()):
            output = '%s.' % (bin_name)
            self.run.info_single(output)


    def merge_bins(self, collection_name, new_bin_name, bin_names_list):
        """Merges a given list of bins in a collection"""

        self.sanity_check(collection_name)

        if not self.db_path:
            raise ConfigError("Something is off. The class does not know which database it is supposed to "
                              "be working with.")

        if not isinstance(bin_names_list, list):
            raise ConfigError("The `bin_names_list` must be of thpe `set` :/")

        bins_info_dict = self.get_bins_info_dict(collection_name)
        collection_dict = self.get_collection_dict(collection_name)

        invalid_bin_names = [b for b in bin_names_list if not b in collection_dict]
        if invalid_bin_names:
            raise ConfigError("Some of the bin names you want to merge is not in the collection %s :/ Here "
                              "is a list of them: %s" % (collection_name, ', '.join(invalid_bin_names)))

        items_in_new_bin = []
        for bin_name in bin_names_list:
            items_in_new_bin.extend(collection_dict[bin_name])

        info_for_new_bin = copy.deepcopy(bins_info_dict[bin_name])
        info_for_new_bin['source'] = 'anvi-merge-bins'

        # time to remove the ones that are merged
        for bin_name in bin_names_list:
            bins_info_dict.pop(bin_name)
            collection_dict.pop(bin_name)

        # add the merged stuff
        bins_info_dict[new_bin_name] = info_for_new_bin
        collection_dict[new_bin_name] = items_in_new_bin

        tables_for_collections = TablesForCollections(self.db_path, run=terminal.Run(verbose=False))
        tables_for_collections.append(collection_name, collection_dict, bins_info_dict)

        self.run.info_single("You did it. Your bins are now merged.. Onward!", nl_before=1, nl_after=1)


    def export_collection(self, collection_name, output_file_prefix=None, include_unbinned=False):
        self.sanity_check(collection_name)

        if not output_file_prefix:
            output_file_prefix = 'collection-%s' % (collection_name.strip().replace(' ', '-'))

        info_file_path = output_file_prefix + '-info.txt'
        items_file_path = output_file_prefix + '.txt'

        self.run.info('Report unbinned items if there are any', include_unbinned)
        self.run.info('Items file path', items_file_path)
        filesnpaths.is_output_file_writable(items_file_path)

        bins_info = self.get_bins_info_dict(collection_name)
        collection = self.get_collection_dict(collection_name)

        if len(bins_info):
            self.run.info('Bins info file path', info_file_path)
            info_file = open(info_file_path, 'w')

            if include_unbinned:
                bins_info['UNBINNED_ITEMS_BIN'] = {'html_color': '#000000', 'source': 'anvi-export-collections'}

            for bin_name in bins_info:
                info_file.write('%s\t%s\t%s\n' % (bin_name, bins_info[bin_name]['source'], bins_info[bin_name]['html_color']))
            info_file.close()

        binned_items = set([])

        items_file = open(items_file_path, 'w')
        for bin_name in collection:
            for item_name in collection[bin_name]:
                items_file.write('%s\t%s\n' % (item_name, bin_name))
                binned_items.add(item_name)

        if include_unbinned:
            all_items = get_all_item_names_from_the_database(self.db_path)

            unbinned_items = all_items.difference(binned_items)

            for item_name in unbinned_items:
                items_file.write('%s\tUNBINNED_ITEMS_BIN\n' % (item_name))

            self.run.warning("As per your request, %d items that were not in any of the bins in the collection '%s' are stored "
                             "in the output file under the bin name 'UNBINNED_ITEMS_BIN'." % (len(unbinned_items), collection_name))

        items_file.close()


class GetSplitNamesInBins:
    def __init__(self, args):
        # we will fill this in and return it
        self.split_names_of_interest = set([])
        self.bins = None

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bin_ids_file_path = A('bin_ids_file')
        self.bin_ids_list = A('bin_ids_list')
        self.bin_id = A('bin_id')
        self.collection_name = A('collection_name')
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.debug = anvio.DEBUG

        if not self.profile_db_path:
            raise ConfigError("You didn't provide a profile database path. When you clearly should have :/ "
                               "This is GetSplitNamesInBins speaking. Has her eyes on you.")

        if self.bin_ids_file_path and self.bin_id:
            raise ConfigError('Either use a file to list all the bin ids (-B), or declare a single bin (-b) '
                               'you would like to focus. Not both :/')

        if not self.collection_name:
            raise ConfigError('This will not work without a collection ID for your bins :/')

        if self.bin_ids_file_path:
            filesnpaths.is_file_exists(self.bin_ids_file_path)
            self.bins = set([b.strip() for b in open(self.bin_ids_file_path).readlines()])
        elif self.bin_id:
            self.bins = set([self.bin_id])

        self.collections = Collections()
        self.collections.populate_collections_dict(self.profile_db_path)

        if self.collection_name not in self.collections.collections_dict:
            progress.reset()
            raise ConfigError('The collection id "%s" does not seem to be in the profile database. These are the '
                               'collections that are available through this profile database: "%s".'\
                                                    % (self.collection_name, ', '.join(self.collections.collections_dict)))

        self.collection_dict = self.collections.get_collection_dict(self.collection_name)

        bins_in_collection = list(self.collection_dict.keys())

        if not self.bins:
            self.bins = bins_in_collection
        else:
            bins_that_do_not_exist_in_collection = [b for b in self.bins if b not in bins_in_collection]
            if len(bins_that_do_not_exist_in_collection):
                progress.reset()
                some_bins_that_exist_in_collection = bins_in_collection if len(bins_in_collection) < 30 else bins_in_collection[:30]
                raise ConfigError('Some of the bins you requested do not appear to have been described in the collection '
                                   '"%s". Here is a list of bins that are missing: "%s". Here is a list of some bins in '
                                   'your collection: "%s"' % (self.collection_name,
                                                              ', '.join(bins_that_do_not_exist_in_collection),
                                                              ', '.join(some_bins_that_exist_in_collection)))

        if not len(self.bins):
            raise ConfigError('There is no bin to work with :/')


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
