# -*- coding: utf-8
# pylint: disable=line-too-long

"""TablesForCollections"""

from itertools import chain

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.tables.tableops import Table


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TablesForCollections(Table):
    """Populates the collections_* tables, where collections of bins of contigs and splits are kept"""
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, utils.get_required_version_for_db(db_path), run, progress)

        # set these dudes so we have access to unique IDs:
        self.set_next_available_id(t.collections_bins_info_table_name)
        self.set_next_available_id(t.collections_contigs_table_name)
        self.set_next_available_id(t.collections_splits_table_name)


    def delete(self, collection_name):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        # remove any pre-existing information for 'collection_name'
        self.delete_entries_for_key('collection_name', collection_name, [t.collections_info_table_name, t.collections_contigs_table_name, t.collections_splits_table_name, t.collections_bins_info_table_name])


    def append(self, collection_name, collection_dict, bins_info_dict={}):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        for bin_name in collection_dict:
            utils.is_this_name_OK_for_database('bin name', bin_name, stringent=False)

        if bins_info_dict:
            if set(collection_dict.keys()) - set(bins_info_dict.keys()):
                raise ConfigError('Bins in the collection dict do not match to the ones in the bins info dict.\
                                    They do not have to be identical, but for each bin id, there must be a unique\
                                    entry in the bins informaiton dict. There is something wrong with your input :/')

        # remove any pre-existing information for 'collection_name'
        self.delete(collection_name)

        num_splits_in_collection_dict = sum([len(splits) for splits in list(collection_dict.values())])
        splits_in_collection_dict = set(list(chain.from_iterable(list(collection_dict.values()))))
        if len(splits_in_collection_dict) != num_splits_in_collection_dict:
            raise ConfigError("TablesForCollections::append: %d of the split or contig IDs appear more than once in\
                                your collections input. It is unclear to anvi'o how did you manage to do this, but we\
                                cannot go anywhere with this :/" % (num_splits_in_collection_dict - len(splits_in_collection_dict)))

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # how many clusters are defined in 'collection_dict'?
        bin_names = list(collection_dict.keys())

        # push information about this search result into serach_info table.
        db_entries = tuple([collection_name, num_splits_in_collection_dict, len(bin_names), ','.join(bin_names)])
        database._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_info_table_name, db_entries)

        if not bins_info_dict:
            colors = utils.get_random_colors_dict(bin_names)
            for bin_name in bin_names:
                bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': 'UNKNOWN'}

        # populate bins info table.
        db_entries = [(self.next_id(t.collections_bins_info_table_name), collection_name, b, bins_info_dict[b]['source'], bins_info_dict[b]['html_color']) for b in bin_names]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.collections_bins_info_table_name, db_entries)

        # populate splits table
        db_entries = []
        for bin_name in collection_dict:
            for split_name in collection_dict[bin_name]:
                db_entries.append(tuple([self.next_id(t.collections_splits_table_name), collection_name, split_name, bin_name]))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_splits_table_name, db_entries)
        num_splits = len(db_entries)


        # FIXME: This function can be called to populate the contigs database (via anvi-populate-collections), or
        # the profile database. when it is contigs database, the superclass Table has the self.splits_info variable
        # set when it is initialized. however, the Table instance is missing self.splis when it is initialized with
        # the profile database. hence some special controls for contigs db (note that collections_contigs_table is
        # only populated in the contigs database):
        if self.db_type == 'contigs':
            splits_only_in_collection_dict = [c for c in splits_in_collection_dict if c not in self.splits_info]
            splits_only_in_db = [c for c in self.splits_info if c not in splits_in_collection_dict]

            if len(splits_only_in_collection_dict):
                self.run.warning('%d of %d splits found in "%s" results are not in the database. This may be OK,\
                                          but you must be the judge of it. If this is somewhat surprising, please use caution\
                                          and make sure all is fine before going forward with you analysis.'\
                                                % (len(splits_only_in_collection_dict), len(splits_in_collection_dict), collection_name))

            if len(splits_only_in_db):
                self.run.warning('%d of %d splits found in the database were missing from the "%s" results. If this\
                                          does not make any sense, please make sure you know why before going any further.'\
                                                % (len(splits_only_in_db), len(self.splits_info), collection_name))

            # then populate contigs table.
            db_entries = self.process_contigs(collection_name, collection_dict)
            database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_contigs_table_name, db_entries)

        database.disconnect()

        num_bins = len(bin_names)
        num_bins_to_report = 50
        if num_bins <= num_bins_to_report:
            bins_to_report = bin_names
            bin_report_msg = "Here is a full list of the bin names in this collection: {}.".format(",".join(bins_to_report))
        else:
            bins_to_report = bin_names[:num_bins_to_report]
            bin_report_msg = "Here is a list of the first {} bin names in this collection: {}.".format(num_bins_to_report, ",".join(bins_to_report))

        self.run.info('Collections', 'The collection "%s" that describes %s splits and %s bins has been successfully added to the\
                                      database at "%s". %s' % (collection_name, pp(num_splits), pp(num_bins), self.db_path, bin_report_msg), mc='green')


    def process_contigs(self, collection_name, collection_dict):
        db_entries_for_contigs = []

        split_to_bin_name = {}
        for bin_name in collection_dict:
            for split_name in collection_dict[bin_name]:
                split_to_bin_name[split_name] = bin_name

        contigs_processed = set([])
        for split_name in split_to_bin_name:
            if split_name not in self.splits_info:
                # which means this split only appears in the input file, but not in the database.
                continue

            contig_name = self.splits_info[split_name]['parent']

            if contig_name in contigs_processed:
                continue
            else:
                contigs_processed.add(contig_name)

            db_entry = tuple([self.next_id(t.collections_contigs_table_name), collection_name, contig_name, split_to_bin_name[split_name]])
            db_entries_for_contigs.append(db_entry)

        return db_entries_for_contigs
