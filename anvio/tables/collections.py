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
P = terminal.pluralize


class TablesForCollections(Table):
    """Populates the collections_* tables, where collections of bins of contigs and items are kept"""
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.version = utils.get_required_version_for_db(db_path)
        self.run = run

        Table.__init__(self, self.db_path, self.version, run, progress)


    def delete(self, collection_name):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        # remove any pre-existing information for 'collection_name'
        self.delete_entries_for_key('collection_name', collection_name, [t.collections_info_table_name, t.collections_contigs_table_name, t.collections_splits_table_name, t.collections_bins_info_table_name])


    def delete_bin(self, collection_name, bin_name):
        database = db.DB(self.db_path, self.version)
        tables_to_clear = [t.collections_contigs_table_name, t.collections_splits_table_name, t.collections_bins_info_table_name]

        for table_name in tables_to_clear:
            database._exec('''DELETE FROM %s WHERE collection_name = "%s" AND \
                                                   bin_name = "%s"''' % (table_name, collection_name, bin_name))

        self.run.warning('All previous entries for "%s" of "%s" is being removed from "%s"'\
                    % (bin_name,collection_name, ', '.join(tables_to_clear)))

        database.disconnect()


    def refresh_collections_info_table(self, collection_name):
        """For a given collection, re-read most up-to-date information from the collection items table and update collections info table"""
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        collection_names_in_db = database.get_single_column_from_table(t.collections_splits_table_name, 'collection_name', unique=True)

        if collection_name not in collection_names_in_db:
            database.disconnect()
            raise ConfigError(f"The collection name '{collection_name}' is not in the collections table :/")

        where_clause = f'collection_name="{collection_name}"'
        # please note that this is not unique yet and it is intentional
        bin_names_in_collection = database.get_single_column_from_table(t.collections_splits_table_name, 'bin_name', where_clause=where_clause)

        num_splits_in_collection = len(bin_names_in_collection)
        bin_names_in_collection = sorted(list(set(bin_names_in_collection)))
        database.disconnect()

        self.delete_entries_for_key('collection_name', collection_name, [t.collections_info_table_name])

        db_entries = tuple([collection_name, num_splits_in_collection, len(bin_names_in_collection), ','.join(bin_names_in_collection)])

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_info_table_name, db_entries)
        database.disconnect()


    def append(self, collection_name, collection_dict, bins_info_dict={}, drop_collection=True):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        for bin_name in collection_dict:
            utils.is_this_name_OK_for_database('bin name', bin_name, stringent=False)

        if bins_info_dict:
            if set(collection_dict.keys()) - set(bins_info_dict.keys()):
                raise ConfigError(f"Bins in the collection dict do not match to the ones in the bins info dict. "
                                  f"They do not have to be identical, but for each bin id, there must be a unique "
                                  f"entry in the bins informaiton dict. There is something wrong with your input :/")

        if drop_collection:
            # remove any pre-existing information for 'collection_name'
            self.delete(collection_name)

        num_splits_in_collection_dict = sum([len(splits) for splits in list(collection_dict.values())])
        splits_in_collection_dict = set(list(chain.from_iterable(list(collection_dict.values()))))
        if len(splits_in_collection_dict) != num_splits_in_collection_dict:
            raise ConfigError(f"TablesForCollections::append: {(num_splits_in_collection_dict - len(splits_in_collection_dict))} "
                              f"split names or contig IDs appear more than once in your input for this collection. This part of "
                              f"the code is unable to predict how you may have ended up here, but check your input file maybe? :/")

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # how many clusters are defined in 'collection_dict'?
        bin_names = list(collection_dict.keys())

        if drop_collection:
            db_entries = tuple([collection_name, num_splits_in_collection_dict, len(bin_names), ','.join(bin_names)])
            database._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_info_table_name, db_entries)

        if not bins_info_dict:
            colors = utils.get_random_colors_dict(bin_names)
            for bin_name in bin_names:
                bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': 'UNKNOWN'}

        # populate bins info table.
        db_entries = [(collection_name, b, bins_info_dict[b]['source'], bins_info_dict[b]['html_color']) for b in bin_names]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_bins_info_table_name, db_entries)

        # populate splits table
        db_entries = []
        for bin_name in collection_dict:
            for split_name in collection_dict[bin_name]:
                db_entries.append(tuple([collection_name, split_name, bin_name]))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?)''' % t.collections_splits_table_name, db_entries)
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
                self.run.warning(f"{len(splits_only_in_collection_dict)} of {len(splits_in_collection_dict)} items found in "
                                 f"collection '{collection_name}' are not known to the contigs database. This may be OK, but "
                                 f"you must be the judge of it. If this surprises you, please use caution and make sure all "
                                 f"is fine before going forward with you analysis.")

            if len(splits_only_in_db):
                self.run.warning('%d of %d items found in the database were missing from the "%s" results. If this '
                                         'does not make any sense, please make sure you know why before going any further.'\
                                                % (len(splits_only_in_db), len(self.splits_info), collection_name))

            # then populate contigs table.
            db_entries = self.process_contigs(collection_name, collection_dict)
            database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_contigs_table_name, db_entries)

        database.disconnect()

        num_bins = len(bin_names)
        num_bins_to_report = 50
        if not drop_collection:
            bins_to_report = bin_names
            bin_report_msg = f"Here is a full list of the bin names added to this collection: {', '.join(bins_to_report)}."
        elif num_bins <= num_bins_to_report:
            bins_to_report = bin_names
            bin_report_msg = f"Here is a full list of the bin names in this collection: {', '.join(bins_to_report)}."
        else:
            bins_to_report = bin_names[:num_bins_to_report]
            bin_report_msg = f"Here is a list of the first {P('bin name', num_bins_to_report)} in this collection: {', '.join(bins_to_report)}."

        if drop_collection:
            self.run.info('Collections', f"The collection '{collection_name}' that describes {P('item', num_splits)} in {P('bin', num_bins)} was successfully "
                                         f"added to the to the database at '{self.db_path}'. {bin_report_msg}", mc='green')
        else:
            self.run.info('Collections', f"The existing collection '{collection_name}' updated and {P('split', num_splits)} in {P('bin', num_bins)} were successfully "
                                         f"added to the to the database at '{self.db_path}'. {bin_report_msg}", mc='green')


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

            db_entries_for_contigs.append(tuple([collection_name, contig_name, split_to_bin_name[split_name]]))

        return db_entries_for_contigs


    def add_default_collection_to_db(self, contigs_db_path=None, collection_name="DEFAULT", bin_name="EVERYTHING", bin_each_item_separately=False):
        """A helper function to add a default collection.

        This function will either add a collection to a given database that describes all items
        in it in a single bin, or each item in it in a separate bin.
        """

        if bin_each_item_separately:
            run.warning("Since you passed the flag '--bin-each-item-separately', anvi'o will create a separate bin for each "
                        "item in your database. This is likely a very bad idea, but anvi'o trusts that you know what you are "
                        "doing.")

        utils.is_pan_or_profile_db(self.db_path)

        if utils.get_db_type(self.db_path) == 'profile' and utils.is_blank_profile(self.db_path):
            if not contigs_db_path:
                raise ConfigError("OK, so anvi'o actually can't add a 'default' collection into a 'blank profile database' such as "
                                  "yours). The main reason behind this is that a blank profile databases do not know about "
                                  "the items names they pretend to know. Why? Well, it is a long story, BUT FORTUNATELY, you can "
                                  "recover from this by providing the contigs database your blank profile is associated with, in "
                                  "which case things will likely work.")
            else:
                utils.is_profile_db_and_contigs_db_compatible(self.db_path, contigs_db_path)
                contigs_db = db.DB(contigs_db_path, t.contigs_db_version)
                all_items = contigs_db.get_single_column_from_table(t.splits_info_table_name, 'split')
                contigs_db.disconnect()
        else:
            if contigs_db_path:
                run.warning("You should provide a contigs database to this script only if you are working with an anvi'o blank "
                            "profile database, which doesn't seem to be the case here. Anvi'o did think about raising an error "
                            "and killing your awesome analysis on its track in the name of being explicit, but then it decided "
                            "to ignore it for this once. Basically, the contigs database you provided will not be utilized for "
                            "anything.")

            all_items = utils.get_all_item_names_from_the_database(self.db_path)

        bins = {}
        if bin_each_item_separately:
            counter = 1
            for item in all_items:
                bins[f"BIN_{counter:07}"] = {item}
                counter += 1
        else:
            bins = {bin_name: all_items}

        self.append(collection_name, bins)
