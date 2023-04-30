# -*- coding: utf-8
# pylint: disable=line-too-long
"""The library to split merged profiles into smaller profiles.

The default client of this library is under bin/anvi-split"""


import os
import sys
import copy
import argparse

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.profiler as profiler
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.panops import Pangenome
from anvio.clusteringconfuguration import ClusteringConfiguration
from anvio.tables.kmers import KMerTablesForContigsAndSplits
from anvio.tables.collections import TablesForCollections
from anvio.tables.genefunctions import TableForGeneFunctions
from anvio.tables.views import TablesForViews


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
run = terminal.Run()
progress = terminal.Progress()


class PanSplitter(summarizer.PanSummarizer):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.pan_db_path = A('pan_or_profile_db')
        args.pan_db = self.pan_db_path

        self.genomes_storage_path = A('genomes_storage')
        self.collection_name = A('collection_name')
        self.bin_name = A('bin_id')
        self.output_directory = A('output_dir')
        self.list_collections = A('list_collections')

        self.collections = ccollections.Collections()
        self.collections.populate_collections_dict(self.pan_db_path)


    def init(self):
        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=True)

        if not self.genomes_storage_path:
            raise ConfigError("You must provide a genomes storage database for this operation.")

        if not self.pan_db_path:
            raise ConfigError("You came all the way here without a pan database. Congratulations! But we "
                              "kinda need it at this stage really :/ GOOD DAY.")

        utils.is_pan_db(self.pan_db_path)

        if self.list_collections:
            self.collections.list_collections()
            sys.exit(0)

        if not self.collection_name:
            raise ConfigError("You must provide a collection name for this to work. If you want to know about "
                              "all the collections in your pan database you can use the program "
                              "`anvi-show-collections-and-bins` or run the same command with the flag "
                              "`--list-collections`.")

        # if this is not set false, the summarizer class attemts to remove the main output directory
        # upon initialization. not doing that is useful in this context since this allows multiple
        # anvi-split runs to work on bins in the same collection in parallel:
        self.args.delete_output_directory_if_exists = False

        self.summary = summarizer.PanSummarizer(self.args, r=self.run, p=self.progress)
        self.summary.load_pan_views()

        self.bin_names_of_interest = sorted(self.summary.bins_info_dict.keys())

        if self.bin_name:
            if self.bin_name not in self.bin_names_of_interest:
                raise ConfigError("The bin name you wish to split from this profile databse is not in the collection. Busted!")
            else:
                self.bin_names_of_interest = [self.bin_name]


    def process(self):
        """This is the function that goes through each bin loaded in the class and proecesses them."""
        self.init()

        filesnpaths.gen_output_directory(self.output_directory)

        self.run.warning("Anvi'o is about to start splitting your bins into individual, self-contained anvi'o profiles. This "
                         "is quite a tricky operation, and even if it finishes successfully, you must double check everyting "
                         "in the resulting profiles to make sure things worked as expected. Although we are doing our best to "
                         "test all these, variation between projects make it impossible to be 100% sure.")

        for bin_name in self.bin_names_of_interest:
            b = PanBinSplitter(bin_name, self.summary, self.args, run=self.run, progress=self.progress)
            b.do_pan_db()

        self.run.info('Num bins processed', len(self.bin_names_of_interest))
        self.run.info("Output directory", self.output_directory)


class ProfileSplitter:
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('pan_or_profile_db')
        args.profile_db = self.profile_db_path

        self.contigs_db_path = A('contigs_db')
        self.collection_name = A('collection_name')
        self.bin_name = A('bin_id')
        self.output_directory = A('output_dir')
        self.skip_variability_tables = A('skip_variability_tables')

        self.collections = ccollections.Collections()
        self.summary = None


    def sanity_check(self):
        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=True)

        if not self.contigs_db_path:
            raise ConfigError("You must provide a contigs database for this operation.")

        if not self.profile_db_path:
            raise ConfigError("No profile db no cookie. Bye.")

        utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        if profile_db.meta['db_type'] != 'profile':
            raise ConfigError("Anvi'o was trying to split this profile, but it just realized that it is not a profile "
                              "database. There is something wrong here.")
        profile_db.disconnect()

        # if this is not set false, the summarizer class attemts to remove the main output directory
        # upon initialization. not doing that is useful in this context since this allows multiple
        # anvi-split runs to work on bins in the same collection in parallel:
        self.args.delete_output_directory_if_exists = False

        self.summary = summarizer.ProfileSummarizer(self.args, r=self.run, p=self.progress)
        self.summary.init()

        self.bin_names_of_interest = sorted(self.summary.bin_ids)
        if self.bin_name:
            if self.bin_name not in self.bin_names_of_interest:
                raise ConfigError("The bin name you wish to split from this profile database is not in the collection. Busted!")
            else:
                self.bin_names_of_interest = [self.bin_name]


    def process(self):
        """This is the function that goes through each bin loaded in the class and proecesses them."""
        self.sanity_check()

        filesnpaths.gen_output_directory(self.output_directory)

        self.run.warning("Anvi'o is about to start splitting your bins into individual, self-contained anvi'o profiles. As of "
                         "2021, we have tested this feature quite extensively and we trust that it will do well. But this is "
                         "still quite a tricky operation and you must double-check things once your split data is ready.",
                         header="ANVI'O TRICKY OPERATIONS DEPARTMENT", lc="green")

        if self.skip_variability_tables:
            self.run.warning("Since you asked so nicely, anvi'o will not migrate variability table data into split profiles.")

        if self.summary.p_meta['blank']:
            self.run.warning("It seems your profile database is a blank one. That's fine. Anvi'o assumes that your actual "
                             "intention is to split your contigs database only. This warning message is here to make sure "
                             "you will not be upset when you realize your split profile missing a profile database :(", lc="yellow")

        for bin_name in self.bin_names_of_interest:
            b = BinSplitter(bin_name, self.summary, self.args, run=self.run, progress=self.progress)
            b.do_contigs_db()

            if not self.summary.p_meta['blank']:
                b.do_profile_db()

                if self.summary.auxiliary_profile_data_available:
                    b.do_auxiliary_profile_data()

        self.run.info('Num bins processed', len(self.bin_names_of_interest))
        self.run.info("Output directory", self.output_directory)


class XSplitter(object):
    """Some common functions both for profile bin and pan bin splitter classes"""

    def __init__(self):
        pass

    def is_hierarchical_clustering_for_bin_OK(self):
        skip_hierarchical_clustering = self.skip_hierarchical_clustering

        if self.num_splits > self.max_num_splits_for_hierarchical_clustering and not self.enforce_hierarchical_clustering:
            self.run.warning("It seems you have more than %s splits in this particular bin. This is the "
                             "soft limit for anvi'o to attempt to create a hierarchical clustering of your splits "
                             "(which becomes the center tree in all anvi'o displays). If you want a hierarchical "
                             "clustering to be done anyway, you can re-run the splitting process only for this bin "
                             "by adding these parameters to your run: '--bin-id %s --enforce-hierarchical-clustering'. "
                             "If you feel like you are lost, don't hesitate to get in touch with anvi'o developers." \
                                                        % (pp(self.max_num_splits_for_hierarchical_clustering), self.bin_id))
            skip_hierarchical_clustering = True

        if self.num_splits > self.max_num_splits_for_hierarchical_clustering and self.enforce_hierarchical_clustering:
            self.run.warning("Becasue you have used the flag `--enforce-hierarchical-clustering`, anvi'o will attempt "
                             "to create a hierarchical clustering of your %s splits for this bin. It may take a bit of "
                             "time, and it is not even anvi'o's fault, you know  :/" \
                                                        % pp(self.max_num_splits_for_hierarchical_clustering))

        return skip_hierarchical_clustering


    def is_dbs_identical(self, source_db_path, target_db_path):
        """Check whether the two dbs have identical table names"""

        source_db = db.DB(source_db_path, None, ignore_version=True, skip_rowid_prepend=True)
        target_db = db.DB(target_db_path, None, ignore_version=True, skip_rowid_prepend=True)

        source_tables = set(source_db.get_table_names())
        target_tables = set(target_db.get_table_names())

        if not source_tables == target_tables:
            self.progress.end()

            diff = source_tables.difference(target_tables)

            raise ConfigError("Something went wrong during subsetting :/ Table names in the parent db (%s) and the child "
                              "db (%s) does not seem to be identical. The following tables are found in the source, but "
                              "missing in the target database: '%s'" % (source_db_path, target_db_path, ', '.join(diff)))

        source_db.disconnect()
        target_db.disconnect()


    def migrate_data(self, tables_dict, source_db_path, target_db_path):
        """Filter data from `source_db_path` into `target_db_path` based on rules defined in `tables_dict`

        Each item in `tables_dict` must be a table name, and should correspond to a tuple with two items. The
        first item of the tuple should be a `key` (of type str) on which the table data should be filtered, and
        the second item should be a set of acceptable values (of type set).
        """

        for table_name in tables_dict:
            self.progress.update("Table '%s' .. setting up the query" % table_name)
            filter_on_key = tables_dict[table_name][0]
            filter_for = ','.join('"%s"' % str(i) for i in tables_dict[table_name][1])
            where_clause = "%s IN (%s)" % (filter_on_key, filter_for)


            self.progress.update("Table '%s' .. reading data" % table_name)
            source_db = db.DB(source_db_path, None, ignore_version=True, skip_rowid_prepend=True)
            data = source_db.get_some_rows_from_table(table_name, where_clause)
            source_db.disconnect()

            if not len(data):
                continue

            self.progress.update("Table '%s' .. writing data" % table_name)
            target_db = db.DB(target_db_path, None, ignore_version=True, skip_rowid_prepend=True)
            target_db._exec_many('''INSERT INTO %s VALUES(%s)''' % (table_name, ','.join(['?'] * len(data[0]))), data)
            target_db.disconnect()

        # make sure things are OK
        self.is_dbs_identical(source_db_path, target_db_path)



class PanBinSplitter(summarizer.PanBin, XSplitter):
    def __init__(self, bin_name, summary_object, args, run=run, progress=progress):
        """A class to split a single bin from its parent.

        The class is not really useful without a summary object, but it makes logistic sense to keep it
        separate since the inheritance from anvio/summarizer.Bin is much easier and sane this way."""
        summarizer.PanBin.__init__(self, summary_object, bin_name, run, progress)

        XSplitter.__init__(self)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.pan_db_path = A('pan_db')
        self.genomes_storage_path = A('genomes_storage')
        self.output_directory = A('output_dir')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default

        # let's remember this for later.
        self.bin_project_name = 'The %s split from "%s"' % (self.bin_id, self.summary.p_meta['project_name'])

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)
        self.clustering_configs = constants.clustering_configs['pan']
        self.database_paths = {'PAN.db': os.path.abspath(self.pan_db_path)}


        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering "
                              "while also asking it to enforce it.")

        # set the output directory, and output file paths
        self.bin_output_directory = os.path.join(self.output_directory, bin_name)
        filesnpaths.gen_output_directory(self.bin_output_directory)

        # let's see whether we are going to do any hierarchical clustering:
        self.max_num_splits_for_hierarchical_clustering = constants.max_num_items_for_hierarchical_clustering
        self.skip_hierarchical_clustering = self.is_hierarchical_clustering_for_bin_OK()

        # set your own db paths
        self.bin_pan_db_path = os.path.join(self.bin_output_directory, 'PAN.db')


    def do_pan_db(self):
        self.progress.new('Splitting "%s"' % self.bin_id)
        self.progress.update('Subsetting the pan database')

        bin_pan_db = dbops.PanDatabase(self.bin_pan_db_path)
        bin_pan_db.touch()

        # copy-paste tables that will largely stay the same from the parent
        bin_pan_db.db.copy_paste(table_name='self', source_db_path=self.pan_db_path)
        bin_pan_db.db.copy_paste(table_name=t.states_table_name, source_db_path=self.pan_db_path)
        bin_pan_db.db.copy_paste(table_name=t.layer_additional_data_table_name, source_db_path=self.pan_db_path)
        bin_pan_db.db.copy_paste(table_name=t.layer_orders_table_name, source_db_path=self.pan_db_path)

        # update some values
        bin_pan_db.db.update_meta_value('project_name', self.bin_project_name)
        bin_pan_db.db.update_meta_value('available_item_orders', None)
        bin_pan_db.db.update_meta_value('items_ordered', None)
        bin_pan_db.db.update_meta_value('num_gene_clusters', self.num_gene_clusters)
        bin_pan_db.db.update_meta_value('num_genes_in_gene_clusters', self.num_genes_in_gene_clusters)

        bin_pan_db.disconnect()

        # summarizer.PanBin already has updated/pruned views dicts. that's why this loop will work.
        for view_name in ['gene_cluster_frequencies', 'gene_cluster_presence_absence']:
            TablesForViews(self.bin_pan_db_path).create_new_view(
                                                view_data=self.views[view_name]['dict'],
                                                table_name=self.views[view_name]['table_name'],
                                                view_name = view_name,
                                                from_matrix_form=True)

        # setup the filtering rules for migrating data:
        tables = {
                    t.item_additional_data_table_name: ('item_name', self.split_names),
                    t.pan_gene_clusters_table_name: ('gene_cluster_id', self.split_names),
                }

        self.migrate_data(tables, self.pan_db_path, self.bin_pan_db_path)

        self.progress.end()

        # add a collection
        collection_dict = {'ALL_SPLITS': self.split_names}
        bins_info_dict = {'ALL_SPLITS': {'html_color': '#FF0000', 'source': 'anvi-split'}}
        collections = TablesForCollections(self.bin_pan_db_path)
        collections.append('DEFAULT', collection_dict, bins_info_dict=bins_info_dict)

        # clustering of items.. this is the most elegant way of doing this:
        p = Pangenome(argparse.Namespace(skip_hierarchical_clustering=self.skip_hierarchical_clustering,
                                         output_dir=self.bin_output_directory,
                                         distance=self.distance,
                                         linkage=self.linkage,
                                         run=self.run,
                                         progress=self.progress,
                                         project_name=self.bin_project_name))
        p.genomes = self.genomes
        p.pan_db_path = self.bin_pan_db_path
        p.gen_hierarchical_clustering_of_gene_clusters()


class DBSplitter:
    """Implements the shittiest factory pattern humankind has ever seen. Behold."""

    def __init__(self, args, run=run, progress=progress):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        if not A('pan_or_profile_db'):
            raise ConfigError("No pan/profile database no cookie.")

        self.mode = None
        if A('contigs_db'):
            self.mode = 'profile'
        elif A('genomes_storage'):
            self.mode = 'pan'
        else:
            raise ConfigError("Well. You are trying to initiate the anvi'o database splitter, but anvi'o has "
                              "no idea what exactly you are trying to do becasue you haven't declared enough "
                              "databases. You should either use a contigs database or a genomes storage among "
                              "your arguments to initiate this class properly.")


    def get(self):
        if self.mode == 'pan':
            return PanSplitter
        elif self.mode == 'profile':
            return ProfileSplitter
        else:
            return None


class BinSplitter(summarizer.Bin, XSplitter):
    def __init__(self, bin_name, summary_object, args, run=run, progress=progress):
        """A class to split a single bin from its parent.

        The class is not really useful without a summary object, but it makes logistic sense to keep it
        separate since the inheritance from anvio/summarizer.Bin is much easier and sane this way."""
        summarizer.Bin.__init__(self, summary_object, bin_name, run, progress)

        XSplitter.__init__(self)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.output_directory = A('output_dir')
        self.skip_variability_tables = A('skip_variability_tables')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.compress_auxiliary_data = A('compress_auxiliary_data')

        self.db_variant = str(utils.get_db_variant(self.profile_db_path))

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)
        self.clustering_configs = constants.clustering_configs['merged']
        self.database_paths = {'CONTIGS.db': os.path.abspath(self.contigs_db_path)}

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering "
                              "while also asking it to enforce it.")

        # set the output directory, and output file paths
        self.bin_output_directory = os.path.join(self.output_directory, bin_name)
        filesnpaths.gen_output_directory(self.bin_output_directory)

        # let's see whether we are going to do any hierarchical clustering:
        self.max_num_splits_for_hierarchical_clustering = constants.max_num_items_for_hierarchical_clustering
        self.skip_hierarchical_clustering = self.is_hierarchical_clustering_for_bin_OK()

        # set your own db paths
        self.bin_contigs_db_path = os.path.join(self.bin_output_directory, 'CONTIGS.db')
        self.bin_profile_db_path = os.path.join(self.bin_output_directory, 'PROFILE.db')


    def do_contigs_db(self):
        self.progress.new('Splitting "%s"' % self.bin_id)
        self.progress.update('Subsetting the contigs database')

        bin_contigs_db = dbops.ContigsDatabase(self.bin_contigs_db_path)
        bin_contigs_db.touch()

        # copy-paste tables that will largely stay the same from the parent
        bin_contigs_db.db.copy_paste(table_name='self', source_db_path=self.contigs_db_path)
        bin_contigs_db.db.copy_paste(table_name='hmm_hits_info', source_db_path=self.contigs_db_path)
        bin_contigs_db.db.copy_paste(table_name='taxon_names', source_db_path=self.contigs_db_path)

        # update some variables in the self table:
        self.contigs_db_hash = bin_contigs_db.get_hash()
        bin_contigs_db.db.update_meta_value('num_contigs', self.num_contigs)
        bin_contigs_db.db.update_meta_value('num_splits', self.num_splits)
        bin_contigs_db.db.update_meta_value('total_length', self.total_length)
        bin_contigs_db.db.update_meta_value('creation_date', bin_contigs_db.get_date())
        bin_contigs_db.db.update_meta_value('contigs_db_hash', self.contigs_db_hash)
        bin_contigs_db.db.update_meta_value('project_name', self.bin_id)

        # the empty contigs db is ready
        bin_contigs_db.disconnect()

        # touch does not create the k-mers tables, so the resulting contigs db is missing them. we
        # will add them to the db here.
        bin_contigs_db = dbops.ContigsDatabase(self.bin_contigs_db_path)
        k = KMerTablesForContigsAndSplits(None, k=bin_contigs_db.meta['kmer_size'])
        for table_name in ['kmer_contigs', 'kmer_splits']:
            bin_contigs_db.db.create_table(table_name, k.kmers_table_structure, k.kmers_table_types)
        bin_contigs_db.disconnect()

        # setup the filtering rules for migrating data:
        tables = {
                    t.contig_sequences_table_name: ('contig', self.contig_names),
                    t.contigs_info_table_name: ('contig', self.contig_names),
                    t.gene_function_calls_table_name: ('gene_callers_id', self.gene_caller_ids),
                    t.gene_amino_acid_sequences_table_name: ('gene_callers_id', self.gene_caller_ids),
                    t.genes_in_contigs_table_name: ('gene_callers_id', self.gene_caller_ids),
                    t.genes_in_splits_table_name: ('gene_callers_id', self.gene_caller_ids),
                    t.genes_taxonomy_table_name: ('gene_callers_id', self.gene_caller_ids),
                    t.hmm_hits_table_name: ('gene_callers_id', self.gene_caller_ids),
                    t.hmm_hits_splits_table_name: ('split', self.split_names),
                    t.splits_info_table_name: ('split', self.split_names),
                    t.splits_taxonomy_table_name: ('split', self.split_names),
                    t.nt_position_info_table_name: ('contig_name', self.contig_names),
                    t.scg_taxonomy_table_name: ('gene_callers_id', self.gene_caller_ids),
                    'kmer_contigs': ('contig', self.split_names),
                    'kmer_splits': ('contig', self.split_names),
                }

        self.migrate_data(tables, self.contigs_db_path, self.bin_contigs_db_path)

        # We're done here in theroy, but there is one more thing to do due to reasons partially explained in
        # issue https://github.com/merenlab/anvio/issues/1593 and PR https://github.com/merenlab/anvio/pull/1595.
        # The solution presented in the PR does not apply to split projects. so here we will calculate
        # what percentage of HMM hits are in splits described in this bin, and remove those that are less
        # than 100%.
        bin_contigs_db = dbops.ContigsDatabase(self.bin_contigs_db_path)
        hmm_hits_in_splits_dict = bin_contigs_db.db.get_table_as_dict(t.hmm_hits_splits_table_name)

        # the purpose of the folloing dict is to keep track of what total percentage of a given HMM hit is
        # described by all contig splits involved in this bin
        hmm_hits_id_percentage_described_dict = Counter({})
        for entry in hmm_hits_in_splits_dict.values():
            hmm_hits_id_percentage_described_dict[entry['hmm_hit_entry_id']] += entry['percentage_in_split']

        # now the `hmm_hits_id_percentage_described_dict` looks like this:
        #
        #   {2: 100, 3: 100.0, 5: 90.86727989487517, 6: 99.99999999999999, 4: 63.99858956276446}
        #
        # HMM hit ids that need to be cleared out from th `hmm_hits_in_splits` table is clear: 5 and 4, in this
        # example. But the problem is, due floating point logistics, in some cases things are not quite 100%,
        # although in reality they are, hence the need for `round`ing the percentages below.
        hmm_hit_ids_to_delete = [hit_id for hit_id in hmm_hits_id_percentage_described_dict if round(hmm_hits_id_percentage_described_dict[hit_id]) < 100]
        where_clause = f"hmm_hit_entry_id IN ({','.join([str(i) for i in hmm_hit_ids_to_delete])})"
        bin_contigs_db.db.remove_some_rows_from_table(t.hmm_hits_splits_table_name, where_clause=where_clause)
        bin_contigs_db.disconnect()

        self.progress.end()


    def do_auxiliary_profile_data(self):
        self.progress.new('Splitting "%s"' % self.bin_id)
        self.progress.update('Subsetting the auxiliary data (for profile db)')

        new_auxiliary_profile_data_path = dbops.get_auxiliary_data_path_for_profile_db(self.bin_profile_db_path)
        parent_auxiliary_profile_data_path = self.summary.auxiliary_data_path

        bin_profile_auxiliary = auxiliarydataops.AuxiliaryDataForSplitCoverages(new_auxiliary_profile_data_path,
                                                                                self.contigs_db_hash,
                                                                                db_variant=self.db_variant,
                                                                                create_new=True)

        parent_profile_auxiliary = auxiliarydataops.AuxiliaryDataForSplitCoverages(parent_auxiliary_profile_data_path,
                                                                                   self.summary.a_meta['contigs_db_hash'],
                                                                                   db_variant=self.db_variant)

        for split_name in self.split_names:
            sample_coverages = parent_profile_auxiliary.get(split_name)
            for sample_name in sample_coverages:
                bin_profile_auxiliary.append(split_name, sample_name, sample_coverages[sample_name])

        bin_profile_auxiliary.store()
        bin_profile_auxiliary.close()
        parent_profile_auxiliary.close()

        if self.compress_auxiliary_data:
            self.progress.update('Compressing the profile db auxiliary data file ...')
            utils.gzip_compress_file(new_auxiliary_profile_data_path)

        self.progress.end()


    def do_profile_db(self):
        # are we working with a merged profile database?
        merged = self.summary.p_meta['merged']
        self.run.info('Merged database', 'True' if merged else 'False')

        self.progress.new('Splitting "%s"' % self.bin_id)
        self.progress.update('Subsetting the %s profile database' % 'merged' if merged else 'single')

        bin_profile_db = dbops.ProfileDatabase(self.bin_profile_db_path)
        bin_profile_db.touch()

        # copy-paste tables that will largely stay the same from the parent
        bin_profile_db.db.copy_paste(table_name='self', source_db_path=self.profile_db_path)
        bin_profile_db.db.copy_paste(table_name='views', source_db_path=self.profile_db_path)
        bin_profile_db.db.copy_paste(table_name='states', source_db_path=self.profile_db_path)
        bin_profile_db.db.copy_paste(table_name='layer_additional_data', source_db_path=self.profile_db_path)
        bin_profile_db.db.copy_paste(table_name='layer_orders', source_db_path=self.profile_db_path)

        # update some values
        bin_profile_db.db.update_meta_value('contigs_db_hash', self.contigs_db_hash)
        bin_profile_db.db.update_meta_value('available_item_orders', None)
        bin_profile_db.db.update_meta_value('sample_id', self.bin_id)

        # setup the filtering rules for migrating data:
        tables = {}

        # dealing with 'view' data tables
        for table_name in constants.essential_data_fields_for_anvio_profiles:
            for target in ['splits', 'contigs']:
                new_table_name = '_'.join([table_name, target])
                new_table_structure = t.view_table_structure
                new_table_types = t.view_table_types
                bin_profile_db.db.create_table(new_table_name, new_table_structure, new_table_types)

                tables[new_table_name] = ('item', self.split_names)

        # we need to migrate these guys, too. unless we don't need to... if we are migrating,
        # the values in the self table are already accurate. if we are skipping, regardless
        # of what the values were, we will set the absolut correct ones.
        if self.skip_variability_tables:
            bin_profile_db.db.update_meta_value('SNVs_profiled', False)
            bin_profile_db.db.update_meta_value('SCVs_profiled', False)
            bin_profile_db.db.update_meta_value('INDELs_profiled', False)
        else:
            tables[t.variable_nts_table_name] = ('split_name', self.split_names)
            tables[t.variable_codons_table_name] = ('corresponding_gene_call', self.gene_caller_ids)
            tables[t.indels_table_name] = ('split_name', self.split_names)

        bin_profile_db.disconnect()

        self.migrate_data(tables, self.profile_db_path, self.bin_profile_db_path)

        self.progress.end()

        if not self.skip_hierarchical_clustering:
            dbops.do_hierarchical_clustering_of_items(self.bin_profile_db_path, constants.clustering_configs['merged' if merged else 'single'], self.split_names, \
                                                      self.database_paths, input_directory=self.bin_output_directory, \
                                                      default_clustering_config=constants.merged_default, distance=self.distance, \
                                                      linkage=self.linkage, run=terminal.Run(verbose=False), progress=self.progress)

        # add a collection
        collection_dict = {'ALL_SPLITS': self.split_names}
        bins_info_dict = {'ALL_SPLITS': {'html_color': '#FF0000', 'source': 'anvi-split'}}
        collections = TablesForCollections(self.bin_profile_db_path)
        collections.append('DEFAULT', collection_dict, bins_info_dict=bins_info_dict)


    def cluster_contigs_anvio(self):
        # clustering of contigs is done for each configuration file under static/clusterconfigs/merged directory;
        # at this point we don't care what those recipes really require because we already merged and generated
        # every data file that may be required.

        self.run.info_single("Anvi'o hierarchical clustering of contigs...", nl_before=1, nl_after=1, mc="blue")

        if not self.skip_hierarchical_clustering:
            for config_name in self.clustering_configs:
                config_path = self.clustering_configs[config_name]

                config = ClusteringConfiguration(config_path, self.output_directory, db_paths=self.database_paths, row_ids_of_interest=self.split_names)

                try:
                    clustering_id, newick = clustering.order_contigs_simple(config, distance=self.distance, linkage=self.linkage, progress=self.progress)
                except Exception as e:
                    self.run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
                    self.progress.end()
                    continue

                _, distance, linkage = clustering_id.split(':')

                dbops.add_items_order_to_db(anvio_db_path=self.profile_db_path,
                                            order_name=config_name,
                                            order_data=newick,
                                            distance=distance,
                                            linkage=linkage,
                                            make_default=config_name == constants.merged_default,
                                            run=self.run)


class LocusSplitter:
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        # following will be filled in during init:
        self.targets = []
        self.num_genes_list = None
        self.gene_caller_ids_of_interest = set([])

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_contigs_db_path = A('contigs_db')
        self.num_genes = A('num_genes')
        self.search_term = A('search_term')
        self.gene_caller_ids = A('gene_caller_ids')
        self.delimiter = A('delimiter')
        self.output_dir = A('output_dir') or os.path.abspath(os.path.curdir)
        self.output_file_prefix = A('output_file_prefix') or 'the_user_provided_no_prefix_so_here_we_go_prefix'
        self.use_hmm = A('use_hmm')
        self.hmm_sources = A('hmm_sources') or set([])
        self.annotation_sources = A('annotation_sources')
        self.remove_partial_hits = A('remove_partial_hits')
        self.reverse_complement_if_necessary = not A('never_reverse_complement')
        self.include_fasta_output = True
        self.is_in_flank_mode = bool(A('flank_mode'))

        if self.annotation_sources:
            self.annotation_sources = self.annotation_sources.split(self.delimiter)

        if A('list_hmm_sources'):
            dbops.ContigsDatabase(self.input_contigs_db_path).list_available_hmm_sources()
            sys.exit()

        # unless we are in debug mode, let's keep things quiet.
        if anvio.DEBUG:
            self.run_object = terminal.Run()
        else:
            self.run_object = terminal.Run(verbose=False)


    def sanity_check(self):
        """Check sanity while straightening some input variables"""

        filesnpaths.is_output_dir_writable(self.output_dir)

        if (not (self.gene_caller_ids or self.search_term)) or (self.gene_caller_ids and self.search_term):
            raise ConfigError("You must specify exactly one of the following: --gene-caller-ids or --search-term")

        if self.use_hmm and not self.search_term:
            raise ConfigError("If you want to use HMMs to find the gene of interest that will define your locus,\
                               you must also specify a --search-term.")

        if self.is_in_flank_mode and self.use_hmm:
            raise ConfigError("Anvi'o currently cannot use hmm search terms in flank-mode. If this "
               "functionality is needed for your analysis, please make a issue on the github "
               "repository page and we will address it.")

        if self.gene_caller_ids and self.is_in_flank_mode:
            num_genes = len(utils.get_gene_caller_ids_from_args(self.gene_caller_ids, delimiter=self.delimiter))
            if num_genes != 2:
                raise ConfigError("You're in flank mode and opted to use gene caller ids to identify the "
                                  "flanking genes. But you provided anvi'o %d gene caller id, and anvi'o "
                                  "needs exactly 2." % num_genes)

        if self.search_term:
            self.search_term = self.search_term.split(self.delimiter)

        utils.is_contigs_db(self.input_contigs_db_path)

        if len(self.hmm_sources):
            self.hmm_sources = set([s.strip() for s in self.hmm_sources.split(self.delimiter)])

        # If user is in default mode, they MUST provide --num-genes
        if not self.is_in_flank_mode:
            if not self.num_genes:
                raise ConfigError("You must provide --num-genes when in default mode.")

        if self.num_genes:
            self.num_genes_list = [int(x) for x in self.num_genes.split(self.delimiter)]
            if len(self.num_genes_list) > 2:
                raise ConfigError("The block size you provided, \"%s\", is not valid. "
                                   "The gene block size is defined by only one or two integers for either "
                                   "a block following the search match or a block preceding and following "
                                   "the search match respectively (e.g., 3,2)." % self.num_genes)

            if len(self.num_genes_list) == 1:
                self.num_genes_list = [0, self.num_genes_list[0]]

            if self.delimiter in self.num_genes:
                self.run.info('Genes to report', '%d genes before the matching gene, and %d that follow' % (self.num_genes_list[0], self.num_genes_list[1]))
            else:
                self.run.info('Genes to report', 'Matching gene, and %d genes after it' % (self.num_genes_list[0]))
        
        utils.check_sample_id(self.output_file_prefix)

        self.run.warning(None, header="Input / Output", lc="cyan")
        self.run.info('Contigs DB', os.path.abspath(self.input_contigs_db_path))
        self.run.info('Output directory', self.output_dir)
        self.run.info('Rev-comp the locus sequence if necessary', self.reverse_complement_if_necessary)


    def init(self):
        """Init calls sanity check and parses all input to get self.gene_caller_ids_of_interest"""

        self.sanity_check()

        self.search_term_to_gene_id_hits_dict = {}

        self.run.warning(None, header="Initialization bleep bloops", lc="cyan")

        if self.gene_caller_ids:
            self.run.info('Mode', 'User-provided gene caller id(s)')

            gene_caller_ids_of_interest = list(utils.get_gene_caller_ids_from_args(self.gene_caller_ids, self.delimiter))
            self.sources = ['gene_caller_ids']
        elif self.use_hmm:
            self.run.info('Mode', 'HMM search')

            s = hmmops.SequencesForHMMHits(self.input_contigs_db_path, sources=self.hmm_sources)

            self.run.info('Search term', self.search_term, mc='green')
            self.run.info('HMM sources being used', ', '.join(s.sources))

            hmm_hits = utils.get_filtered_dict(s.hmm_hits, 'gene_name', set(self.search_term))
            gene_caller_ids_of_interest = [entry['gene_callers_id'] for entry in hmm_hits.values()]

            self.targets.append('HMMs')
            self.sources = s.sources
        else:
            self.run.info('Mode', 'Function search')

            contigs_db = dbops.ContigsSuperclass(self.args, r=self.run_object)
            contigs_db.init_functions()

            # use functional annotation
            # if self.search_term:
            gene_caller_ids_of_interest = []
            counter = 1
            for term in self.search_term:
                self.run.info('Search term %d of %d' % (counter,len(self.search_term)), term, mc='green')
                self.run.info('Function calls being used', ', '.join((contigs_db.gene_function_call_sources
                                                                      if not self.annotation_sources
                                                                      else self.annotation_sources)))

                foo, search_report = contigs_db.search_for_gene_functions([term], requested_sources=self.annotation_sources, verbose=True)
                # gene id's of genes with the searched function
                genes_that_hit = [i[0] for i in search_report]
                gene_caller_ids_of_interest.extend(genes_that_hit)

                self.search_term_to_gene_id_hits_dict[term] = set(genes_that_hit)

                self.targets.append('functions')
                self.sources = contigs_db.gene_function_call_sources
                counter += 1

        # Multiple sources could annotate the same gene, so make sure the list is unique
        self.gene_caller_ids_of_interest = set(gene_caller_ids_of_interest)

        if len(self.gene_caller_ids_of_interest):
            run.info('Matching genes',
                     '%d genes matched your search' % len(self.gene_caller_ids_of_interest),
                     mc='green', nl_after=1)


    def process(self, skip_init=False):
        if not skip_init:
            self.init()

        if not len(self.gene_caller_ids_of_interest):
            self.run.warning("There aren't any gene calls that match to the criteria you provided to anvi'o "
                             "export locus magic. Is this yet another case of you did everything right "
                             "yet anvi'o failed you? If that's the case, let us know :( This class will quietly "
                             "kill this process without reporting any error since a lack of hit may be the "
                             "expected outcome of some weird processes somewhere.")
            return

        self.contigs_db = dbops.ContigsSuperclass(self.args, r=self.run_object)
        self.contigs_db.init_functions()

        # Here we will differentiate between being in default-mode OR flank-mode. If in
        # default-mode, we will iterate through self.gene_caller_ids_of_interest and cut out the
        # locus X amount of genes above and below the gene_callers_id based on the --num-genes given
        # by the user. If in flank-mode, we will only export 1 locus based on the flanking gene
        # caller ids provided (--gene-caller-ids).

        # default
        if not self.is_in_flank_mode:
            counter = 1
            for gene_callers_id in self.gene_caller_ids_of_interest:
                self.run.warning(None,
                                 header="Exporting locus %d of %d" % (counter, len(self.gene_caller_ids_of_interest)),
                                 nl_after=0)

                output_path_prefix = os.path.join(self.output_dir, "%s_%.4d" % (self.output_file_prefix, counter))

                self.export_locus(output_path_prefix, gene_callers_id)

                counter += 1

            ##################
            # This is where anvi'o should print a warning message that the exported loci
            # are overlapping. Please see this issue for a reproducible example:
            # https://github.com/merenlab/anvio/issues/1228

            # Print warning if exported loci overlap on contig
            # from itertools import combinations

            # counter = 0
            # for combo in combinations(list(self.gene_caller_ids_of_interest), 2):  # 2 for pairs
            #     gene_1_1 = combo[0] - self.num_genes_list[0]
            #     gene_2_1 = combo[0] + self.num_genes_list[1]
            #     first_gene_of_the_block_1 = min(gene_1_1, gene_2_1)
            #     last_gene_of_the_block_1 = max(gene_1_1, gene_2_1)

            #     gene_1_2 = combo[1] - self.num_genes_list[0]
            #     gene_2_2 = combo[1] + self.num_genes_list[1]
            #     first_gene_of_the_block_2 = min(gene_1_2, gene_2_2)
            #     last_gene_of_the_block_2 = max(gene_1_2, gene_2_2)

            #     if first_gene_of_the_block_2 < first_gene_of_the_block_1 < last_gene_of_the_block_2:
            #         self.run.warning(None,
            #                         header="OVERLAPPPING!",
            #                         nl_after=0)

            #     elif first_gene_of_the_block_2 < last_gene_of_the_block_1 < last_gene_of_the_block_2:
            #         self.run.warning(None,
            #                         header="OVERLAPPPING!",
            #                         nl_after=0)
            #     counter += 1

            ##################

        # flank-mode
        elif self.is_in_flank_mode:
            self.run.warning(None,
                             header="Exporting locus 1 of 1",
                             nl_after=0)

            output_path_prefix = os.path.join(self.output_dir, self.output_file_prefix)
            gene_caller_ids_flank_pair = list(self.gene_caller_ids_of_interest)
            self.export_locus(output_path_prefix, None, gene_caller_ids_flank_pair)


    def export_locus(self, output_path_prefix, gene_callers_id=None, gene_caller_ids_flank_pair=None):
        """
        This function takes gene_callers_id or gene_caller_ids_flank_pair, and exports a contigs database.

        If gene_callers_id is provided, export_locus use --num-genes to  cut X above and below the
        gene_callers_id. If gene_caller_ids_flank_pair is provided, export_locus will cut out the locus
        between the pair provided. It is REQUIRED that gene_caller_ids_flank_pair is a pair of gen-callers-ID.

        Output path prefix should be unique for every export locus call. If the prefix you provide
        looks like this:

            >>> output_path_prefix = '/path/to/dir/file_name_prefix'

        the output files will be stored as this:

            >>> '/path/to/dir/file_name_prefix.fa'
            >>> '/path/to/dir/file_name_prefix.db'

        """

        if gene_callers_id is not None and not isinstance(gene_callers_id, int):
            raise ConfigError("The gene_caller_id must be an integer.")
        if gene_callers_id is not None and gene_caller_ids_flank_pair is not None:
            raise ConfigError("You can only provide the gene_callers_id or gene_caller_id_pair (with a , delimiter).")
        elif gene_callers_id is None and gene_caller_ids_flank_pair is None:
            raise ConfigError("You must provide at least 1 of the following: gene_callers_id or gene_caller_id_pair (with a , delimiter).")

        if self.is_in_flank_mode:
            if not isinstance(gene_caller_ids_flank_pair, list):
                raise ConfigError("The gene_caller_ids_flank_pair must be integers.")
            if [g for g in gene_caller_ids_flank_pair if not isinstance(g, int) or g < 0]:
                raise ConfigError("Both gene-caller_ids inputs must be integers!")

            if len(gene_caller_ids_flank_pair) == 1:
                raise ConfigError("You are in flank-mode, and anvi'o only found %d gene caller id(s). "
                                  "Anvi'o cannot handle this because flank-mode needs a pair of gene caller id's "
                                  "to cut out a locus (i.e., only a pair of flanking genes)! This most likely occured because 1 of your "
                                  "search-terms was not found the functions of the CONTIGS.db. Please try again with another "
                                  "search-term :)" % (len(self.gene_caller_ids_of_interest)))
            elif len(gene_caller_ids_flank_pair) > 2:
                raise ConfigError("You are in flank-mode, and anvi'o found %d total gene caller id's from the search-terms provided. "
                                  "Anvi'o cannot handle this because flank-mode needs a pair of gene caller id's "
                                  "to cut out a locus (i.e., only a pair of flanking genes)! Here are the gene caller ids anvi'o found "
                                  "from the search-terms %s: %s and %s: %s. Please use `anvi-export-functions` on your CONTIGS.db, locate "
                                  "these gene caller id's, then confirm the correct flanking gene caller ids. Anvi'o recommends you "
                                  "use the --gene caller ids flag to specify the specific pair gene caller ids you need to cut out the locus "
                                  "so there are no more mix ups. On the other hand, if you are trying to extract multiple loci from a genome "
                                  "using the same flanking genes, anvi'o cannot currently handle this in --flank-mode. If this functionality "
                                  "is necessary for your analysis, please make an issue on github and we will address it." % (len(self.gene_caller_ids_of_interest),
                                                                       str(self.search_term[0]),
                                                                       str(self.search_term_to_gene_id_hits_dict[self.search_term[0]]),
                                                                       str(self.search_term[1]),
                                                                       str(self.search_term_to_gene_id_hits_dict[self.search_term[1]]),
                                                                       ))

            gene_caller_ids = list(gene_caller_ids_flank_pair)
            gene_callers_id = gene_caller_ids[0] # just for getting contig name from contigDB


        if os.path.isdir(output_path_prefix):
            raise ConfigError("Output path prefix can't be a directory name...")

        filesnpaths.is_output_file_writable(output_path_prefix + '.fa')

        # if not already initiated, re-initiate contigsDB
        if not self.contigs_db:
            self.contigs_db = dbops.ContigsSuperclass(self.args, r=self.run_object)
            self.contigs_db.init_functions()

        # Query for gene_call, contig_name, and genes_in_contig_sorted
        gene_call = self.contigs_db.genes_in_contigs_dict[gene_callers_id]

        # if in flank mode check that both search terms are NOT on the same contig
        if self.is_in_flank_mode:
            contig_name_1 = self.contigs_db.genes_in_contigs_dict[gene_caller_ids[0]]['contig']
            contig_name_2 = self.contigs_db.genes_in_contigs_dict[gene_caller_ids[1]]['contig']
            if contig_name_1 != contig_name_2:
                raise ConfigError(f"Soooooo it turns out that the flanking genes you picked "
                                  f"are found on two separate contigs: {contig_name_1} and {contig_name_2}. That means we can't prepare "
                                  f"a smaller piece of DNA for you :/")
            
        contig_name = self.contigs_db.genes_in_contigs_dict[gene_callers_id]['contig']

        # Sort by gene start position
        genes_in_contig_sorted = sorted(list(self.contigs_db.contig_name_to_genes[contig_name]), key=lambda tup: tup[1])

        # Here we run Default-mode and cut out the locus based on the anchor gene_callers_id (index) and
        # cut + self.num_genes_list[0] and - self.num_genes_list[1] around it
        D = lambda: 1 if gene_call['direction'] == 'f' else -1

        if not self.is_in_flank_mode:
            counter = 0
            for gcid,start,stop in genes_in_contig_sorted:
                if gcid == gene_callers_id:
                    anchor_gene_index = counter
                counter = counter + 1


            gene_1 = anchor_gene_index - self.num_genes_list[0] * D()
            gene_2 = anchor_gene_index + self.num_genes_list[1] * D()
            first_gene_of_the_block = min(gene_1, gene_2)
            last_gene_of_the_block = max(gene_1, gene_2)

        if self.is_in_flank_mode:
            counter = 0
            anchor_gene_index_flank_pair = []
            for gene in genes_in_contig_sorted:
                if gene[0] in gene_caller_ids:
                    anchor_gene_index_flank_pair.append(counter)
                counter = counter + 1
            first_gene_of_the_block, last_gene_of_the_block = sorted(anchor_gene_index_flank_pair)
            
        # Print out locus info for user
        self.run.info("Contig name", contig_name)
        self.run.info("Contig length", self.contigs_db.contigs_basic_info[contig_name]['length'])
        self.run.info("Num genes in contig", len(genes_in_contig_sorted))
        self.run.info("Target gene call", gene_callers_id)
        self.run.info("Target gene direction", "Forward" if D() == 1 else "Reverse", mc = 'green' if D() == 1 else 'red')

        # getting gene indexes for the first and last genes in the contig
        last_gene_in_contig = len(genes_in_contig_sorted) - 1
        first_gene_in_contig = 0

        premature = False
        if last_gene_of_the_block > last_gene_in_contig:
            last_gene_of_the_block = last_gene_in_contig
            premature = True

        if first_gene_of_the_block < first_gene_in_contig:
            first_gene_of_the_block = first_gene_in_contig
            premature = True

        if premature and self.remove_partial_hits:
            self.run.info_single("A premature locus is found .. the current configuration says 'skip'. Skipping.", mc="red", nl_before=1)
            return
        elif premature and not self.remove_partial_hits:
            self.run.info_single("A premature locus is found .. the current configuration says 'whatevs'. Anvi'o will continue.", mc="yellow", nl_before=1, nl_after=1)

        self.run.info("First and last gene of the locus ", "%d and %d" % (first_gene_of_the_block, last_gene_of_the_block))

        # convert gene position indexes to gene-callers-ids
        first_gene_of_the_block_gene_callers_id = genes_in_contig_sorted[first_gene_of_the_block][0]
        last_gene_of_the_block_gene_callers_id = genes_in_contig_sorted[last_gene_of_the_block][0]

        locus_start = self.contigs_db.genes_in_contigs_dict[first_gene_of_the_block_gene_callers_id]['start']
        locus_stop = self.contigs_db.genes_in_contigs_dict[last_gene_of_the_block_gene_callers_id]['stop']

        # being a performance nerd here yes
        contig_sequence = db.DB(self.input_contigs_db_path, None, ignore_version=True, skip_rowid_prepend=True) \
                            .get_some_rows_from_table(t.contig_sequences_table_name,
                                                      where_clause="contig='%s'" % contig_name)[0][1]

        # Extract the locus!
        locus_sequence = contig_sequence[locus_start:locus_stop]

        # here we will create a gene calls dict for genes that are specific to our locus. since we trimmed
        # the contig sequence to the locus of interest, we will have to adjust start and stop positions of
        # genes in the gene calls dict.
        locus_gene_calls_dict = {}
        for g in range(first_gene_of_the_block, last_gene_of_the_block + 1):
            gci = genes_in_contig_sorted[g][0]
            locus_gene_calls_dict[gci] = copy.deepcopy(self.contigs_db.genes_in_contigs_dict[gci])
            excess = self.contigs_db.genes_in_contigs_dict[first_gene_of_the_block_gene_callers_id]['start']
            locus_gene_calls_dict[gci]['start'] -= excess
            locus_gene_calls_dict[gci]['stop'] -= excess

        self.run.info("Locus gene call start/stops excess (nts)", excess)

        if D() != 1 and self.reverse_complement_if_necessary:
            reverse_complement = True
        else:
            reverse_complement = False

        self.run.info('Reverse complementing everything', reverse_complement, mc='green')

        # report a stupid FASTA file.
        if self.include_fasta_output:
            fasta_file_path = output_path_prefix + ".fa"

            self.run.info("Output FASTA file", fasta_file_path)
            with open(fasta_file_path, 'w') as f:
                locus_header = contig_name + ' ' + \
                               '|'.join(['target:%s' % ','.join(self.targets),
                                         'sources:%s' % ','.join(self.sources),
                                         'query:%s' % self.search_term or 'None',
                                         'hit_contig:%s' % contig_name,
                                         'hit_gene_callers_id:%s' % str(gene_callers_id),
                                         'project_name:%s' % self.contigs_db.a_meta['project_name'].replace(' ', '_').replace("'", '_').replace('"', '_'),
                                         'locus:%s,%s' % (str(first_gene_of_the_block), str(last_gene_of_the_block)),
                                         'nt_positions_in_contig:%s:%s' % (str(locus_start), str(locus_stop)),
                                         'premature:%s' % str(premature),
                                         'reverse_complemented:%s' % str(reverse_complement)])

                f.write('>%s\n' % locus_header)
                f.write('%s\n' % utils.rev_comp(locus_sequence) if reverse_complement else locus_sequence)

        # report a fancy anvi'o contigs database
        self.store_locus_as_contigs_db(contig_name,
                                       locus_sequence,
                                       locus_gene_calls_dict,
                                       output_path_prefix,
                                       reverse_complement)


    def store_locus_as_contigs_db(self, contig_name, sequence, gene_calls, output_path_prefix, reverse_complement=False):
        """Generates a contigs database and a blank profile for a given locus"""

        temporary_files = []

        # dealing with some output file business.
        E = lambda e: output_path_prefix + e
        locus_output_db_path = E(".db")
        locus_sequence_fasta = E("_sequence.fa")
        locus_external_gene_calls = E("_external_gene_calls.txt")
        temporary_files.extend([locus_external_gene_calls, locus_sequence_fasta])

        # we will generate a blank profile database at the end of this. let's get the directory
        # business sorted.
        profile_output_dir = output_path_prefix + '-PROFILE'
        filesnpaths.check_output_directory(profile_output_dir)

        # sort out the contigs database output path
        filesnpaths.is_output_file_writable(locus_output_db_path, ok_if_exists=False)

        # do we need to reverse complement this guy? if yes, we will take care of the contigs sequence and
        # gene calls here, and remember this for later.
        gene_calls_list = list(gene_calls.keys())
        if reverse_complement:
            sequence = utils.rev_comp(sequence)
            gene_calls, gene_caller_id_conversion_dict = utils.rev_comp_gene_calls_dict(gene_calls, sequence)
        else:
            gene_caller_id_conversion_dict = dict([(gene_calls_list[g], g) for g in range(0, len(gene_calls_list))])
            new_gene_calls = {}
            for g in range(0, len(gene_calls_list)):
                gene_call = copy.deepcopy(gene_calls[gene_calls_list[g]])
                new_gene_calls[g] = gene_call
            gene_calls = new_gene_calls


        # write the sequence as a temporary FASTA file since the design of ContigsDatabase::create
        # will work seamlessly with this approach:
        with open(locus_sequence_fasta, 'w') as f:
            f.write('>%s\n%s\n' % (contig_name, sequence))

        # similarly, here we will store external gene calls so there will be no gene calling during
        # the generation of the contigs database
        headers = ['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version']
        utils.store_dict_as_TAB_delimited_file(gene_calls, locus_external_gene_calls, headers=headers)

        # this is where magic happens. we ask anvi'o to create a contigs database for us.
        args = argparse.Namespace(contigs_fasta=locus_sequence_fasta,
                                  project_name=os.path.basename(output_path_prefix),
                                  split_length=sys.maxsize,
                                  kmer_size=4,
                                  external_gene_calls=locus_external_gene_calls,
                                  ignore_internal_stop_codons=True)

        dbops.ContigsDatabase(locus_output_db_path, run=self.run_object).create(args)

        # while we are at it, here we generate a blank profile, too. so visualization of the
        # new contigs database for debugging or other purposes through anvi'o.
        args = argparse.Namespace(blank_profile=True,
                                  contigs_db=locus_output_db_path,
                                  skip_hierarchical_clustering=False,
                                  output_dir=profile_output_dir,
                                  sample_name=os.path.basename(output_path_prefix))
        profiler.BAMProfiler(args, r=terminal.Run(verbose=False), p=self.progress)._run()

        # while we are at it, let's add a default collection to the resulting blank profile
        TablesForCollections(os.path.join(profile_output_dir, 'PROFILE.db'), run=terminal.Run(verbose=False), progress=self.progress).add_default_collection_to_db(contigs_db_path=locus_output_db_path)
 

        # so we have a contigs database! but there isn't much in it. the following where clause will
        # help us read from the tables of the original contigs database, and store it into the
        # new one throughout the following sections of the code.
        where_clause = "gene_callers_id in (%s)" % ', '.join(['"%d"' % g for g in gene_caller_id_conversion_dict])

        # a lousy anonymous function to read data from tables given the gene calls of interest
        R = lambda table_name: db.DB(self.input_contigs_db_path, None, ignore_version=True) \
                                              .get_some_rows_from_table_as_dict(table_name,
                                                                                where_clause=where_clause,
                                                                                error_if_no_data=False)

        G = lambda g: gene_caller_id_conversion_dict[g]

        ############################################################################################
        # DO FUNCTIONS
        ###########################################################################################
        function_calls = R(t.gene_function_calls_table_name)

        for entry_id in function_calls:
            function_calls[entry_id]['gene_callers_id'] = G(function_calls[entry_id]['gene_callers_id'])

        gene_function_calls_table = TableForGeneFunctions(locus_output_db_path, run=self.run_object)
        gene_function_calls_table.create(function_calls)

        self.run.info("Output contigs DB path", locus_output_db_path)
        self.run.info("Output blank profile DB path", os.path.join(profile_output_dir, 'PROFILE.db'))

        ############################################################################################
        # DO AMINO ACID SEQUENCES -- we are using external gene calls to generate the new contigs
        #                            database, but amino acid sequences are kept in a different table
        #                            and anvi'o checks whether provided gene calls resolve to amino
        #                            acid sequences with proper starts and stops. if not, it skips
        #                            them. but amino acid sequences for each gene call was stored
        #                            in the original contigs database, and the best practice is to
        #                            carry them into the new one. so here we will remove all data
        #                            from the amino acid sequences table in the new database, and
        #                            copy the contents from the original one.
        ############################################################################################
        amino_acid_sequences = R(t.gene_amino_acid_sequences_table_name)

        entries = [(gene_caller_id_conversion_dict[g], amino_acid_sequences[g]['sequence']) for g in amino_acid_sequences]

        locus_db = db.DB(locus_output_db_path, None, ignore_version=True, skip_rowid_prepend=True)
        locus_db._exec("DELETE FROM %s" % t.gene_amino_acid_sequences_table_name)
        locus_db.insert_many(t.gene_amino_acid_sequences_table_name, entries=entries)
        locus_db.disconnect()

        ############################################################################################
        # REMOVE TEMP FILES
        ###########################################################################################
        if anvio.DEBUG:
            self.run.info_single("Temp output files were kept for inspection due to --debug")
        else:
            [os.remove(f) for f in temporary_files]
