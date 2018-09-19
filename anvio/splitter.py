# -*- coding: utf-8
# pylint: disable=line-too-long
"""The library to split merged profiles into smaller profiles.

The default client of this library is under bin/anvi-split"""


import os

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.clusteringconfuguration import ClusteringConfiguration
from anvio.tables.kmers import KMerTablesForContigsAndSplits
from anvio.tables.collections import TablesForCollections


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


class ProfileSplitter:
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
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
        if profile_db.meta['blank']:
            raise ConfigError("The anvi-split workflow is not prepared to deal with blank profiles :/ Sorry!")

        if profile_db.meta['db_type'] != 'profile':
            raise ConfigError("Anvi'o was trying to split this profile, but it just realized that it is not a profile\
                               database. There is something wrong here.")
        profile_db.disconnect()

        # if this is not set false, the summarizer class attemts to remove the main output directory
        # upon initialization. not doing that is useful in this context since this allows multiple
        # anvi-split runs to work on bins in the same collection in parallel:
        self.args.delete_output_directory_if_exists = False

        self.summary = summarizer.ProfileSummarizer(self.args)
        self.summary.init()

        self.bin_names_of_interest = sorted(self.summary.bin_ids)
        if self.bin_name:
            if self.bin_name not in self.bin_names_of_interest:
                raise ConfigError("The bin name you wish to split from this profile databse is not in the collection. Busted!")
            else:
                self.bin_names_of_interest = [self.bin_name]


    def process(self):
        """This is the function that goes through each bin loaded in the class and proecesses them."""
        self.sanity_check()

        filesnpaths.gen_output_directory(self.output_directory)

        self.run.warning("Anvi'o is about to start splitting your bins into individual, self-contained anvi'o profiles. This\
                          is quite a tricky operation, and even if it finishes successfully, you must double check everyting\
                          in the resulting profiles to make sure things worked as expected. Although we are doing our best to\
                          test all these, variation between projects make it impossible to be 100% sure.")

        if self.skip_variability_tables:
            self.run.warning("Since you asked so nicely, anvi'o will not migrate variability table data into split profiles.")

        for bin_name in self.bin_names_of_interest:
            b = BinSplitter(bin_name, self.summary, self.args, run=self.run, progress=self.progress)
            b.do_contigs_db()
            b.do_profile_db()

            if self.summary.auxiliary_profile_data_available:
                b.do_auxiliary_profile_data()

        self.run.info('Num bins processed', len(self.bin_names_of_interest))
        self.run.info("Output directory", self.output_directory)


class BinSplitter(summarizer.Bin):
    def __init__(self, bin_name, summary_object, args, run=run, progress=progress):
        """A class to split a single bin from its parent.

        The class is not really useful without a summary object, but it makes logistic sense to keep it
        separate since the inheritance from anvio/summarizer.Bin is much easier and sane this way."""
        summarizer.Bin.__init__(self, summary_object, bin_name, run, progress)

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

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)
        self.clustering_configs = constants.clustering_configs['merged']
        self.database_paths = {'CONTIGS.db': os.path.abspath(self.contigs_db_path)}

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering\
                               while also asking it to enforce it.")

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
                    'kmer_contigs': ('contig', self.split_names),
                    'kmer_splits': ('contig', self.split_names),
                }

        self.migrate_data(tables, self.contigs_db_path, self.bin_contigs_db_path)

        self.progress.end()


    def do_auxiliary_profile_data(self):
        self.progress.new('Splitting "%s"' % self.bin_id)
        self.progress.update('Subsetting the auxiliary data (for profile db)')

        new_auxiliary_profile_data_path = dbops.get_auxiliary_data_path_for_profile_db(self.bin_profile_db_path)
        parent_auxiliary_profile_data_path = self.summary.auxiliary_data_path

        bin_profile_auxiliary = auxiliarydataops.AuxiliaryDataForSplitCoverages(new_auxiliary_profile_data_path,
                                                                                self.contigs_db_hash,
                                                                                create_new=True)

        parent_profile_auxiliary = auxiliarydataops.AuxiliaryDataForSplitCoverages(parent_auxiliary_profile_data_path,
                                                                                   self.summary.a_meta['contigs_db_hash'])

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

        # update some values
        bin_profile_db.db.update_meta_value('contigs_db_hash', self.contigs_db_hash)
        bin_profile_db.db.update_meta_value('available_clusterings', None)
        bin_profile_db.db.update_meta_value('sample_id', self.bin_id)

        # setup the filtering rules for migrating data:
        tables = {}

        # this is to deal with merge atomic data tables that are stored in merged profiles.
        # they are being created on the fly during merge, so bin_profile_db.touch() did not
        # create them, and we have to do it here ourselves. while creating them in the target
        # db, we will also populate the tables dictionary for data migration::
        sample_names = self.summary.p_meta['samples']
        if merged:
            for table_name in t.atomic_data_table_structure[1:-1]:
                for target in ['splits', 'contigs']:
                    new_table_name = '_'.join([table_name, target])
                    new_table_structure = ['contig'] + sample_names + ['__parent__']
                    new_table_types = ['text'] + ['numeric'] * len(sample_names) + ['text']
                    bin_profile_db.db.create_table(new_table_name, new_table_structure, new_table_types)

                    tables[new_table_name] = ('contig', self.split_names)
        else:
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            table_structure = profile_db.db.get_table_structure('atomic_data_contigs')
            table_types = profile_db.db.get_table_column_types('atomic_data_contigs')
            for table_name in ['atomic_data_splits', 'atomic_data_contigs']:
                new_table_structure = profile_db.db.get_table_structure(table_name)
                bin_profile_db.db.create_table(table_name, table_structure, table_types)

                tables[table_name] = ('contig', self.split_names)


        # we need to migrate these guys, too. unless we don't need to... if we are migrating,
        # the values in the self table are already accurate. if we are skipping, regardless
        # of what the values were, we will set the absolut correct ones.
        if self.skip_variability_tables:
            bin_profile_db.db.update_meta_value('SNVs_profiled', False)
            bin_profile_db.db.update_meta_value('SCVs_profiled', False)
        else:
            tables[t.variable_nts_table_name] = ('split_name', self.split_names)
            tables[t.variable_codons_table_name] = ('corresponding_gene_call', self.gene_caller_ids)

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


    def is_dbs_identical(self, source_db_path, target_db_path):
        """Check whether the two dbs have identical table names"""

        source_db = db.DB(source_db_path, None, ignore_version=True)
        target_db = db.DB(target_db_path, None, ignore_version=True)

        source_tables = set(source_db.get_table_names())
        target_tables = set(target_db.get_table_names())

        if not source_tables == target_tables:
            self.progress.end()

            diff = source_tables.difference(target_tables)

            raise ConfigError("Something went wrong during subsetting :/ Table names in the parent db (%s) and the child\
                               db (%s) does not seem to be identical. The following tables are found in the source, but\
                               missing in the target database: '%s'" % (source_db_path, target_db_path, ', '.join(diff)))

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
            source_db = db.DB(source_db_path, None, ignore_version=True)
            data = source_db.get_some_rows_from_table(table_name, where_clause)
            source_db.disconnect()

            if not len(data):
                continue

            self.progress.update("Table '%s' .. writing data" % table_name)
            target_db = db.DB(target_db_path, None, ignore_version=True)
            target_db._exec_many('''INSERT INTO %s VALUES(%s)''' % (table_name, ','.join(['?'] * len(data[0]))), data)
            target_db.disconnect()

        # make sure things are OK
        self.is_dbs_identical(source_db_path, target_db_path)


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


    def is_hierarchical_clustering_for_bin_OK(self):
        skip_hierarchical_clustering = self.skip_hierarchical_clustering

        if self.num_splits > self.max_num_splits_for_hierarchical_clustering and not self.enforce_hierarchical_clustering:
            self.run.warning("It seems you have more than %s splits in this particular bin. This is the\
                              soft limit for anvi'o to attempt to create a hierarchical clustering of your splits\
                              (which becomes the center tree in all anvi'o displays). If you want a hierarchical\
                              clustering to be done anyway, you can re-run the splitting process only for this bin\
                              by adding these parameters to your run: '--bin-id %s --enforce-hierarchical-clustering'.\
                              If you feel like you are lost, don't hesitate to get in touch with anvi'o developers." \
                                                        % (pp(self.max_num_splits_for_hierarchical_clustering), self.bin_id))
            skip_hierarchical_clustering = True

        if self.num_splits > self.max_num_splits_for_hierarchical_clustering and self.enforce_hierarchical_clustering:
            self.run.warning("Becasue you have used the flag `--enforce-hierarchical-clustering`, anvi'o will attempt\
                              to create a hierarchical clustering of your %s splits for this bin. It may take a bit of\
                              time, and it is not even anvi'o's fault, you know  :/" \
                                                        % pp(self.max_num_splits_for_hierarchical_clustering))

        return skip_hierarchical_clustering


