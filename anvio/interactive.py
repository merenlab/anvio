# -*- coding: utf-8
# pylint: disable=line-too-long
"""The module that curates data for the interactive interface"""

import os
import sys
import numpy
import textwrap
from ete3 import Tree

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.clusteringconfuguration import ClusteringConfiguration
from anvio.dbops import ProfileSuperclass, ContigsSuperclass, PanSuperclass, SamplesInformationDatabase, TablesForStates, ProfileDatabase
from anvio.dbops import is_profile_db_and_contigs_db_compatible, is_profile_db_and_samples_db_compatible, get_description_in_db
from anvio.dbops import get_default_clustering_id, get_split_names_in_profile_db
from anvio.completeness import Completeness
from anvio.errors import ConfigError, RefineError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


progress = terminal.Progress()
run = terminal.Run()


class Interactive(ProfileSuperclass, PanSuperclass, ContigsSuperclass):
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args, external_clustering=None):
        self.args = args
        self.views = {}
        self.states_table = None
        self.p_meta = {}
        self.title = 'Unknown Project'

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.mode = A('mode')
        self.pan_db_path = A('pan_db')
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.genomes_storage_path = A('genomes_storage')
        self.collection_name = A('collection_name')
        self.manual_mode = A('manual_mode')
        self.split_hmm_layers = A('split_hmm_layers')
        self.taxonomic_level = A('taxonomic_level') or 't_genus'
        self.additional_layers_path = A('additional_layers')
        self.additional_view_path = A('additional_view')
        self.samples_information_db_path = A('samples_information_db')
        self.view = A('view')
        self.fasta_file = A('fasta_file')
        self.view_data_path = A('view_data')
        self.tree = A('tree')
        self.title = A('title')
        self.output_dir = A('output_dir')
        self.show_views = A('show_views')
        self.state_autoload = A('state_autoload')
        self.collection_autoload = A('collection_autoload')
        self.read_only = A('read_only')
        self.show_states = A('show_states')
        self.skip_check_names = A('skip_check_names')
        self.list_collections = A('list_collections')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.skip_init_functions = A('skip_init_functions')
        self.skip_auto_ordering = A('skip_auto_ordering')
        self.debug = A('debug')
        self.bin_ids_file_path = A('bin_ids_file')
        self.bin_id = A('bin_id')
        self.collection_name = A('collection_name')

        if self.pan_db_path and self.profile_db_path:
            raise ConfigError("You can't set both a profile database and a pan database in arguments\
                                you send to this class. What are you doing?")

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        self.displayed_item_names_ordered = None
        self.auxiliary_profile_data_available = False

        self.samples_information_dict = {}
        self.samples_order_dict = {}
        self.samples_information_default_layer_order = {}

        self.additional_layers_dict = {}
        self.additional_layers_headers = []

        # make sure the mode will be set properly
        if self.collection_name and self.manual_mode:
            raise ConfigError("You can't anvi-interactive in manual mode with a collection name.")

        self.external_clustering = external_clustering

        self.collections = ccollections.Collections()

        ContigsSuperclass.__init__(self, self.args)
        self.init_splits_taxonomy(self.taxonomic_level)

        if self.samples_information_db_path:
            samples_information_db = SamplesInformationDatabase(self.samples_information_db_path)
            self.samples_information_dict, self.samples_order_dict = samples_information_db.get_samples_information_and_order_dicts()
            self.samples_information_default_layer_order = samples_information_db.get_samples_information_default_layer_order()
            samples_information_db.disconnect()

        if self.contigs_db_path:
            self.completeness = Completeness(self.contigs_db_path)
            self.collections.populate_collections_dict(self.contigs_db_path)
        else:
            self.completeness = None

        # make sure we are not dealing with apples and oranges here.
        if self.contigs_db_path and self.profile_db_path:
            is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.P = lambda x: os.path.join(self.p_meta['output_dir'], x)
        self.cwd = os.getcwd()

        # here is where the big deal stuff takes place:
        if not self.mode and self.manual_mode:
            self.mode = 'manual'
            self.run.info('Mode', self.mode, mc='red')
            self.load_manual_mode()
        elif self.mode == 'refine':
            self.load_full_mode()
            self.load_refine_mode()
        elif self.mode == 'pan':
            self.load_pan_mode()
        elif self.collection_name or self.list_collections:
            self.mode = 'collection'
            self.run.info('Mode', self.mode, mc='green')
            self.load_collection_mode()
        else:
            self.mode = 'full'
            self.load_full_mode()

        # make sure the samples information database, if there is one, is in fact compatible with the profile database
        # the reason we are doing this here is because when we are in 'self.manual_mode', the self.p_meta['samples'] is
        # being filled within the self.load_manual_mode function based on the headers of the view data.
        if self.profile_db_path and self.samples_information_db_path:
            is_profile_db_and_samples_db_compatible(self.profile_db_path, self.samples_information_db_path, manual_mode_exception=self.manual_mode)

        if self.external_clustering:
            self.p_meta['clusterings'] = self.clusterings = self.external_clustering['clusterings']
            self.p_meta['available_clusterings'] = list(self.clusterings.keys())
            self.p_meta['default_clustering'] = self.external_clustering['default_clustering']

        if not self.state_autoload and 'default' in self.states_table.states:
            self.state_autoload = 'default'

        if not self.collection_autoload and 'default' in self.collections.collections_dict:
            self.collection_autoload = 'default'

        if not self.p_meta['clusterings']:
            if self.p_meta['db_type'] == 'pan':
                raise ConfigError("This pangenome (which you gracefully named as '%s') does not seem to have any hierarchical\
                                   clustering of protein clusters (PCs) in it. Maybe you skipped the clustering step, maybe\
                                   anvi'o skipped it on your behalf because you had too many PCs or something. Regardless of\
                                   who did what, you don't get to display your pangenome at this particular instance. Sorry :/" \
                                                            % (self.p_meta['project_name']))
            else:
                if self.p_meta['merged']:
                    raise ConfigError("This merged profile database does not seem to have any hierarchical clustering\
                                       of splits that is required by the interactive interface. It may have been generated\
                                       by anvi-merge with the `--skip-hierarchical-clustering` flag, or hierarchical\
                                       clustering step may have been skipped by anvi-merge because you had too many stplits\
                                       to get the clustering in a reasonable amount of time. Please read the help menu for\
                                       anvi-merge, and/or refer to the tutorial: \
                                       http://merenlab.org/2015/05/01/anvio-tutorial/#clustering-during-merging")
                else:
                    raise ConfigError("This single profile database does not seem to have any hierarchical clustering\
                                       that is required by the interactive interface. You must use `--cluster-contigs`\
                                       flag for single profiles to access to this functionality. Please read the help\
                                       menu for anvi-profile, and/or refer to the tutorial.")

        # self.displayed_item_names_ordered is going to be the 'master' names list. everything else is going to
        # need to match these names:
        default_clustering = self.p_meta['clusterings'][self.p_meta['default_clustering']]
        if 'newick' in default_clustering:
            self.displayed_item_names_ordered = utils.get_names_order_from_newick_tree(default_clustering['newick'], reverse=True)
        elif 'basic' in default_clustering:
            self.displayed_item_names_ordered = default_clustering['basic']
        else:
            raise ConfigError("There is something wrong here, and anvi'o needs and adult :( Something that should\
                               never happen happened. The default clustering does not have a basic or newick type.")

        # now we knot what splits we are interested in (self.displayed_item_names_ordered), we can get rid of all the
        # unnecessary splits stored in views dicts.
        self.prune_view_dicts()

        # we are going to iterate the newick trees, and make sure that internal nodes have labels
        for clustering_name in self.p_meta['clusterings']:
            if 'newick' in self.p_meta['clusterings'][clustering_name]:
                tree = Tree(self.p_meta['clusterings'][clustering_name]['newick'], format=1)

                node_counter = 0
                for node in tree.traverse():
                    if node.name == "":
                        node.name = "Int_%d" % node_counter
                        node_counter += 1

                if node_counter > 0:
                    # if we did not changed any branch name there is no need to spend time for 
                    # serialization back to newick
                    self.p_meta['clusterings'][clustering_name]['newick'] = tree.write(format=1) 

        # if there are any HMM search results in the contigs database other than 'singlecopy' sources,
        # we would like to visualize them as additional layers. following function is inherited from
        # Contigs DB superclass and will fill self.hmm_searches_dict if appropriate data is found in
        # search tables:
        if self.mode == 'full':
            self.init_non_singlecopy_gene_hmm_sources(self.displayed_item_names_ordered, return_each_gene_as_a_layer=self.split_hmm_layers)

        # take care of additional layers, and update ordering information for items
        if self.additional_layers_path:
            filesnpaths.is_file_tab_delimited(self.additional_layers_path)

            self.additional_layers_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_layers_path, dict_to_append=self.additional_layers_dict, assign_none_for_missing=True)
            self.additional_layers_headers = self.additional_layers_headers + utils.get_columns_of_TAB_delim_file(self.additional_layers_path)

        self.check_names_consistency()
        self.auto_order_items_based_on_additional_layers_data()
        self.convert_view_data_into_json()


    def auto_order_items_based_on_additional_layers_data(self):
        if self.skip_auto_ordering:
            return

        self.progress.new('Processing additional layer data for ordering of splits (to skip: --skip-auto-ordering)')
        skipped_additional_data_layers = []
        # go through additional layers that are not of type `bar`.
        for layer in [additional_layer for additional_layer in self.additional_layers_headers if '!' not in additional_layer]:
            self.progress.update('for "%s" ...' % layer)
            layer_type = utils.get_predicted_type_of_items_in_a_dict(self.additional_layers_dict, layer)

            if layer_type == None:
                skipped_additional_data_layers.append(layer)
                continue

            item_layer_data_tuple = [(layer_type(self.additional_layers_dict[item][layer]), item) for item in self.additional_layers_dict]

            self.p_meta['available_clusterings'].append('>> %s:none:none' % layer)
            self.p_meta['clusterings']['>> %s' % layer] = {'basic': [i[1] for i in sorted(item_layer_data_tuple)]}

            self.p_meta['available_clusterings'].append('>> %s_(reverse):none:none' % layer)
            self.p_meta['clusterings']['>> %s_(reverse)' % layer] = {'basic': [i[1] for i in sorted(item_layer_data_tuple, reverse=True)]}

        self.progress.end()

        if len(skipped_additional_data_layers):
            self.run.warning("One or more of your additional data columns were completely empty. Like, they didn't have any data at all :/\
                              In the best case scenario you will see completely blank layers in your display. In the worst case scenario\
                              other things will break. Since you are a curious person, anvi'o thought you would like to know. These are\
                              the empty variables: %s." % ', '.join(['"%s"' % s for s in skipped_additional_data_layers]))


    def load_manual_mode(self):
        if self.contigs_db_path:
            raise ConfigError("When you want to use the interactive interface in manual mode, you must\
                                not use a contigs database.")

        if not self.profile_db_path:
            raise ConfigError("Even when you want to use the interactive interface in manual mode, you need\
                                to provide a profile database path. But you DO NOT need an already existing\
                                profile database, since anvi'o will generate an empty one for you. The profile\
                                database in this mode only used to read or store the 'state' of the display\
                                for visualization purposes, or to allow you to create and store collections.")

        # if the user is using an existing profile database, we need to make sure that it is not associated
        # with a contigs database, since it would mean that it is a full anvi'o profile database and should
        # not be included in manual operations.
        if filesnpaths.is_file_exists(self.profile_db_path, dont_raise=True):
            profile_db = ProfileDatabase(self.profile_db_path)
            if profile_db.meta['contigs_db_hash']:
                raise ConfigError("Well. It seems the profile database is associated with a contigs database,\
                                    which means using it in manual mode is not the best way to use it. Probably\
                                    what you wanted to do is to let the manual mode create a new profile database\
                                    for you. Simply type in a new profile database path (it can be a file name\
                                    that doesn't exist).")

        if not self.tree and not self.view_data_path:
            raise ConfigError("You must be joking Mr. Feynman. No tree file, and no data file? What is it that\
                               anvi'o supposed to visualize? :(")

        if not self.tree:
            self.run.warning("You haven't declared a tree file. Anvi'o will do its best to come up with an\
                              organization of your items.")

        if self.view:
            raise ConfigError("You can't use '--view' parameter when you are running the interactive interface\
                                in manual mode")

        if self.show_views:
            raise ConfigError("Sorry, there are no views to show in manual mode :/")

        if self.show_states:
            raise ConfigError("Sorry, there are no states to show in manual mode :/")

        if self.tree:
            filesnpaths.is_file_exists(self.tree)
            newick_tree_text = ''.join([l.strip() for l in open(os.path.abspath(self.tree)).readlines()])
            item_names = utils.get_names_order_from_newick_tree(newick_tree_text)
        else:
            item_names = utils.get_column_data_from_TAB_delim_file(self.view_data_path, column_indices=[0])[0][1:]

        # try to convert item names into integer values for proper sorting later. it's OK if it does
        # not work.
        try:
            item_names = [int(n) for n in item_names]
        except:
            pass

        view_data_path = os.path.abspath(self.view_data_path) if self.view_data_path else None
        self.p_meta['splits_fasta'] = os.path.abspath(self.fasta_file) if self.fasta_file else None
        self.p_meta['output_dir'] = None
        self.p_meta['views'] = {}
        self.p_meta['merged'] = True
        self.p_meta['default_view'] = 'single'
        self.default_view = self.p_meta['default_view']

        # set some default organizations of data:
        self.p_meta['clusterings'] = {'Alphabetical_(reverse):none:none': {'basic': [str(i) for i in sorted(item_names)]},
                                      'Alphabetical:none:none': {'basic': [str(i) for i in sorted(item_names, reverse=True)]}}
        self.p_meta['available_clusterings'] = ['Alphabetical_(reverse):none:none', 'Alphabetical:none:none']
        self.p_meta['default_clustering'] = self.p_meta['available_clusterings'][0]

        # if we have a tree, let's make arrangements for it:
        if self.tree:
            clustering_id = '%s:unknown:unknown' % filesnpaths.get_name_from_file_path(self.tree)
            self.p_meta['default_clustering'] = clustering_id
            self.p_meta['available_clusterings'].append(clustering_id)
            self.p_meta['clusterings'][clustering_id] = {'newick': newick_tree_text}

        if self.view_data_path:
            # sanity of the view data
            filesnpaths.is_file_tab_delimited(view_data_path)
            view_data_columns = utils.get_columns_of_TAB_delim_file(view_data_path, include_first_column=True)

            # load view data as the default view:
            self.views[self.default_view] = {'header': view_data_columns[1:],
                                             'dict': utils.get_TAB_delimited_file_as_dictionary(view_data_path)}
        else:
            # no view data is provided... it is only the tree we have. we will creaet a mock 'view data dict'
            # here using what is in the tree.
            ad_hoc_dict = {}
            for item in item_names:
                ad_hoc_dict[item] = {'names': str(item)}

            self.views[self.default_view] = {'header': ['names'],
                                             'dict': ad_hoc_dict}

        self.displayed_item_names_ordered = list(self.views[self.default_view]['dict'].keys())

        # we assume that the sample names are the header of the view data, so we might as well set it up:
        self.p_meta['samples'] = self.views[self.default_view]['header']

        # if we have an input FASTA file, we will set up the split_sequences and splits_basic_info dicts,
        # otherwise we will leave them empty
        self.splits_basic_info = {}
        self.split_sequences = None
        if self.p_meta['splits_fasta']:
            filesnpaths.is_file_fasta_formatted(self.p_meta['splits_fasta'])
            self.split_sequences = utils.get_FASTA_file_as_dictionary(self.p_meta['splits_fasta'])

            names_missing_in_FASTA = set(self.displayed_item_names_ordered) - set(self.split_sequences.keys())
            num_names_missing_in_FASTA = len(names_missing_in_FASTA)
            if num_names_missing_in_FASTA:
                raise ConfigError('Some of the names in your view data does not have corresponding entries in the\
                                    FASTA file you provided. Here is an example to one of those %d names that occur\
                                    in your data file, but not in the FASTA file: "%s"' % (num_names_missing_in_FASTA, names_missing_in_FASTA.pop()))

            # setup a mock splits_basic_info dict
            for split_id in self.displayed_item_names_ordered:
                self.splits_basic_info[split_id] = {'length': len(self.split_sequences[split_id]),
                                                    'gc_content': utils.get_GC_content_for_sequence(self.split_sequences[split_id])}

        # create a new, empty profile database for manual operations
        if not os.path.exists(self.profile_db_path):
            profile_db = ProfileDatabase(self.profile_db_path)
            profile_db.create({'db_type': 'profile', 'merged': True, 'contigs_db_hash': None, 'samples': ','.join(self.p_meta['samples'])})

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path)

        # also populate collections, if there are any
        self.collections.populate_collections_dict(self.profile_db_path)

        # read description from self table, if it is not available get_description function will return placeholder text
        self.p_meta['description'] = get_description_in_db(self.profile_db_path)

        if self.title:
            self.title = self.title


    def cluster_splits_of_interest(self):
        # clustering of contigs is done for each configuration file under static/clusterconfigs/merged directory;
        # at this point we don't care what those recipes really require because we already merged and generated
        # any data that may be required.
        clusterings = {}

        for config_name in self.clustering_configs:
            config_path = self.clustering_configs[config_name]

            config = ClusteringConfiguration(config_path, 
                                            self.input_directory, 
                                            db_paths={'CONTIGS.db': self.contigs_db_path, 
                                                      'PROFILE.db': self.profile_db_path},
                                            row_ids_of_interest=self.split_names_of_interest)

            try:
                clustering_id, newick = clustering.order_contigs_simple(config, progress=progress)
            except Exception as e:
                run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
                progress.end()
                continue

            clusterings[clustering_id] = {'newick': newick}

        run.info('available_clusterings', list(clusterings.keys()))

        return clusterings


    def load_refine_mode(self):
        self.split_names_of_interest = set([])
        self.is_merged = self.p_meta['merged']
        self.is_blank = self.p_meta['blank']
        self.clustering_configs = constants.clustering_configs['merged' if self.is_merged else 'blank' if self.is_blank else 'single']

        progress.new('Initializing')
        progress.update('Getting split names')

        split_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.bins = list(split_dict.keys())

        for split_names in list(split_dict.values()):
            self.split_names_of_interest.update(split_names)

        progress.end()

        # if the user updates the refinement of a single bin or bins, there shouldn't be multiple copies
        # of that stored in the database. so everytime 'store_refined_bins' function is called,
        # it will check this varlable and, (1) if empty, continue updating stuff in db store updates
        # in it, (2) if not empty, remove items stored in this variable from collections dict, and continue
        # with step (1). the starting point is of course self.bins. when the store_refined_bins function is
        # called the first time, it will read collection data for collection_name, and remove the bin(s) in
        # analysis from it before it stores the data:
        self.ids_for_already_refined_bins = self.bins

        self.input_directory = os.path.dirname(os.path.abspath(self.profile_db_path))

        self.run.info('Input directory', self.input_directory)
        self.run.info('Collection ID', self.collection_name)
        self.run.info('Number of bins', len(self.bins))
        self.run.info('Number of splits', len(self.split_names_of_interest))

        clusterings = self.cluster_splits_of_interest()
        default_clustering_class = constants.merged_default if self.is_merged else constants.single_default

        default_clustering_id = dbops.get_default_clustering_id(default_clustering_class, clusterings)

        self.clusterings = clusterings
        self.p_meta['clusterings'] = clusterings
        self.p_meta['available_clusterings'] = list(self.clusterings.keys())
        self.p_meta['default_clustering'] = default_clustering_id

        self.collections = ccollections.Collections()
        self.collections.populate_collections_dict(self.profile_db_path)

        bins = sorted(list(self.bins))

        if not self.args.title:
            self.title = textwrap.fill('Refining %s%s from "%s"' % (', '.join(bins[0:3]),
                                                      ' (and %d more)' % (len(bins) - 3) if len(bins) > 3 else '',
                                                      self.collection_name))

    def load_collection_mode(self):
        self.collections.populate_collections_dict(self.profile_db_path)

        if self.list_collections:
            self.collections.list_collections()
            sys.exit()

        if self.collection_name not in self.collections.collections_dict:
            raise ConfigError("%s is not a valid collection name. See a list of available ones with '--list-collections' flag" % self.collection_name)

        # learn whether HMMs were run and we have access to completion estimates, and initialize the hmm_access if they
        # did
        completion_redundancy_available = True
        completeness = Completeness(self.contigs_db_path)

        if not len(completeness.sources):
            self.run.warning('HMMs for single-copy core genes were not run for this contigs database. So you will not\
                              see completion / redundancy estimates in the collection mode as additional layers. SAD.')
            completion_redundancy_available = False

            self.progress.new('Accessing HMM hits')
            self.progress.update('...')
            self.hmm_access = hmmops.SequencesForHMMHits(self.contigs_db_path, sources=set(completeness.sources))
            self.progress.end()

        # we are about to request a collections dict that contains only split names that appear in the
        # profile database along with other info:
        self.collection, bins_info_dict, split_names_in_db_but_missing_in_collection = self.collections.get_trimmed_dicts(self.collection_name, get_split_names_in_profile_db(self.profile_db_path))

        # we will do something quite tricky here. first, we will load the full mode to get the self.views
        # data structure fully initialized based on the profile database. Then, we using information about
        # bins in the selected collection, we will create another views data structure, and replace it with
        # the one we have. that will be LOVELY.
        self.load_full_mode()

        # FIXME: `clusterings` should become `orderings` thorughout the code.
        self.p_meta['available_clusterings'] = []
        self.p_meta['clusterings'] = {}

        # we just cleared out all orderings the full mode added, let's make sure to add the
        # user tree if there is one.
        self.add_user_tree()

        # setting up a new view:
        views_for_collection = {}
        for view in self.views:
            v = self.views[view]

            d = {}
            d['table_name'] = v['table_name']
            d['header'] = [h for h in v['header'] if not h == '__parent__']
            d['dict'] = {}

            for bin_id in self.collection:
                d['dict'][bin_id] = {}
                for header in d['header']:
                     d['dict'][bin_id][header] = numpy.mean([v['dict'][split_name][header] for split_name in self.collection[bin_id]])

            clustering_id = ':'.join([view, self.distance, self.linkage])
            self.p_meta['available_clusterings'].append(clustering_id)
            self.p_meta['clusterings'][clustering_id] = {'newick': clustering.get_newick_tree_data_for_dict(d['dict'], distance=self.distance, linkage=self.linkage)}

            # clustering is done, we can get prepared for the expansion of the view dict
            # with new layers. Note that these layers are going to be filled later.
            if completion_redundancy_available:
                d['header'].extend(['percent_completion', 'percent_redundancy'])
            d['header'].extend(['bin_name', 'source'])

            views_for_collection[view] = d

        default_clustering_class = 'mean_coverage' if self.p_meta['merged'] else 'single'
        self.p_meta['default_clustering'] = get_default_clustering_id(default_clustering_class, self.p_meta['clusterings'])

        # replace self.views with the new view:
        self.views = views_for_collection

        # preparing a new 'splits_basic_info'
        basic_info_for_collection = {}
        for bin_id in self.collection:
            basic_info_for_collection[bin_id] = {'length': sum([self.splits_basic_info[s]['length'] for s in self.collection[bin_id]]),
                                                 'gc_content': numpy.mean([self.splits_basic_info[s]['gc_content'] for s in self.collection[bin_id]])}


        # replace it with the new one!
        self.splits_basic_info = basic_info_for_collection

        # additional layers for each bin, INCLUDING completion and redundancy estimates:
        self.progress.new('Additional layers for bins')
        num_bins = len(self.collection)
        current_bin = 1
        for bin_id in self.collection:
            self.progress.update('%d of %d :: %s ...' % (current_bin, num_bins, bin_id))

            if completion_redundancy_available:
                # get completeness estimate
                p_completion, p_redundancy, domain, domain_confidence, results_dict = completeness.get_info_for_splits(set(self.collection[bin_id]))

            for view in self.views:
                self.views[view]['dict'][bin_id]['bin_name'] = bin_id

                if completion_redundancy_available:
                    self.views[view]['dict'][bin_id]['percent_completion'] = p_completion
                    self.views[view]['dict'][bin_id]['percent_redundancy'] = p_redundancy

                self.views[view]['dict'][bin_id]['source'] = bins_info_dict[bin_id]['source']

            current_bin += 1
        self.progress.end()

        # let empty some variables to avoid confusion downstream
        self.hmm_sources_info = {}
        self.split_sequences = None
        self.splits_taxonomy_dict = {}
        self.genes_in_splits_summary_dict = {}
        self.displayed_item_names_ordered = sorted(self.views[self.default_view]['dict'].keys())

        # set the title:
        R = lambda x: x.replace('-', ' ').replace('_', ' ')
        self.title = "Collection '%s' for %s" % (R(self.collection_name), R(self.p_meta['sample_id']))


    def load_pan_mode(self):
        if not self.pan_db_path:
            raise ConfigError("So you want to display a pan genome without a pan database? Anvi'o is\
                                confused :/")

        PanSuperclass.__init__(self, self.args)

        self.init_protein_clusters()

        if not self.skip_init_functions:
            self.init_protein_clusters_functions()

        self.init_additional_layer_data()

        self.p_meta['clusterings'] = self.clusterings

        self.load_pan_views()
        self.default_view = self.p_meta['default_view']

        self.collections.populate_collections_dict(self.pan_db_path)

        # set title for the pangenome
        if self.title:
            self.title = self.title
        else:
            self.title = self.p_meta['project_name'].replace('-', ' ').replace('_', ' ')


    def load_full_mode(self):
        if not self.contigs_db_path:
            raise ConfigError("Anvi'o needs the contigs database to make sense of this run (or maybe you\
                                should use the `--manual` flag if that's what your intention).")

        if not self.profile_db_path:
            raise ConfigError("So you want to run anvi'o in full mode, but without a profile database?\
                                Well. This does not make any sense.")

        if not self.skip_init_functions:
            self.init_functions()

        ProfileSuperclass.__init__(self, self.args)

        # this is a weird place to do it, but we are going to ask ContigsSuperclass function to load
        # all the split sequences since only now we know the mun_contig_length that was used to profile
        # this stuff
        self.init_split_sequences(self.p_meta['min_contig_length'])

        self.collections.populate_collections_dict(self.profile_db_path)

        self.p_meta['self_path'] = self.profile_db_path
        self.p_meta['output_dir'] = os.path.join(os.getcwd(), os.path.dirname(self.profile_db_path))

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path)

        # load views from the profile database
        if self.p_meta['blank']:
            blank_dict = {}
            for split_name in self.splits_basic_info:
                blank_dict[split_name] = {'blank_view': 0}

            self.views['blank_view'] = {'header': ['blank_view'],
                                        'dict': blank_dict}
            self.default_view = 'blank_view'

        else:
            self.load_views()
            self.default_view = self.p_meta['default_view']

        # if the user wants to see available views, show them and exit.
        if self.show_views:
            run.warning('', header='Available views (%d)' % len(self.views), lc='green')
            for view in self.views:
                run.info(view,
                         'Via "%s" table' % self.views[view]['table_name'],
                         lc='crimson',
                         mc='green' if view == self.default_view else 'crimson')
            print()
            sys.exit()

        if self.show_states:
            run.warning('', header='Available states (%d)' % len(self.states_table.states), lc='green')
            for state in self.states_table.states:
                run.info(state,
                         'Last modified %s' % self.states_table.states[state]['last_modified'],
                         lc='crimson',
                         mc='crimson')
            print()
            sys.exit()

        # if the user has an additional view data, load it up into the self.views dict.
        if self.additional_view_path:
            filesnpaths.is_file_tab_delimited(self.additional_view_path)
            additional_view_columns = utils.get_columns_of_TAB_delim_file(self.additional_view_path)

            if not additional_view_columns[-1] == '__parent__':
                raise ConfigError("The last column of the additional view must be '__parent__' with the proper\
                                    parent information for each split.")

            column_mapping = [str] + [float] * (len(additional_view_columns) - 1) + [str]

            self.views['user_view'] = {'table_name': 'NA',
                                       'header': additional_view_columns,
                                       'dict': utils.get_TAB_delimited_file_as_dictionary(self.additional_view_path, column_mapping=column_mapping)}

        # if the user specifies a view, set it as default:
        if self.view:
            if not self.view in self.views:
                raise ConfigError("The requested view ('%s') is not available for this run. Please see\
                                          available views by running this program with --show-views flag." % self.view)

            self.default_view = self.view

        self.p_meta['clusterings'] = self.clusterings

        # add user tree if there is one
        self.add_user_tree()

        # set title
        if self.title:
            self.title = self.title
        else:
            self.title = self.p_meta['sample_id'].replace('-', ' ').replace('_', ' ')

        # do we have auxiliary data available?
        if not self.auxiliary_profile_data_available:
            summary_cp_available = os.path.exists(os.path.join(os.path.dirname(self.profile_db_path), 'SUMMARY.cp'))
            self.run.warning("Auxiliary data is not available; which means you will not be able to perform\
                              certain operations (i.e., the inspect menu in the interactive interface will\
                              not work, etc). %s" % ('' if not summary_cp_available else "Although, you have\
                              a SUMMARY.cp file in your work directory, which means you are working with an\
                              outdated anvi'o run. You can convert your SUMMARY.cp into an auxiliary data file\
                              by using `anvi-script-generate-auxiliary-data-from-summary-cp` script."))

        if self.state_autoload:
            if not self.state_autoload in self.states_table.states:
                raise ConfigError("The requested state ('%s') is not available for this run. Please see\
                                          available states by running this program with --show-states flag." % self.state_autoload)


    def add_user_tree(self):
        if self.tree:
            clustering_id = '%s:unknown:unknown' % filesnpaths.get_name_from_file_path(self.tree)
            if not self.p_meta['clusterings']:
                self.p_meta['default_clustering'] = clustering_id
                self.p_meta['available_clusterings'] = [clustering_id]
                self.p_meta['clusterings'] = {clustering_id: {'newick': open(os.path.abspath(self.tree)).read()}}
                run.info('Additional Tree', "Splits will be organized based on '%s'." % clustering_id)
            else:
                self.p_meta['clusterings'][clustering_id] = {'newick': open(os.path.abspath(self.tree)).read()}
                run.info('Additional Tree', "'%s' has been added to available trees." % clustering_id)


    def check_names_consistency(self):
        if self.skip_check_names:
            return

        splits_in_tree = set(self.displayed_item_names_ordered)
        splits_in_view_data = set(self.views[self.default_view]['dict'].keys())
        splits_in_database = set(self.split_sequences) if self.split_sequences else None
        splits_in_additional_view = set(self.views['user_view']['dict'].keys()) if self.additional_view_path else None

        splits_in_tree_but_not_in_view_data = splits_in_tree - splits_in_view_data
        splits_in_tree_but_not_in_database = splits_in_tree - splits_in_database if splits_in_database else set([])
        splits_in_additional_view_but_not_in_tree = splits_in_additional_view - splits_in_tree if splits_in_additional_view else set([])

        if splits_in_additional_view_but_not_in_tree:
            raise ConfigError("There are some split names in your additional view data file ('%s') that are missing from\
                                split names characterized in the database. There are in fact %d of them. For instance,\
                                here is a random split name that is in your additional view data, yet not in the database:\
                                '%s'. This is not going to work for anvi'o :/" \
                                    % (self.additional_view_path, len(splits_in_additional_view_but_not_in_tree), splits_in_additional_view_but_not_in_tree.pop()))

        if splits_in_tree_but_not_in_view_data:
            num_examples = 5 if len(splits_in_tree_but_not_in_view_data) >= 5 else len(splits_in_tree_but_not_in_view_data)
            example_splits_missing_in_view = [splits_in_tree_but_not_in_view_data.pop() for _ in range(0, num_examples)]
            raise ConfigError('Some split names found in your tree are missing in your view data. Hard to\
                                know what cuased this, but essentially your tree and your view does not\
                                seem to be compatible. Here is a couple of splits that appear in the tree\
                                but not in the view data: %s.' % ', '.join(example_splits_missing_in_view))

        if splits_in_tree_but_not_in_database:
            raise ConfigError('Some split names found in your tree are missing from your database. Hard to\
                                know why is this the case, but here is a couple of them: %s'\
                                    % ', '.join(list(splits_in_tree_but_not_in_database)[0:5]))

        if self.additional_layers_path:
            splits_in_additional_layers = set(sorted([l.split('\t')[0] for l in open(self.additional_layers_path).readlines()[1:]]))
            splits_only_in_additional_layers = []
            for split_name in splits_in_additional_layers:
                if split_name not in splits_in_tree:
                    splits_only_in_additional_layers.append(split_name)
            if len(splits_only_in_additional_layers):
                one_example = splits_only_in_additional_layers[-1]
                num_all = len(splits_only_in_additional_layers)
                run.warning("Some of the contigs in your addtional view data file does not appear to be in anywhere else.\
                             Additional layers file is not required to have data for all items (which means, there may be\
                             items in view dictionaries that are not in the additional view data file), however, finding\
                             items that are only in the additional layers file usually means trouble. Anvi'o will continue,\
                             but please go back and check your files if you think there may be something wrong. Here is a\
                             random item name that was only in your file: '%s'. And there were %d of them in total. So you\
                             are warned!" % (one_example, num_all))


    def prune_view_dicts(self):
        self.progress.new('Pruning view dicts')
        self.progress.update('...')
        splits_in_views = set(list(self.views.values())[0]['dict'].keys())
        splits_to_remove = splits_in_views - set(self.displayed_item_names_ordered)

        for view in self.views:
            self.progress.update('processing view "%s"' % view)
            for split_name in splits_to_remove:
                self.views[view]['dict'].pop(split_name)

        self.progress.end()


    def convert_view_data_into_json(self):
        '''This function's name must change to something more meaningful.'''

        for view in self.views:
            # here we will populate runinfo['views'] with json objects.
            view_dict = self.views[view]['dict']
            view_headers = self.views[view]['header']

            json_object = []

            # (1) set the header line with the first entry:
            json_header = ['contigs']

            # (2) add taxonomy, if exitsts:
            if len(self.splits_taxonomy_dict):
                json_header.extend(['taxonomy'])

            # (3) then add split summaries from contigs db, if exists
            if len(self.genes_in_splits_summary_dict):
                json_header.extend(self.genes_in_splits_summary_headers[1:])

            # (4) then add length and GC content IF we have sequences available
            if self.splits_basic_info:
                basic_info_headers = ['length', 'gc_content']
                json_header.extend(basic_info_headers)

            # (5) then add the view!
            json_header.extend(view_headers)

            # (6) then add 'additional' headers as the outer ring:
            if self.additional_layers_headers:
                json_header.extend(self.additional_layers_headers)

            # (7) finally add hmm search results
            if self.hmm_searches_header:
                json_header.extend([tpl[0] for tpl in self.hmm_searches_header])

            # (8) and finalize it (yay):
            json_object.append(json_header)

            for split_name in view_dict:
                # (1)
                json_entry = [split_name]

                # (2)
                if self.splits_taxonomy_dict:
                    if split_name in self.splits_taxonomy_dict:
                        json_entry.extend([self.splits_taxonomy_dict[split_name]])
                    else:
                        json_entry.extend([None])

                # (3)
                if self.genes_in_splits_summary_dict:
                    json_entry.extend([self.genes_in_splits_summary_dict[split_name][header] for header in self.genes_in_splits_summary_headers[1:]])

                # (4)
                if self.splits_basic_info:
                    json_entry.extend([self.splits_basic_info[split_name][header] for header in basic_info_headers])

                # (5) adding essential data for the view
                json_entry.extend([view_dict[split_name][header] for header in view_headers])

                # (6) adding additional layers
                json_entry.extend([self.additional_layers_dict[split_name][header] if split_name in self.additional_layers_dict else None for header in self.additional_layers_headers])

                # (7) adding hmm stuff
                if self.hmm_searches_dict:
                    if self.split_hmm_layers:
                        json_entry.extend([self.hmm_searches_dict[split_name][header] if split_name in self.hmm_searches_dict else None for header in [tpl[0] for tpl in self.hmm_searches_header]])
                    else:
                        json_entry.extend([len(self.hmm_searches_dict[split_name][header]) if split_name in self.hmm_searches_dict else 0 for header in [tpl[1] for tpl in self.hmm_searches_header]])

                # (8) send it along!
                json_object.append(json_entry)

            self.views[view] = json_object

    def store_refined_bins(self, refined_bin_data, refined_bins_info_dict):
        if 0 in [len(b) for b in list(refined_bin_data.values())]:
            raise RefineError('One or more of your bins have zero splits. If you are trying to remove this bin from your collection,\
                                this is not the right way to do it.')

        progress.new('Storing refined bins')
        progress.update('accessing to collection "%s" ...' % self.collection_name)

        collection_dict = self.collections.get_collection_dict(self.collection_name)
        bins_info_dict = self.collections.get_bins_info_dict(self.collection_name)
        
        progress.end()

        bad_bin_names = [b for b in collection_dict if (b in refined_bin_data and b not in self.ids_for_already_refined_bins)]
        if len(bad_bin_names):
            raise RefineError('%s of your bin names %s NOT unique, and already exist%s in the database. You must rename\
                                %s to something else: %s' % (
                                                              'One' if len(bad_bin_names) == 1 else len(bad_bin_names),
                                                              'is' if len(bad_bin_names) == 1 else 'are',
                                                              's' if len(bad_bin_names) == 1 else '',
                                                              'this one' if len(bad_bin_names) == 1 else 'these',
                                                              ', '.join(bad_bin_names)
                                                             ))

        # remove bins that should be updated in the database:
        for bin_id in self.ids_for_already_refined_bins:
            collection_dict.pop(bin_id)
            bins_info_dict.pop(bin_id)

        # zero it out
        self.ids_for_already_refined_bins = set([])

        if self.debug:
            run.info('collection from db', collection_dict)
            run.info('bins info from db', bins_info_dict)
            run.info_single('')

            run.info('incoming collection data', refined_bin_data)
            run.info('incoming bins info', refined_bins_info_dict)
            run.info_single('')

        for bin_id in refined_bin_data:
            collection_dict[bin_id] = refined_bin_data[bin_id]
            bins_info_dict[bin_id] = refined_bins_info_dict[bin_id]
            self.ids_for_already_refined_bins.add(bin_id)


        if self.debug:
            run.info('resulting collection', collection_dict)
            run.info('resulting bins info', bins_info_dict)
            run.info_single('')

        collections = dbops.TablesForCollections(self.profile_db_path)
        collections.append(self.collection_name, collection_dict, bins_info_dict)

        run.info_single('"%s" collection is updated!' % self.collection_name)


    def end(self):
        # FIXME: remove temp files and stuff
        pass
