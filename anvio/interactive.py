# -*- coding: utf-8
# pylint: disable=line-too-long
"""The module that curates data for the interactive interface"""

import os
import sys
import copy
import json
import numpy
import argparse
import textwrap
import pandas as pd

import anvio
import anvio.tables as t
import anvio.kegg as kegg
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.summarizer as summarizer
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.structureops as structureops
import anvio.variabilityops as variabilityops

from anvio.completeness import Completeness
from anvio.dbops import get_description_in_db
from anvio.tables.views import TablesForViews
from anvio.variabilityops import VariabilitySuper
from anvio.variabilityops import variability_engines
from anvio.dbops import get_default_item_order_name
from anvio.genomedescriptions import AggregateFunctions
from anvio.errors import ConfigError, RefineError, GenesDBError
from anvio.clusteringconfuguration import ClusteringConfiguration
from anvio.dbops import ProfileSuperclass, ContigsSuperclass, PanSuperclass, TablesForStates, ProfileDatabase

from anvio.tables.miscdata import (
    TableForItemAdditionalData,
    TableForLayerAdditionalData,
    TableForLayerOrders,
    TableForAminoAcidAdditionalData,
)
from anvio.tables.collections import TablesForCollections


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


progress = terminal.Progress()
run = terminal.Run()
pp = terminal.pretty_print


class Interactive(ProfileSuperclass, PanSuperclass, ContigsSuperclass):
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args, external_clustering=None):
        self.args = args
        self.views = {}
        self.states_table = None
        self.p_meta = {}
        self.title = 'Unknown Project'
        self.anvio_news = None

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.mode = A('mode')
        self.pan_db_path = A('pan_db')
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.genomes_storage_path = A('genomes_storage')
        self.collection_name = A('collection_name')
        self.manual_mode = A('manual_mode')
        self.split_hmm_layers = A('split_hmm_layers')
        self.show_all_layers = A('show_all_layers')
        self.taxonomic_level = A('taxonomic_level') or 't_genus'
        self.additional_layers_path = A('additional_layers')
        self.additional_view_path = A('additional_view')
        self.view = A('view')
        self.fasta_file = A('fasta_file')
        self.view_data_path = A('view_data')
        self.tree = A('tree')
        self.item_order_path = A('items_order')
        self.title = A('title')
        self.output_dir = A('output_dir')
        self.show_views = A('show_views')
        self.state_autoload = A('state_autoload')
        self.collection_autoload = A('collection_autoload')
        self.read_only = A('read_only')
        self.show_states = A('show_states')
        self.skip_check_names = A('skip_check_names')
        self.list_collections = A('list_collections')
        self.load_full_state = A('load_full_state')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.skip_init_functions = A('skip_init_functions')
        self.skip_auto_ordering = A('skip_auto_ordering')
        self.bin_ids_file_path = A('bin_ids_file')
        self.bin_id = A('bin_id')
        self.collection_name = A('collection_name')
        self.gene_mode = A('gene_mode')
        self.inspect_split_name = A('split_name')
        self.just_do_it = A('just_do_it')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.skip_news = A('skip_news')
        self.annotation_source_for_per_split_summary = A('annotation_source_for_per_split_summary')

        if self.pan_db_path and self.profile_db_path:
            raise ConfigError("You can't set both a profile database and a pan database in arguments "
                               "you send to this class. What are you doing?")

        if self.profile_db_path and filesnpaths.is_file_exists(self.profile_db_path, dont_raise=True):
            utils.is_profile_db(self.profile_db_path)

        if self.additional_layers_path:
            filesnpaths.is_file_tab_delimited(self.additional_layers_path)

        if self.gene_mode:
            if self.collection_name is None or self.bin_id is None:
                raise ConfigError("Gene view requires a collection and a bin to be specified. To SEE all collections and bins in "
                                  "your profile databse, you can use the program `anvi-show-collections-and-bins`. If you want to "
                                  "ADD a collection to your profile database, you can either use `anvi-import-collection` or "
                                  "`anvi-cluster-contigs` if appropriate, or you can go lazy and use the program "
                                  "`anvi-script-add-default-collection`, which would add a single collection that describes all contigs "
                                  "in your contigs database (in an ideal world you should use the last one only if you are working with "
                                  "a single genome).")

        if self.collection_name and (not self.gene_mode) and (self.bin_id or self.bin_ids_file_path) and self.mode != 'refine':
            raise ConfigError("There is something confusing here. On the one hand you provide a collection name, telling "
                              "anvi'o that you wish to run the interactive display with a collection focus. It would have "
                              "been fine if you stopped there, because then anvi'o would have fired up in 'collection "
                              "mode'. But you *also* provided a bin name. Now that's another ball game. Well, maybe you "
                              "wish to display this bin in 'gene mode'. But for that, you would need to add the flag "
                              "`--gene-mode`. Or maybe you actually want to 'refine' this bin, but in that case you "
                              "would need to run the program `anvi-refine` instead of `anvi-interactive` with the same "
                              "list of parameters and flags. Either of these would have made things much less confusing "
                              "here :(")

        if self.annotation_source_for_per_split_summary and self.collection_name:
            raise ConfigError("Unfortunately the flag `--annotation-source-for-per-split-summary` is not compatible with the collection "
                              "mode :/ This is because while this flag instructs anvi'o to prepare a split-level summary of functions, "
                              "in reality the splits in the collection mode are not actual splits, but 'bins'. Thus, the use of this "
                              "flag in collection mode requires additional logic to consolidate all information in splits in "
                              "a single bin to 're-summarize' data to the level of individual bins. This logic is currently not in place, "
                              "hence the flag is useless in this mode :( Apologies for the inconvenience. Push us through GitHub if "
                              "you really want us to address this.")

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        self.displayed_item_names_ordered = None
        self.auxiliary_profile_data_available = False

        # get additional data for items and layers, and get layer orders data.
        try:
            a_db_is_found = (os.path.exists(self.pan_db_path) if self.pan_db_path else False) or (os.path.exists(self.profile_db_path) if self.profile_db_path else False)
            self.items_additional_data_keys, self.items_additional_data_dict = TableForItemAdditionalData(self.args).get() if a_db_is_found else ([], {})
            self.layers_additional_data_keys, self.layers_additional_data_dict = TableForLayerAdditionalData(self.args).get_all() if a_db_is_found else ([], {})
            self.layers_order_data_dict = TableForLayerOrders(self.args).get() if a_db_is_found else {}
        except GenesDBError as e:
            self.items_additional_data_keys, self.items_additional_data_dict = [], {}
            self.layers_additional_data_keys, self.layers_additional_data_dict = [], {}
            self.layers_order_data_dict = {}

            run.warning("Most likely the misc data module is complaining about a missing genes database. If that's the case,\
                         you can ignore it as this will sort itself out in a second (this is a workaround for a design\
                         bottleneck sadface.png): %s" % e.clear_text(), header="EXCEPTION OMMITTED (BUT PROBABLY YOU'RE FINE)", lc='yellow')

        for group_name in self.layers_additional_data_keys:
            layer_orders = TableForLayerOrders(self.args).update_orders_dict_using_additional_data_dict({},
                self.layers_additional_data_keys[group_name], self.layers_additional_data_dict[group_name]) if a_db_is_found else {}
            for order_name in layer_orders:
                self.layers_order_data_dict['%s :: %s' % (group_name, order_name)] = layer_orders[order_name]

        # make sure the mode will be set properly
        if self.collection_name and self.manual_mode:
            raise ConfigError("You can't anvi-interactive in manual mode with a collection name.")

        self.external_clustering = external_clustering

        self.collections = ccollections.Collections()

        # we will set the split names of interest early on to avoid
        # reading everything in the contigs database.
        self.split_names_of_interest = set([])
        self.bin_names_of_interest = []
        if self.collection_name and (self.bin_ids_file_path or self.bin_id):
            progress.new('Initializing split names of interest')
            progress.update('...')

            splits_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
            self.bin_names_of_interest = list(splits_dict.keys())

            for split_names in list(splits_dict.values()):
                self.split_names_of_interest.update(split_names)

            progress.end()
        elif self.inspect_split_name:
            self.split_names_of_interest = set([self.inspect_split_name])
            self.skip_check_names = True

        if self.contigs_db_path:
            self.contigs_db_variant = utils.get_db_variant(self.contigs_db_path)
            self.completeness = Completeness(self.contigs_db_path)
            self.collections.populate_collections_dict(self.contigs_db_path)
        else:
            self.contigs_db_variant = None
            self.completeness = None

        # if the mode has not been set from within the arguments, we will set something up here:
        if not self.mode:
            if self.manual_mode:
                self.mode = 'manual'
            elif self.gene_mode:
                # collection mode and gene view mode both uses collection_name
                # so gene_mode needs to be placed before collection view
                self.mode = 'gene'
            elif self.collection_name or self.list_collections:
                self.mode = 'collection'
            elif self.contigs_db_variant == 'trnaseq':
                self.mode = 'trnaseq'
            else:
                self.mode = 'full'

        if self.mode in ['full', 'collection', 'trnaseq', 'gene'] and not self.profile_db_path:
            raise ConfigError("You must declare a profile database for this to work :(")

        ContigsSuperclass.__init__(self, self.args)
        self.init_splits_taxonomy(self.taxonomic_level)

        # make sure we are not dealing with apples and oranges here.
        if self.contigs_db_path and self.profile_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.P = lambda x: os.path.join(self.p_meta['output_dir'], x)
        self.cwd = os.getcwd()

        # here is where the big deal stuff takes place. depending on the mode, we will call
        # the appropriate function for initializing the interafce class.
        self.run.info('Interactive mode', self.mode, mc='green')
        if self.mode == 'manual':
            self.load_manual_mode()
        elif self.mode == 'gene':
            self.load_gene_mode()
        elif self.mode == 'refine':
            self.load_full_mode()
            self.load_refine_mode()
        elif self.mode == 'trnaseq':
            self.load_full_mode()
            self.load_trnaseq_mode()
        elif self.mode == 'functional':
            self.load_functional_mode()
        elif self.mode == 'pan':
            self.load_pan_mode()
        elif self.mode == 'collection':
            self.load_collection_mode()
        elif self.mode == 'full':
            self.load_full_mode()
        elif self.mode == 'inspect':
            self.load_full_mode()
            self.load_inspect_mode()
        else:
            raise ConfigError("The interactive class is called with a mode that no one knows anything "
                              "about. '%s'... What kind of a mode is that anyway :/" % self.mode)

        if self.external_clustering:
            self.p_meta['clusterings'] = self.clusterings = self.external_clustering['clusterings']
            self.p_meta['available_clusterings'] = list(self.clusterings.keys())
            self.p_meta['default_clustering'] = self.external_clustering['default_clustering']

        if not self.state_autoload and 'default' in self.states_table.states:
            self.state_autoload = 'default'

        if not self.collection_autoload and 'default' in self.collections.collections_dict:
            self.collection_autoload = 'default'

        self.check_for_clusterings()
        self.set_displayed_item_names()

        # now we know what splits we are interested in (self.displayed_item_names_ordered), we can get rid of all the
        # unnecessary splits stored in views and other additional dicts.
        self.prune_view_and_additional_data_dicts()

        self.process_external_item_order()
        self.gen_item_auto_order_length()
        self.gen_item_auto_order_alphabetical()

        if not self.p_meta['default_item_order'] and len(self.p_meta['available_item_orders']):
            self.p_meta['default_item_order'] = self.p_meta['available_item_orders'][0]

        # if there are any HMM search results in the contigs database other than 'singlecopy' sources,
        # we would like to visualize them as additional layers. following function is inherited from
        # Contigs DB superclass and will fill self.hmm_searches_dict if appropriate data is found in
        # search tables:
        if self.mode == 'full' or self.mode == 'refine':
            self.init_non_singlecopy_gene_hmm_sources(self.displayed_item_names_ordered, return_each_gene_as_a_layer=self.split_hmm_layers)

        # take care of additional layers, and update ordering information for items
        if self.additional_layers_path:
            self.items_additional_data_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_layers_path, dict_to_append=self.items_additional_data_dict, assign_none_for_missing=True)
            self.items_additional_data_keys = self.items_additional_data_keys + utils.get_columns_of_TAB_delim_file(self.additional_layers_path)

        self.check_names_consistency()
        self.gen_orders_for_items_based_on_additional_layers_data()
        self.finalize_view_data()
        self.get_anvio_news()


    def set_displayed_item_names(self):
        """Sets the master list of names .. UNLESS we are in manual or functional mode, in which case names
           will be set within corresponding functions"""

        if self.mode == 'manual' or self.mode == 'functional':
            return

        if not self.p_meta['item_orders']:
            raise ConfigError("Apologies :( The code managed to came all the way down here without any item "
                              "orders set for your project, which makes it impossible for "
                              "anvi'o intearctive interface to display anything :/ Please check warning "
                              "messages you may have on your screen.")

        # self.displayed_item_names_ordered is going to be the 'master' names list. everything else is going to
        # need to match these names:
        if self.p_meta['default_item_order']:
            default_item_order = self.p_meta['item_orders'][self.p_meta['default_item_order']]
        else:
            default_item_order = self.p_meta['item_orders'][self.p_meta['available_item_orders'][0]]

        if default_item_order['type'] == 'newick':
            # anvi'o complains about newick trees with items names of which are composed of
            # integers. but in the case of gene mode, all items will have integer names as
            # gene caller ids in anvi'o are always integers
            names_with_only_digits_ok = self.mode == 'gene'

            self.displayed_item_names_ordered = utils.get_names_order_from_newick_tree(default_item_order['data'], reverse=True, names_with_only_digits_ok=names_with_only_digits_ok)
        elif default_item_order['type'] == 'basic':
            self.displayed_item_names_ordered = default_item_order['data']
        else:
            raise ConfigError("There is something wrong here, and anvi'o needs and adult :( Something that should "
                              "never happen happened. The default clustering does not have a basic or newick type.")


    def check_for_clusterings(self):
        if self.mode == 'manual':
            return

        if not self.p_meta['item_orders'] or not len([o for o in self.p_meta['item_orders'].values() if o['type'] == 'newick']):
            if self.p_meta['db_type'] == 'pan':
                self.run.warning("This pangenome (which you gracefully named as '%s') does not seem to have any hierarchical "
                                 "clustering of gene clusters it contains. Maybe you skipped the clustering step, maybe "
                                 "anvi'o skipped it on your behalf because you had too many gene clusters or something. Anvi'o "
                                 "will do its best to recover from this situation. But you will very likely not be able to see "
                                 "a hierarchical dendrogram in the resulting display even if everything goes alright.. In some "
                                 "cases using a parameter like `--min-occurrence 2`, which would reduce the number of gene clusters by "
                                 "removing singletons that appear in only one genome can help solve this issue and make your "
                                 "pangenome great again." % (self.p_meta['project_name']))
            else:
                if self.item_order_path:
                    self.run.warning("This merged profile database does not seem to have any hierarchical clustering "
                                     "of splits that is required by the interactive interface. But it seems you did provide "
                                     "an items order file. So anvi'o will try to use that and display your data.")
                else:
                    if self.p_meta['merged']:
                        self.run.warning("This merged profile database does not seem to have any hierarchical clustering "
                                         "of splits :/ It may be the case becasue it may have been generated by anvi-merge "
                                         "with the `--skip-hierarchical-clustering` flag, or hierarchical "
                                         "clustering step may have been skipped by anvi-merge because you had too many splits "
                                         "to get the clustering in a reasonable amount of time. Basically you will not see a "
                                         "clustering dendrogram in the center of your display. Please read the help menu for "
                                         "anvi-merge, and/or refer to the tutorial: "
                                         "http://merenlab.org/2015/05/01/anvio-tutorial/#clustering-during-merging")
                    else:
                        self.run.warning("This single profile database does not seem to have any hierarchical clustering "
                                         "that is required by the interactive interface :/ Anvi'o does not automatically compute "
                                         "a hierarchical clustering of contigs for single profiles, but you can ask for it "
                                         "explicilty by passing the `--cluster-contigs` flag to the program `anvi-profile`. "
                                         "This way you will be able to run `anvi-interactive` on your single profile, and see "
                                         "a dendrogram at the center of your display. Please take a look at the help menu for "
                                         "`anvi-profile`.")


    def gen_orders_for_items_based_on_additional_layers_data(self):
        if self.skip_auto_ordering or not self.items_additional_data_keys:
            return

        self.progress.new('Processing additional data to order items (to skip: --skip-auto-ordering)')
        skipped_additional_data_layers = []

        sum_stackbar_items = {}
        for layer in [additional_layer for additional_layer in self.items_additional_data_keys]:
            if '!' in layer:
                stackbar_name = layer.split('!')[0]
                if stackbar_name not in sum_stackbar_items:
                    sum_stackbar_items[stackbar_name] = {}

                for item in self.displayed_item_names_ordered:
                    if item not in sum_stackbar_items[stackbar_name]:
                        sum_stackbar_items[stackbar_name][item] = 0.0

                    if item in self.items_additional_data_dict:
                        sum_stackbar_items[stackbar_name][item] += float(self.items_additional_data_dict[item][layer])

        for layer in [additional_layer for additional_layer in self.items_additional_data_keys]:
            self.progress.update('for "%s" ...' % layer)
            layer_type = utils.get_predicted_type_of_items_in_a_dict(self.items_additional_data_dict, layer)

            if layer_type == None:
                skipped_additional_data_layers.append(layer)
                continue

            item_layer_data_tuple = []
            items_for_which_we_put_zeros_for_missing_values = set([])
            for item in self.displayed_item_names_ordered:
                if item not in self.items_additional_data_dict:
                    if layer_type != str:
                        item_layer_data_tuple.append((0.0, item))
                        items_for_which_we_put_zeros_for_missing_values.add(item)
                    else:
                        item_layer_data_tuple.append(('', item))
                else:
                    if self.items_additional_data_dict[item][layer] == None:
                        if layer_type != str:
                            items_for_which_we_put_zeros_for_missing_values.add(item)
                            item_layer_data_tuple.append((0.0, item))
                        else:
                            item_layer_data_tuple.append(('', item))
                    else:
                        if '!' in layer:
                            stackbar_name = layer.split('!')[0]
                            if float(sum_stackbar_items[stackbar_name][item]) != 0.0:
                                item_layer_data_tuple.append((float(self.items_additional_data_dict[item][layer]) / (1.0 * float(sum_stackbar_items[stackbar_name][item])), item))
                            else:
                                item_layer_data_tuple.append((0.0, item))
                        else:
                            item_layer_data_tuple.append((layer_type(self.items_additional_data_dict[item][layer]), item))

            if len(items_for_which_we_put_zeros_for_missing_values):
                self.progress.end()
                self.run.warning("OK. While working on the layer '%s', which actually looked like a numerical layer, anvi'o realized "
                                 "that %d of your items (for instance '%s' was one of them) did not have a value for this layer. To "
                                 "make sure things will continue working in the interface, anvi'o took the liberty of adding zeros "
                                 "as values for these items. Which is not the smartest thing to do, but we unfortunately do not "
                                 "support empty values for numerical layers yet. In MetalBeard's voice: things shall continue to "
                                 "work, but ye here be warned. Back to anvi'o regular voice: Please keep this in mind while you "
                                 "are studying the interactive interface be extra careful how to interpret your analysis when "
                                 "you see zero values in the layer '%s'." % (layer,
                                                                             len(items_for_which_we_put_zeros_for_missing_values),
                                                                             items_for_which_we_put_zeros_for_missing_values.pop(),
                                                                             layer))
                self.progress.new('Processing additional data to order items (to skip: --skip-auto-ordering)')


            # try to fancify the layer names that will appear in items order combo for stacked
            # bar data type before adding them in:
            if '!' in layer:
                stacked_bar_name, item_name = layer.split('!')
                layer_name = '%s [%s]' % (stacked_bar_name, item_name)
            else:
                layer_name = layer

            self.p_meta['available_item_orders'].append('>> %s:none:none' % layer_name)
            self.p_meta['item_orders']['>> %s' % layer_name] = {'type': 'basic', 'data': [i[1] for i in sorted(item_layer_data_tuple)]}
            self.p_meta['available_item_orders'].append('>> %s_(reverse):none:none' % layer_name)
            self.p_meta['item_orders']['>> %s_(reverse)' % layer_name] = {'type': 'basic', 'data': [i[1] for i in sorted(item_layer_data_tuple, reverse=True)]}

        self.progress.end()

        if len(skipped_additional_data_layers):
            self.run.warning("One or more of your additional data columns were completely empty. Like, they didn't have any data at all :/ "
                             "In the best case scenario you will see completely blank layers in your display. In the worst case scenario "
                             "other things will break. Since you are a curious person, anvi'o thought you would like to know. These are "
                             "the empty variables: %s." % ', '.join(['"%s"' % s for s in skipped_additional_data_layers]))


    def gen_item_auto_order_length(self):
        """Generate an items order based on item length when a contigs-db is present.

        Since anvi'o soft-splits long contigs, each split that belongs to the same parent
        contig will have a different length. Ordering solely based on item lengths will
        break that linkage, thus, this function respets the 'parent' information and
        first orders contigs based on length, and then orders splits with respect to the
        length of contigs to which they belong.
        """

        if self.skip_auto_ordering or not self.contigs_db_path:
            return

        self.progress.new('Making items order: Length')
        self.progress.update('...')

        # get a dictionary to translate between contig names and split names
        contig_name_to_splits_dict = utils.get_contig_name_to_splits_dict(self.contigs_db_path)

        # get a list of contig names based on their lengths
        contig_names_sorted_by_length = [tpl[0] \
                for tpl in sorted([(k, self.contigs_basic_info[k]['length']) \
                for k in self.contigs_basic_info], key = lambda x: x[1], reverse=True)]

        split_names = set(self.displayed_item_names_ordered)

        split_names_ordered_by_size = []
        for contig_name in contig_names_sorted_by_length:
            for split_name in sorted(contig_name_to_splits_dict[contig_name]):
                if split_name in split_names:
                    split_names_ordered_by_size.append(split_name)

        self.p_meta['item_orders']['<> Length:none:none'] = {'type': 'basic', 'data': split_names_ordered_by_size[::-1]}
        self.p_meta['available_item_orders'].append('<> Length:none:none')

        self.progress.end()


    def gen_item_auto_order_alphabetical(self):
        """This function populates self.p_meta with additional organizations of data, such as alphabetical ordering\
           of data items, etc. In the interface these additional orders appear in the 'items order' combo box"""

        if self.skip_auto_ordering:
            return

        self.progress.new('Making items order: Alphabetical')
        self.progress.update('...')

        # add an alphabetical order:
        self.p_meta['item_orders']['<> Alphabetical:none:none'] = {'type': 'basic', 'data': sorted(self.displayed_item_names_ordered[::-1], reverse=True)}
        self.p_meta['available_item_orders'].append('<> Alphabetical:none:none')

        # and the reverse-alphabetical, too:
        self.p_meta['item_orders']['<> Alphabetical_(reverse):none:none'] = {'type': 'basic', 'data': sorted(self.displayed_item_names_ordered)}
        self.p_meta['available_item_orders'].append('<> Alphabetical_(reverse):none:none')

        self.progress.end()


    def process_external_item_order(self):
        """This function processes self.item_order_path to update available item orders"""

        if not self.item_order_path:
            return

        filesnpaths.is_file_exists(self.item_order_path)

        item_order = [l.strip() for l in open(self.item_order_path, 'r').readlines()]
        self.run.info('Items order', 'An items order with %d items is found at %s.' % (len(item_order), self.item_order_path), mc='cyan')

        self.progress.new('External items order')
        self.progress.update('...')

        # make sure all items we will be working with is in items order.
        item_order_set = set(item_order)
        missing_items = [i for i in self.displayed_item_names_ordered if i not in item_order_set]
        if (missing_items):
            self.progress.end()
            raise ConfigError("While processing your items order file, anvi'o realized that some of the items in your view data are not "
                              "in your items order file. In fact there are like %d of them missing, and one of the missing items look "
                              "like this if it makes any sense: '%s'" % (len(missing_items), missing_items.pop()))

        if len(item_order) != len(self.displayed_item_names_ordered):
            self.progress.end()
            raise ConfigError("While processing your items order file, anvi'o realized that the number of items described in your file "
                              "(%s) is not equal to the number of items you have in your view data (%s). This is totally a deal "
                              "breaker :/" % (pp(len(item_order)), pp(len(self.displayed_item_names_ordered))))

        # because of the special case of items order, the clusterings and available_clusterings items
        # may not be initialized. check them first.
        if not self.p_meta['item_orders']:
            self.p_meta['item_orders'] = {}
        if not self.p_meta['available_item_orders']:
            self.p_meta['available_item_orders'] = []

        # add an alphabetical order:
        self.p_meta['item_orders']['<> User order:none:none'] = { 'type': 'basic', 'data': item_order[::-1]}
        self.p_meta['available_item_orders'].append('<> User order:none:none')

        # and the reverse-alphabetical, too:
        self.p_meta['item_orders']['<> User order (reverse):none:none'] = {'type': 'basic', 'data': item_order}
        self.p_meta['available_item_orders'].append('<> User order (reverse):none:none')

        self.progress.end()


    def load_manual_mode(self):
        if self.contigs_db_path:
            raise ConfigError("When you want to use the interactive interface in manual mode, you must "
                               "not use a contigs database.")

        if not self.profile_db_path:
            raise ConfigError("Even when you want to use the interactive interface in manual mode, you need "
                               "to provide a profile database path. But you DO NOT need an already existing "
                               "profile database, since anvi'o will generate an empty one for you. The profile "
                               "database in this mode only used to read or store the 'state' of the display "
                               "for visualization purposes, or to allow you to create and store collections.")

        # we don't expect ad hoc profile databases to have views tables,
        # but if the profile db has one, we don't want to ignore that either.
        # we will define this variable here as none, and if there is a profile
        # db path, we will fill it in if it contains view tables.
        views_table = None

        # if the user is using an existing profile database, we need to make sure that it is not associated
        # with a contigs database, since it would mean that it is a full anvi'o profile database and should
        # not be included in manual operations.
        tree_order_found_in_db = False
        item_orders_in_db = None
        if filesnpaths.is_file_exists(self.profile_db_path, dont_raise=True):
            profile_db = ProfileDatabase(self.profile_db_path)
            if profile_db.meta['contigs_db_hash']:
                raise ConfigError("Well. It seems the profile database is associated with a contigs database,\
                                    which means using it in manual mode is not the best way to use it. Probably\
                                    what you wanted to do is to let the manual mode create a new profile database\
                                    for you. Simply type in a new profile database path (it can be a file name\
                                    that doesn't exist).")

            # if there is a profile database, let's check whether there are any item orders available in it.
            # we will later put this in p_meta down below:
            item_orders_in_db = profile_db.db.get_table_as_dict(t.item_orders_table_name)
            for item_order in item_orders_in_db:
                if item_orders_in_db[item_order]['type'] == 'basic':
                    try:
                        item_orders_in_db[item_order]['data'] = item_orders_in_db[item_order]['data'].split(',')
                    except:
                        self.progress.end()
                        raise ConfigError("Something is wrong with the basic order `%s` in this profile database :(" % (item_order))
                elif item_orders_in_db[item_order]['type'] == 'newick':
                    tree_order_found_in_db = item_order

            if t.views_table_name in profile_db.db.get_table_names():
                # surprise surprise!
                views_table = profile_db.db.get_table_as_dict(t.views_table_name)

            profile_db.disconnect()


        if not self.tree and not self.view_data_path and not tree_order_found_in_db:
            raise ConfigError("You must be joking Mr. Feynman. No tree file, and no data file, and no tree order in the "
                              "database? What is it that anvi'o supposed to visualize? :(")

        if not self.tree and not tree_order_found_in_db:
            self.run.warning("You haven't declared a tree file and there was no tree order in the database. Anvi'o will "
                             "do its best to come up with an organization of your items.")

        if self.view:
            raise ConfigError("You can't use '--view' parameter when you are running the interactive interface "
                              "in manual mode")

        if self.show_views:
            raise ConfigError("Sorry, there are no views to show in manual mode :/")

        if self.show_states:
            raise ConfigError("Sorry, there are no states to show in manual mode :/")

        if self.tree:
            filesnpaths.is_file_exists(self.tree)
            newick_tree_text = ''.join([l.strip() for l in open(os.path.abspath(self.tree)).readlines()])
            self.displayed_item_names_ordered = sorted(utils.get_names_order_from_newick_tree(newick_tree_text))
        elif tree_order_found_in_db:
            self.displayed_item_names_ordered = sorted(utils.get_names_order_from_newick_tree(item_orders_in_db[tree_order_found_in_db]['data']))
        else:
            self.displayed_item_names_ordered = sorted(utils.get_column_data_from_TAB_delim_file(self.view_data_path, column_indices=[0])[0][1:])

        # try to convert item names into integer values for proper sorting later. it's OK if it does
        # not work.
        try:
            self.displayed_item_names_ordered = [str(s) for s in sorted([int(n) for n in self.displayed_item_names_ordered])]
        except:
            pass

        view_data_path = os.path.abspath(self.view_data_path) if self.view_data_path else None
        self.p_meta['splits_fasta'] = os.path.abspath(self.fasta_file) if self.fasta_file else None
        self.p_meta['output_dir'] = None
        self.p_meta['views'] = {}
        self.p_meta['db_type'] = 'profile'
        self.p_meta['merged'] = True
        self.p_meta['blank'] = True

        # set some default organizations of data:
        if not item_orders_in_db:
            self.p_meta['item_orders'] = {}
            self.p_meta['available_item_orders'] = []
            self.p_meta['default_item_order'] = []
        else:
            self.p_meta['item_orders'] = item_orders_in_db
            self.p_meta['available_item_orders'] = sorted(list(item_orders_in_db.keys()))
            self.p_meta['default_item_order'] = tree_order_found_in_db or self.p_meta['available_item_orders'][0]

        # if we have a tree, let's make arrangements for it:
        if self.tree:
            item_order_name = '%s:unknown:unknown' % filesnpaths.get_name_from_file_path(self.tree)
            self.p_meta['default_item_order'] = item_order_name
            self.p_meta['available_item_orders'].append(item_order_name)
            self.p_meta['item_orders'][item_order_name] = {'type': 'newick', 'data': newick_tree_text}


        if views_table and filesnpaths.is_file_exists(self.profile_db_path, dont_raise=True):
            # we have view tables
            profile_db = ProfileDatabase(self.profile_db_path)

            for view in views_table:
                table_name = views_table[view]['target_table']

                data, header = profile_db.db.get_view_data(table_name)

                self.views[view] = {'table_name': table_name,
                                    'header': header,
                                    'dict': data}

            self.p_meta['default_view'] = sorted(list(self.views.keys()), reverse=True)[0]
            self.default_view = self.p_meta['default_view']
            profile_db.disconnect()

        elif self.view_data_path:
            # sanity of the view data
            filesnpaths.is_file_tab_delimited(view_data_path)
            view_data_columns = utils.get_columns_of_TAB_delim_file(view_data_path, include_first_column=True)

            utils.check_misc_data_keys_for_format(view_data_columns)

            # load view data as the default view:
            self.p_meta['default_view'] = 'single'
            self.default_view = self.p_meta['default_view']
            self.views[self.default_view] = {'header': view_data_columns[1:],
                                             'dict': utils.get_TAB_delimited_file_as_dictionary(view_data_path)}
        else:
            # no view data is provided... it is only the tree we have. we will creaet a mock 'view data dict'
            # here using what is in the tree.
            ad_hoc_dict = {}
            for item in self.displayed_item_names_ordered:
                ad_hoc_dict[item] = {'names': str(item)}

            self.p_meta['default_view'] = 'single'
            self.default_view = self.p_meta['default_view']
            self.views[self.default_view] = {'header': ['names'],
                                             'dict': ad_hoc_dict}

        # we assume that the sample names are the header of the view data, so we might as well set it up:
        sample_names = [self.title.replace(' ', '_')] if self.title else self.views[self.default_view]['header']
        self.p_meta['samples'] = sample_names
        self.p_meta['sample_id'] = 'AD HOC DISPLAY'

        # if we have an input FASTA file, we will set up the split_sequences and splits_basic_info dicts,
        # otherwise we will leave them empty
        self.splits_basic_info = {}
        self.split_sequences = None
        if self.p_meta['splits_fasta']:
            filesnpaths.is_file_fasta_formatted(self.p_meta['splits_fasta'])

            # the utils.get_FASTA_file_as_dictionary returns a dictionary, but anvi'o expects
            # a different format when bottleroutes accesses interactive.split_sequences. so
            # here we will first turn the naive dictionary from get_FASTA_file_as_dictionary
            # into a dictionary with 'sequence' keys:
            self.split_sequences = {}
            split_sequences = utils.get_FASTA_file_as_dictionary(self.p_meta['splits_fasta'])
            for split_name in split_sequences:
                self.split_sequences[split_name] = {'sequence': split_sequences[split_name]}

            names_missing_in_FASTA = set(self.displayed_item_names_ordered) - set(self.split_sequences.keys())
            num_names_missing_in_FASTA = len(names_missing_in_FASTA)
            if num_names_missing_in_FASTA:
                raise ConfigError('Some of the names in your view data does not have corresponding entries in the '
                                   'FASTA file you provided. Here is an example to one of those %d names that occur '
                                   'in your data file, but not in the FASTA file: "%s"' % (num_names_missing_in_FASTA, names_missing_in_FASTA.pop()))

            # setup a mock splits_basic_info dict
            for split_id in self.displayed_item_names_ordered:
                self.splits_basic_info[split_id] = {'length': len(self.split_sequences[split_id]['sequence']),
                                                    'gc_content': utils.get_GC_content_for_sequence(self.split_sequences[split_id]['sequence'])}

        # create a new, empty profile database for manual operations
        if not os.path.exists(self.profile_db_path):
            profile_db = ProfileDatabase(self.profile_db_path)
            profile_db.create({'db_type': 'profile',
                               'db_variant': 'ad-hoc-display',
                               'blank': True,
                               'merged': True,
                               'contigs_db_hash': None,
                               'items_ordered': False,
                               'samples': ', '.join(self.p_meta['samples']),
                               'sample_id': self.p_meta['sample_id']})

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path)

        # also populate collections, if there are any
        self.collections.populate_collections_dict(self.profile_db_path)

        # read description from self table, if it is not available get_description function will return placeholder text
        self.p_meta['description'] = get_description_in_db(self.profile_db_path)

        self.title = self.title or self.p_meta['sample_id']

        # did the user ask anvi'o to add an additional layer to summarize proportion of functions per split?
        self.update_items_additional_data_with_functions_per_split_summary()


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

            clusterings[clustering_id] = {'type': 'newick', 'data': newick}

        run.info('Available clusterings', list(clusterings.keys()))

        return clusterings


    def load_inspect_mode(self):
        self.displayed_item_names_ordered = sorted(self.views[self.default_view]['dict'].keys())
        alphabetical_order = {'type': 'basic', 'data': sorted(self.displayed_item_names_ordered[::])}
        self.p_meta['item_orders']['alphabetical'] = alphabetical_order
        self.p_meta['available_item_orders'].append('alphabetical')

        if not self.inspect_split_name or self.inspect_split_name not in self.displayed_item_names_ordered:
            if self.just_do_it:
                self.inspect_split_name = alphabetical_order['data'][0]
                self.run.warning("Since you have asked so kindly, anvi'o decided to start the interactive "
                                 "interface with this split: %s" % self.inspect_split_name)
            else:
                raise ConfigError("Either you forgot to provide a split name to `anvi-inspect` or the split name "
                                  "you have provided does not exist. If you don't care and want to start the "
                                  "interactive inteface with a random split from the profile database, please use "
                                  "the flag `--just-do-it`")


    def load_functional_mode(self):
        # these are all going to be filled later nicely.
        self.p_meta['item_orders'] = {}
        self.p_meta['available_item_orders'] = []
        self.p_meta['default_item_order'] = 'presence_absence'

        # these are irrelevant to this mode, but necessary to fill in.
        self.p_meta['splits_fasta'] = None
        self.p_meta['output_dir'] = None
        self.p_meta['views'] = {}
        self.p_meta['db_type'] = 'profile'
        self.p_meta['merged'] = True
        self.p_meta['blank'] = True
        self.splits_basic_info = {}
        self.split_sequences = None

        if not self.profile_db_path:
            raise ConfigError("It is OK if you don't have a profile database, but you still should provide a "
                              "profile database path so anvi'o can generate a new one for you. The profile "
                              "database in this mode only used to read or store the 'state' of the display "
                              "for visualization purposes, or to allow you to create and store collections. "
                              "Still very useful.")

        if os.path.exists(self.profile_db_path):
            raise ConfigError(f"So. Everytime you run the anvi'o interactive interface in display functions mode "
                              f"you must provide a profile database that actually does not exist. Yeah, we know "
                              f"it doesn't make sense, but actually it does, but really very difficult to "
                              f"explain. If you really want to re-display your existing profile database with "
                              f"functions, you should run `anvi-interactive -p {self.profile_db_path} --manual`. "
                              f"Now, that should work like a charm.")

        # initialize all functions from wherever.
        facc = AggregateFunctions(self.args, r=self.run, p=self.progress)

        num_facc_keys = len(facc.functions_across_layers_frequency)
        num_genomes = len(facc.layer_names_considered)
        if num_facc_keys < 25:
            self.run.warning(f"Among all of your {num_genomes} genomes, anvi'o found only {num_facc_keys} unique "
                             f"{facc.K()}s for {facc.function_annotation_source}. This doesn't look very good, but "
                             f"it is difficult to know whether you find this normal or even pleasing... Things will "
                             f"continue to work, but we just thought you should know about this.")
        if num_facc_keys > 2500 and num_facc_keys < 10000:
            self.run.warning(f"Among all of your {num_genomes} genomes, anvi'o found as many as {num_facc_keys} unique "
                             f"{facc.K()}s for {facc.function_annotation_source}. It looks like you have quite a lot of "
                             f"them, and this may cause some delays during the calculation of your dendrogram. Anvi'o wiil "
                             f"do her best, but if things crash downstream, this will likely be the reason.")
        if num_facc_keys > 10000:
            raise ConfigError(f"Among all of your {num_genomes} genomes, anvi'o found as many as {num_facc_keys} unique "
                              f"{facc.K()}s for {facc.function_annotation_source} :/ This is just way too much for anvi'o to try to "
                              f"put together a hierarchical clustering for and display. If you think this had to work, let us know and "
                              f"we will see if we can do something to help you.")

        run.warning(None, header="SUMMARY OF WHAT IS GOING ON", lc="green")
        self.run.info('Num genomes', num_genomes)
        self.run.info('Function annotation source', facc.function_annotation_source)
        self.run.info('Num unique keys', num_facc_keys)
        self.run.info('Keys correspond to', f"'{facc.K()}s' rather than '{facc.V()}s'")
        self.run.info('Only the best hits are considered', "False" if facc.aggregate_using_all_hits else "True", nl_after=1)

        # now we will work on views and clustering our data to get newick trees.
        self.p_meta['default_view'] = 'functions_presence_absence_view'
        self.default_view = self.p_meta['default_view']

        # our 'samples' in this context are individual genomes we are about to display
        self.p_meta['samples'] = list(facc.layer_names_considered)
        self.p_meta['sample_id'] = f"{facc.function_annotation_source} DISPLAY"

        # setup the views dict
        self.views = {'functions_frequency_view'       : {'header': sorted(list(facc.layer_names_considered)),
                                                          'dict': facc.functions_across_layers_frequency},
                      'functions_presence_absence_view': {'header': sorted(list(facc.layer_names_considered)),
                                                          'dict': facc.functions_across_layers_presence_absence}}

        # create a new, empty profile database for manual operations
        if not os.path.exists(self.profile_db_path):
            profile_db = ProfileDatabase(self.profile_db_path)
            profile_db.create({'db_type': 'profile',
                               'blank': True,
                               'merged': True,
                               'db_variant': 'functions-display',
                               'description': facc.summary_markdown,
                               'function_annotation_source': facc.function_annotation_source,
                               'min_occurrence_of_function': facc.min_occurrence,
                               'aggregate_based_on_accession': facc.aggregate_based_on_accession,
                               'aggregate_using_all_hits': facc.aggregate_using_all_hits,
                               'contigs_db_hash': None,
                               'items_ordered': False,
                               'samples': ', '.join(self.p_meta['samples']),
                               'sample_id': self.p_meta['sample_id']})

        # now we start working on our views and item orders.
        for view in self.views:
            # first, generate an items order for a given view:
            items_order = clustering.get_newick_tree_data_for_dict(self.views[view]['dict'], zero_fill_missing=True, distance=self.distance, linkage=self.linkage)
            item_order_name = f"{view[10:-5]}"
            self.p_meta['available_item_orders'].append(item_order_name)
            self.p_meta['item_orders'][item_order_name] = {'type': 'newick', 'data': copy.deepcopy(items_order)}

            # then add the items order to the database
            dbops.add_items_order_to_db(self.profile_db_path, item_order_name, items_order, order_data_type_newick=True,
                                        distance=self.distance, linkage=self.linkage, make_default=True if view == 'functions_presence_absence_view' else False,
                                        check_names_consistency=False)

            # if there are more than one genome in the analysis (which should almost always be the case),
            # let's get a layer's order, too, and immediately add it to the database:
            if num_genomes > 1:
                layer_order = clustering.get_newick_tree_data_for_dict(self.views[view]['dict'], transpose=True, zero_fill_missing=True, distance=self.distance, linkage=self.linkage)
                args = argparse.Namespace(profile_db=self.profile_db_path, target_data_table="layer_orders", just_do_it=True)
                TableForLayerOrders(args, r=terminal.Run(verbose=False)).add({f"{view.upper()}": {'data_type': 'newick', 'data_value': layer_order}}, skip_check_names=True)
                self.layers_order_data_dict = TableForLayerOrders(args, r=terminal.Run(verbose=False)).get()

            # add vew tables to the database
            TablesForViews(self.profile_db_path).create_new_view(
                                            view_data=self.views[view]['dict'],
                                            table_name=f"{view}",
                                            view_name=f"{view}",
                                            from_matrix_form=True)

        # let's do this here as well so our dicts are not pruned.
        self.displayed_item_names_ordered = sorted(utils.get_names_order_from_newick_tree(items_order))

        # here we will store a layer for function names in the additional data tables of the profile databse and
        args = argparse.Namespace(profile_db=self.profile_db_path, target_data_table="items", just_do_it=True)
        TableForItemAdditionalData(args, r=terminal.Run(verbose=False)).add(facc.hash_to_function_dict, [facc.function_annotation_source], skip_check_names=True)

        # if we have functional enrichment analysis results for these genomes, let's add that into
        # the database as well!
        if facc.functional_enrichment_stats_dict:
            TableForItemAdditionalData(args, r=terminal.Run(verbose=False)).add(facc.functional_enrichment_stats_dict, ['enrichment_score', 'unadjusted_p_value', 'adjusted_q_value', 'associated_groups'], skip_check_names=True)

        # here we will read the items additional data back so it is both in there for future `anvi-interactive` ops,
        # AND in here in the interactive class to visualize the information.
        self.items_additional_data_keys, self.items_additional_data_dict = TableForItemAdditionalData(args, r=terminal.Run(verbose=False)).get()

        # everything we need is in the database now. time to add a mini state (note that we
        # replace the function layer name template with the annotation source on the fly):
        mini_state = open(os.path.join(os.path.dirname(anvio.__file__), 'data/mini-states/display-functions.json')).read().replace('__FUNCTIONS_LAYER_NAME__', facc.function_annotation_source)
        TablesForStates(self.profile_db_path).store_state('default', mini_state)

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path)

        # also populate collections, if there are any
        self.collections.populate_collections_dict(self.profile_db_path)

        # read description from self table, if it is not available get_description function will return placeholder text
        self.p_meta['description'] = get_description_in_db(self.profile_db_path)

        self.title = self.args.title or self.p_meta['sample_id']
        if facc.min_occurrence > 1:
            self.title += f" (MIN OCCURRENCE: {facc.min_occurrence})"


    def load_refine_mode(self):
        self.is_merged = self.p_meta['merged']
        self.is_blank = self.p_meta['blank']
        self.clustering_configs = constants.clustering_configs['merged' if self.is_merged else 'blank' if self.is_blank else 'single']

        # if the user updates the refinement of a single bin or bins, there shouldn't be multiple copies
        # of that stored in the database. so everytime 'store_refined_bins' function is called,
        # it will check this varlable and, (1) if empty, continue updating stuff in db store updates
        # in it, (2) if not empty, remove items stored in this variable from collections dict, and continue
        # with step (1). the starting point is of course self.bin_names_of_interest. when the store_refined_bins function is
        # called the first time, it will read collection data for collection_name, and remove the bin(s) in
        # analysis from it before it stores the data:
        self.ids_for_already_refined_bins = self.bin_names_of_interest

        self.input_directory = os.path.dirname(os.path.abspath(self.profile_db_path))

        self.run.info('Input directory', self.input_directory)
        self.run.info('Collection ID', self.collection_name)
        self.run.info('Number of bins', len(self.bin_names_of_interest))
        self.run.info('Number of splits', len(self.split_names_of_interest))

        if not len(self.split_names_of_interest):
            raise ConfigError("Something must have gone wrong :/ You are knee deep in refine mode, and have "
                              "zero split names to work with. You have to go back, Kate. Like, ALL the way "
                              "back to complaining to the developers on Discord or something :(")

        if not self.skip_hierarchical_clustering:
            item_orders = self.cluster_splits_of_interest()

            default_clustering_class = constants.merged_default if self.is_merged else constants.single_default

            default_item_order = dbops.get_default_item_order_name(default_clustering_class, item_orders)
        else:
            # even though we do not generate any hierarchical clustering,
            # we need to clear these, otherwise full mode trees will appear.
            item_orders = {}
            default_item_order = None

        self.item_orders = item_orders
        self.p_meta['item_orders'] = item_orders
        self.p_meta['available_item_orders'] = list(self.item_orders.keys())
        self.p_meta['default_item_order'] = default_item_order

        self.add_user_tree()

        if self.skip_hierarchical_clustering and not self.tree:
            raise ConfigError("You wanted to skip hierarchical clustering but did not provided any newick file.")

        self.collections = ccollections.Collections()
        self.collections.populate_collections_dict(self.profile_db_path)

        bins = sorted(list(self.bin_names_of_interest))

        if not self.args.title:
            self.title = textwrap.fill('Refining %s%s from "%s"' % (', '.join(bins[0:3]),
                                                      ' (and %d more)' % (len(bins) - 3) if len(bins) > 3 else '',
                                                      self.collection_name))


    def load_trnaseq_mode(self):
        # Item additional data will have already been loaded from table item_additional_data in the
        # profile databases.
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        trna_taxonomy_table = contigs_db.db.get_table_as_dataframe(t.trna_taxonomy_table_name,
                                                                   error_if_no_data=False).set_index('gene_callers_id')
        items_additional_data_dict = self.items_additional_data_dict
        items_additional_data_keys = self.items_additional_data_keys
        if len(trna_taxonomy_table):
            ranks = t.taxon_names_table_structure[1: ]
            for split_name, gene_callers_id in contigs_db.db.get_table_as_dataframe(t.genes_in_splits_table_name).loc[:, ['split', 'gene_callers_id']].values:
                try:
                    # Check that the split is taxonomically annotated.
                    taxonomy_series = trna_taxonomy_table.loc[gene_callers_id, :]
                except KeyError:
                    try:
                        # Check that the split already has additional data.
                        # Add placeholder entries for the taxonomic data.
                        item_dict = items_additional_data_dict[split_name]
                        for rank in ranks:
                            item_dict[rank[2: ]] = None
                    except KeyError:
                        continue
                    continue

                try:
                    # Check that the split already has additional data.
                    item_dict = items_additional_data_dict[split_name]
                except KeyError:
                    items_additional_data_dict[split_name] = item_dict = {}.fromkeys(items_additional_data_keys)

                for rank in ranks:
                    item_dict[rank[2: ]] = taxonomy_series.loc[rank] # convert "t_domain" to "domain", etc.

            items_additional_data_keys.extend([rank[2:] for rank in ranks])


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
            self.run.warning('HMMs for single-copy core genes were not run for this contigs database. So you will not '
                             'see completion / redundancy estimates in the collection mode as additional layers. SAD.')
            completion_redundancy_available = False
            self.hmm_access = None
        else:
            self.progress.new('Accessing HMM hits')
            self.progress.update('...')
            self.hmm_access = hmmops.SequencesForHMMHits(self.contigs_db_path, sources=set(completeness.sources))
            self.progress.end()

        # we are about to request a collections dict that contains only split names that appear in the
        # profile database along with other info:
        self.collection, bins_info_dict, split_names_in_db_but_missing_in_collection = \
                                        self.collections.get_trimmed_dicts(self.collection_name,
                                                                           utils.get_all_item_names_from_the_database(self.profile_db_path))

        # we will do something quite tricky here. first, we will load the full mode to get the self.views
        # data structure fully initialized based on the profile database. Then, we using information about
        # bins in the selected collection, we will create another views data structure, and replace it with
        # the one we have. that will be LOVELY.
        self.load_full_mode()

        self.p_meta['available_item_orders'] = []
        self.p_meta['item_orders'] = {}

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

            item_order_name = ':'.join([view, self.distance, self.linkage])
            self.p_meta['available_item_orders'].append(item_order_name)

            if len(d['dict']) == 1:
                clustering_data = '(%s);' % list(d['dict'].keys())[0]
            else:
                clustering_data = clustering.get_newick_tree_data_for_dict(d['dict'],
                                                                           distance=self.distance,
                                                                           linkage=self.linkage)

            self.p_meta['item_orders'][item_order_name] = {'type': 'newick', 'data': clustering_data}

            # clustering is done, we can get prepared for the expansion of the view dict
            # with new layers. Note that these layers are going to be filled later.
            if completion_redundancy_available:
                d['header'].extend(['percent_completion', 'percent_redundancy', 'matching_domain'])
            d['header'].extend(['bin_name', 'source'])

            views_for_collection[view] = d

        default_clustering_class = 'mean_coverage' if self.p_meta['merged'] else 'single'
        self.p_meta['default_item_order'] = get_default_item_order_name(default_clustering_class, self.p_meta['item_orders'])

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
                p_completion, p_redundancy, domain, domain_probabilities, info_text, _ = completeness.get_info_for_splits(set(self.collection[bin_id]))

            for view in self.views:
                self.views[view]['dict'][bin_id]['bin_name'] = bin_id

                if completion_redundancy_available:
                    matching_domain, matching_domain_confidence = sorted(domain_probabilities.items(), key = lambda x: x[1], reverse=True)[0]
                    self.views[view]['dict'][bin_id]['percent_completion'] = p_completion
                    self.views[view]['dict'][bin_id]['percent_redundancy'] = p_redundancy
                    self.views[view]['dict'][bin_id]['matching_domain'] = '%s (%.1f)' % (matching_domain, matching_domain_confidence)

                self.views[view]['dict'][bin_id]['source'] = bins_info_dict[bin_id]['source']

            current_bin += 1
        self.progress.end()

        # let empty some variables to avoid confusion downstream
        self.hmm_sources_info = {}
        self.split_sequences = None
        self.splits_taxonomy_dict = {}
        self.displayed_item_names_ordered = sorted(self.views[self.default_view]['dict'].keys())

        # set the title:
        R = lambda x: x.replace('-', ' ').replace('_', ' ')
        self.title = "Collection '%s' for %s" % (R(self.collection_name), R(self.p_meta['sample_id']))


    def load_pan_mode(self):
        if not self.pan_db_path:
            raise ConfigError("So you want to display a pan genome without a pan database? Anvi'o is "
                               "confused :/")

        PanSuperclass.__init__(self, self.args)

        self.init_gene_clusters()

        if not self.skip_init_functions:
            self.init_gene_clusters_functions()
            self.init_gene_clusters_functions_summary_dict()

        PanSuperclass.init_items_additional_data(self)

        self.p_meta['item_orders'] = self.item_orders

        self.load_pan_views()
        self.default_view = self.p_meta['default_view']

        self.collections.populate_collections_dict(self.pan_db_path)

        # set title for the pangenome
        if self.title:
            self.title = self.title
        else:
            self.title = self.p_meta['project_name'].replace('-', ' ').replace('_', ' ')

        # add user tree if there is one
        self.add_user_tree()


    def load_full_mode(self):
        if not self.contigs_db_path:
            raise ConfigError("Anvi'o needs the contigs database to make sense of this run (or maybe you "
                               "should use the `--manual` flag if that's what your intention).")

        if not self.profile_db_path:
            raise ConfigError("So you want to run anvi'o in full mode, but without a profile database? "
                               "Well. This does not make any sense.")

        if not self.skip_init_functions:
            self.init_functions()

        if self.mode == 'refine':
            # the refine trick: initialize profile super without a collection
            # name since we already know our split names:
            collection_name = self.collection_name
            self.collection_name = None
            self.args.collection_name = None
            ProfileSuperclass.__init__(self, self.args)
            self.collection_name = collection_name
            self.args.collection_name = collection_name
        else:
            ProfileSuperclass.__init__(self, self.args)

        # init item additional data
        ProfileSuperclass.init_items_additional_data(self)

        # this is a weird place to do it, but we are going to ask ContigsSuperclass function to load
        # all the split sequences since only now we know the min_contig_length that was used to profile
        # this stuff
        self.init_split_sequences(min_contig_length=self.p_meta['min_contig_length'])

        self.collections.populate_collections_dict(self.profile_db_path)

        self.p_meta['self_path'] = self.profile_db_path
        self.p_meta['output_dir'] = os.path.join(os.getcwd(), os.path.dirname(self.profile_db_path))

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path)

        # load views from the profile database
        if self.p_meta['blank']:
            blank_dict = {}
            for split_name in self.splits_basic_info:
                blank_dict[split_name] = {'blank_view': 0, '__parent__': self.splits_basic_info[split_name]['parent']}

            self.views['blank_view'] = {'header': ['blank_view', '__parent__'],
                                        'dict': blank_dict}
            self.default_view = 'blank_view'

        else:
            # Here we will ensure we are loading the views only for those splits that we need.
            # This may come across as odd: why would anyone need anything but all splits
            # in full mode? Well, there _is_ a reason: to avoid any redundancy in the code,
            # the refine mode is initiated right after the full mode. in refine mode anvi'o
            # always focuses on a subset of split names, and by checking for that here we are
            # essentially avoiding the loading of all splits in a profile database to build and
            # take care of yet another performance bottleneck Xabier Vázquez-Campos identified
            # in #1458:
            if len(self.split_names_of_interest):
                self.load_views(split_names_of_interest=self.split_names_of_interest)
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
                raise ConfigError("The last column of the additional view must be '__parent__' with the proper "
                                   "parent information for each split.")

            column_mapping = [str] + [float] * (len(additional_view_columns) - 1) + [str]

            self.views['user_view'] = {'table_name': 'NA',
                                       'header': additional_view_columns,
                                       'dict': utils.get_TAB_delimited_file_as_dictionary(self.additional_view_path, column_mapping=column_mapping)}

        # if the user specifies a view, set it as default:
        if self.view:
            if not self.view in self.views:
                raise ConfigError("The requested view ('%s') is not available for this run. Please see "
                                         "available views by running this program with --show-views flag." % self.view)

            self.default_view = self.view

        self.p_meta['item_orders'] = self.item_orders

        # add user tree if there is one (but skip the collection mode as it adds its own trees)
        if not self.mode == 'collection':
            self.add_user_tree()

        # set title
        if self.title:
            self.title = self.title
        else:
            self.title = self.p_meta['sample_id'].replace('-', ' ').replace('_', ' ')

        # do we have auxiliary data available?
        if not self.auxiliary_profile_data_available:
            summary_cp_available = os.path.exists(os.path.join(os.path.dirname(self.profile_db_path), 'SUMMARY.cp'))
            self.run.warning("Auxiliary data is not available; which means you will not be able to perform "
                             "certain operations (i.e., the inspect menu in the interactive interface will "
                             "not work, etc). %s" % ('' if not summary_cp_available else "Although, you have "
                             "a SUMMARY.cp file in your work directory, which means you are working with an "
                             "outdated anvi'o run. You can convert your SUMMARY.cp into an auxiliary data file "
                             "by using `anvi-script-generate-auxiliary-data-from-summary-cp` script."))

        if self.state_autoload:
            if not self.state_autoload in self.states_table.states:
                raise ConfigError("The requested state ('%s') is not available for this run. Please see "
                                         "available states by running this program with --show-states flag." % self.state_autoload)

        # did the user ask anvi'o to add an additional layer to summarize proportion of functions per split?
        self.update_items_additional_data_with_functions_per_split_summary()


    def load_gene_mode(self):
        if not self.skip_init_functions:
            self.init_functions()

        ProfileSuperclass.__init__(self, self.args)

        self.init_gene_level_coverage_stats_dicts()

        # the gene_level_coverage_stats_dict contains a mixture of data, some of which are not relevant to
        # our purpose of generating views for the interactive interface. here we explicitly list keys that
        # correspond to views we wish to generate:
        if self.inseq_stats:
            views_of_interest = ['mean_coverage', 'insertions', 'insertions_normalized', 'mean_disruption', 'below_disruption']
        else:
            views_of_interest = ['mean_coverage', 'detection', 'non_outlier_mean_coverage', 'non_outlier_coverage_std']

        for view in views_of_interest:
            self.views[view] = {
                'table_name': 'genes',
                'header': self.p_meta['samples'],
                'dict': {}
                }

        self.collections.populate_collections_dict(self.profile_db_path)
        all_gene_callers_ids = []

        # NOTE: here we ask anvi'o to not worry about names consistency. because we are initializing
        # split sequences in gene mode, downstream code will be confused regarding how are there items
        # in this object (gene caller IDs) that do not match split sequences (that are coming from the
        # contigs database). Bypassing those checks helps us avoid headaches later.
        self.skip_check_names = True

        self.init_split_sequences(self.p_meta['min_contig_length'], split_names_of_interest=self.split_names_of_interest)

        gene_caller_ids_missing_in_gene_level_cov_stats_dict = set([])
        gene_caller_sources_for_missing_gene_caller_ids_in_gene_level_cov_stats_dict = set([])
        for split_name in self.split_names_of_interest:
            genes_in_splits_entries = self.split_name_to_genes_in_splits_entry_ids[split_name]

            for genes_in_splits_entry in genes_in_splits_entries:
                e = self.genes_in_splits[genes_in_splits_entry]
                gene_callers_id = e['gene_callers_id']
                gene_caller_source = self.genes_in_contigs_dict[gene_callers_id]['source']

                if gene_callers_id not in self.gene_level_coverage_stats_dict:
                    gene_caller_ids_missing_in_gene_level_cov_stats_dict.add(gene_callers_id)
                    gene_caller_sources_for_missing_gene_caller_ids_in_gene_level_cov_stats_dict.add(gene_caller_source)
                    continue

                all_gene_callers_ids.append(gene_callers_id)

                for view in views_of_interest:
                    self.views[view]['dict'][str(gene_callers_id)] = {}
                    for sample_name in self.gene_level_coverage_stats_dict[gene_callers_id]:
                        self.views[view]['dict'][str(gene_callers_id)][sample_name] = self.gene_level_coverage_stats_dict[gene_callers_id][sample_name][view]

        if len(gene_caller_ids_missing_in_gene_level_cov_stats_dict):
            self.run.warning("Anvi'o observed something weird while it was processing gene level coverage statistics. "
                             "Some of the gene calls stored in your contigs database (%d of them, precisely) did not "
                             "have any information in gene level coverage stats dictionary. It is likely they were added "
                             "to the contigs database *after* these gene level coverage stats were computed and stored. "
                             "One way to address this is to remove the database file for gene coverage stats and re-run "
                             "this step. Another way to do it is the good ol' way of ignoring these warnings as you usually "
                             "do. If you do the latter, you will probably be fine. But we wanted to keep you in the loop "
                             "just in case you are not fine and the code does not realize that yet. The sources for gene "
                             "calls that will not be included in your gene view include these ones: '%s'." % \
                                    (len(gene_caller_ids_missing_in_gene_level_cov_stats_dict),
                                     ', '.join(gene_caller_sources_for_missing_gene_caller_ids_in_gene_level_cov_stats_dict)))

        # this is a bit tricky. we already populated the collections dict above, but it was to find out
        # genes caller ids of interest. we are now resetting the collections dict, and populating it
        # from scratch from the genes database. similarly, in store collection operations in bottleroutes
        # and elsewhere will target the gene database tables for gene mode.
        self.collections = ccollections.Collections()
        self.collections.populate_collections_dict(self.genes_db_path)

        # here we set the states table to the genes database to make sure the information from the interface
        # goes into the right table
        self.states_table = TablesForStates(self.genes_db_path)

        # dealing with item ordes in the genes database
        self.p_meta['available_item_orders'], self.p_meta['item_orders'] = dbops.get_item_orders_from_db(self.genes_db_path)

        # FIXME: These clustering data should be stored in the genes database :/ what is here is from the ad hoc
        #        attempt to run things in gene mode before we had the genes database, and we can do better than this
        #        now. But for that we also need to remove the 'gene_level_coverage_stats' table, and instead sotre
        #        these information into view tables. both pan and profile databases are nicely following that design,
        #        and there is no need to make genes database so different.
        for view in views_of_interest:
            item_order_name = view
            newick_tree_text = clustering.get_newick_tree_data_for_dict(self.views[view]['dict'],
                                                                        linkage=self.linkage,
                                                                        distance=self.distance)

            self.p_meta['available_item_orders'].append(item_order_name)
            self.p_meta['item_orders'][item_order_name] = {'type': 'newick', 'data': newick_tree_text}

        self.p_meta['item_orders']['synteny'] = {'type': 'basic', 'data': list(map(str, sorted(all_gene_callers_ids)))}

        self.p_meta['default_item_order'] = 'mean_coverage'
        self.default_view = 'mean_coverage'
        self.title = "Genes in '%s' [mode: %s]" % (self.bin_id, 'IN-Seq/Tn-Seq' if self.inseq_stats else 'Standard')

        # FIXME: When we are in gene-mode mode, our item names are no longer split names, hence the
        # following dictionaries are useless. Until we find a better way to fill them up with
        # potentially useful information, we can nullify them
        self.split_lengths_info = dict([(split_name, self.splits_basic_info[split_name]['length']) for split_name in self.splits_basic_info])
        self.splits_basic_info = {}
        self.splits_taxonomy_dict = {}
        self.p_meta['description'] = 'None'

        for view in views_of_interest:
            # don't bother if this is a single profile
            if not self.p_meta['merged']:
                continue

            data_value = clustering.get_newick_tree_data_for_dict(self.views[view]['dict'],
                                                                  distance=self.distance,
                                                                  linkage=self.linkage,
                                                                  transpose=True)

            # additional layers orders are added with the prefix "genes_"
            # this is because there are similar orders available from the contigs data
            #(for example contigs detections vs. genes detections)
            self.layers_order_data_dict['genes_' + view] = {'newick': data_value, 'basic': None}


    def add_user_tree(self):
        if self.tree:
            clustering_id = '%s:unknown:unknown' % filesnpaths.get_name_from_file_path(self.tree)

            if not self.p_meta['item_orders']:
                self.p_meta['default_item_order'] = clustering_id
                self.p_meta['available_item_orders'] = [clustering_id]
                self.p_meta['item_orders'] = {clustering_id: {'type': 'newick', 'data': open(os.path.abspath(self.tree)).read()}}
                run.info('Additional Tree', "Splits will be organized based on '%s'." % clustering_id)
            else:
                self.p_meta['item_orders'][clustering_id] = {'type': 'newick', 'data': open(os.path.abspath(self.tree)).read()}
                run.info('Additional Tree', "'%s' has been added to available trees." % clustering_id)

            # make sure item names in collection mode matches to the items in the user-provided tree
            if self.mode == 'collection':
                item_names_in_collection = set(self.collection.keys())
                item_names_in_user_tree = set(utils.get_names_order_from_newick_tree(self.tree, names_with_only_digits_ok=True))
                if not item_names_in_collection == item_names_in_user_tree:
                    if anvio.FORCE_USE_MY_TREE:
                        run.warning(f"Anvi'o told you that the tree file you provided had {len(item_names_in_user_tree)} names, which "
                                    f"the {len(item_names_in_collection)} names in the collection '{self.collection_name}', but you said "
                                    f"--force-use-my-tree. So that's what will happen now (and you will thank Obama if things go wrong)",
                                    header="YOU USED YOUR AUTHORITY ◔_◔")
                    else:
                        raise ConfigError(f"You are attempting to run the anvi'o interactive in collection mode, and you ALSO provide a "
                                          f"tree file to organize your items. Which is all great and this is exactly what anvi'o is for. "
                                          f"But it seems the {len(item_names_in_user_tree)} items in your tree file do not match to the "
                                          f"{len(item_names_in_collection)} item names the '{self.collection_name}' collection  describes :/ "
                                          f"In case it helps you solve this puzzle, here is a name that appears in your tree file: "
                                          f"'{item_names_in_user_tree.pop()}'. And here is a name that appears in your collection: "
                                          f"'{item_names_in_collection.pop()}'. If they look good to you, then the problem may be related to "
                                          f"items that are only in your collection and not in your tree, or vice versa. There should be a one "
                                          f"to one match between the two.")
            elif self.mode == 'pan':
                item_names_in_user_tree = set(utils.get_names_order_from_newick_tree(self.tree, names_with_only_digits_ok=True))
                if not item_names_in_user_tree == set(self.gene_cluster_names):
                    raise ConfigError(f"You are passing a tree to the interactive interface, which is a mechanism for you to organize your "
                                      f"gene clusters the way you like. More power to you! But there is a problem: the names in your tree do "
                                      f"not match the names of the gene clusters. For instance, here is a gene cluster name anvi'o knows about "
                                      f"your pangenome: '{self.gene_cluster_names.pop()}'. In contrast, here is an example to the kinds of names "
                                      f"anvi'o finds in your tree: '{item_names_in_user_tree.pop()}' (pro tip: if your intention with this tree "
                                      f"is to organize your GENOMES rather than your GENE CLUSTERS, than you need to use the program "
                                      f"`anvi-import-misc-data` to import a new layer order data into your pan-db. If this error message doesn't "
                                      f"make sense to you at all, don't despair: find us online and ask for help. You're on the right track).")
            else:
                pass


    def update_items_additional_data_with_functions_per_split_summary(self):
        """Adds a layer of stacked bar chart for proportion of functions in split for a given source"""

        if self.annotation_source_for_per_split_summary:
            self.items_additional_data_dict, self.items_additional_data_keys = \
                self.get_items_additional_data_for_functions_per_split_summary(self.annotation_source_for_per_split_summary,
                                                                               split_names_of_interest=self.split_sequences.keys(),
                                                                               data_dict=self.items_additional_data_dict,
                                                                               keys_list=self.items_additional_data_keys)


    def search_for_functions(self, search_terms, requested_sources=None):
        search_terms = [s.strip() for s in search_terms.split(',')]
        full_report = None

        if self.mode == 'full' or self.mode == 'gene' or self.mode == 'refine':
            items, full_report = ContigsSuperclass.search_for_gene_functions(self, search_terms, verbose=False, requested_sources=requested_sources)
            if self.mode == 'gene' or self.mode == 'refine':
                # otherwise gene mode report functions from other splits are not the bin interactive initialized.
                for item in items:
                    items[item] = [i for i in items[item] if i in self.split_names_of_interest]

                full_report = [i for i in full_report if i[5] in self.split_names_of_interest]
        elif self.mode == 'pan':
            items, full_report = PanSuperclass.search_for_gene_functions(self, search_terms, verbose=False, requested_sources=requested_sources)
        else:
            raise ConfigError("Searching functions are not supported for this mode.")

        return items, full_report


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
            raise ConfigError("There are some split names in your additional view data file ('%s') that are missing from "
                               "split names characterized in the database. There are in fact %d of them. For instance, "
                               "here is a random split name that is in your additional view data, yet not in the database: "
                               "'%s'. This is not going to work for anvi'o :/" \
                                    % (self.additional_view_path, len(splits_in_additional_view_but_not_in_tree), splits_in_additional_view_but_not_in_tree.pop()))

        if splits_in_tree_but_not_in_view_data:
            num_examples = 5 if len(splits_in_tree_but_not_in_view_data) >= 5 else len(splits_in_tree_but_not_in_view_data)
            example_splits_missing_in_view = [splits_in_tree_but_not_in_view_data.pop() for _ in range(0, num_examples)]
            raise ConfigError('Some split names found in your tree are missing in your view data. Hard to '
                               'know what cuased this, but essentially your tree and your view does not '
                               'seem to be compatible. Here is a couple of splits that appear in the tree '
                               'but not in the view data: %s.' % ', '.join(example_splits_missing_in_view))

        if splits_in_tree_but_not_in_database:
            raise ConfigError('Some split names found in your tree are missing from your database. Hard to '
                               'know why is this the case, but here is a couple of them: %s'\
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
                run.warning("Some of the contigs in your addtional view data file does not appear to be in anywhere else. "
                            "Additional layers file is not required to have data for all items (which means, there may be "
                            "items in view dictionaries that are not in the additional view data file), however, finding "
                            "items that are only in the additional layers file usually means trouble. Anvi'o will continue, "
                            "but please go back and check your files if you think there may be something wrong. Here is a "
                            "random item name that was only in your file: '%s'. And there were %d of them in total. So you "
                            "are warned!" % (one_example, num_all))


    def prune_view_and_additional_data_dicts(self):
        self.progress.new('Pruning view dicts')
        self.progress.update('...')

        splits_in_views = set(list(self.views.values())[0]['dict'].keys())
        splits_to_remove = splits_in_views - set(self.displayed_item_names_ordered)

        for view in self.views:
            self.progress.update('processing view "%s"' % view)
            for split_name in splits_to_remove:
                self.views[view]['dict'].pop(split_name)

        self.progress.update("working on items additional data")
        splits_in_items_additional_data_dict = set(list(self.items_additional_data_dict.keys()))
        splits_to_remove = splits_in_items_additional_data_dict - set(self.displayed_item_names_ordered)

        for split_name in splits_to_remove:
            self.items_additional_data_dict.pop(split_name)

        # if you remove all splits from the additional data dict (say because you are in)
        # collection mode and your `items` are bin names and not split names), make sure
        # the items additional data keys variable reflects that fact (reported by
        # Florentin Constancias / @fconstancias in #1705):
        if not len(self.items_additional_data_dict):
            self.items_additional_data_keys = []

        self.progress.end()


    def get_layer_names_with_at_least_one_hit_for_splits(self, layer_names_list, splits_dict):
        """This is a funciton to remove layers with no hits from default displays."""

        if self.show_all_layers:
            return layer_names_list

        layer_names_with_hits = []
        layer_names_to_consider = copy.deepcopy(layer_names_list)
        layer_names_considered = set([])

        layer_types = {}
        for layer_name in layer_names_to_consider:
            layer_types[layer_name] = utils.get_predicted_type_of_items_in_a_dict(splits_dict, layer_name)

        # while this looks extremely inefficient, it is not that bad, in fact. This loop simply goes thorugh
        # every layer name that is in `layer_names_to_consider` until it hits the first non-None value that
        # has a length or larger than 0 given the type of the layer. The fate of the well-populated layers
        # is determined very rapidly. What takes much longer are those that have no values for any splits
        # or values that are very sparse, and they will not be common to most profiles.
        for hits in splits_dict.values():
            for layer_name in layer_names_to_consider:
                if layer_name in layer_names_considered:
                    continue

                if layer_types[layer_name] == None:
                    # all items in this layer has type None, there can't possibly
                    # be anything worth showing here.
                    layer_names_considered.add(layer_name)
                    continue

                if layer_name in hits:
                    if layer_types[layer_name] in [int, float]:
                        if hits[layer_name] and float(hits[layer_name]) > 0:
                            layer_names_with_hits.append(layer_name)
                            layer_names_considered.add(layer_name)
                    else:
                        # which means we have layer type of str, dict, or list
                        if hits[layer_name] and len(hits[layer_name]):
                            layer_names_with_hits.append(layer_name)
                            layer_names_considered.add(layer_name)

        # update our class variable so we know what is not going to be shown
        self.layers_that_will_not_be_shown.update(set(layer_names_list) - set(layer_names_with_hits))

        return layer_names_with_hits


    def finalize_view_data(self):
        '''Finalize what will be shown in the interactive interface.

        This includes expanding available 'view's with additional layers to be shown,
        including HMM hits, or user-provided data through the command line or misc
        data tables.
        '''

        # we will fill this one up through `self.get_layer_names_with_at_least_one_hit_for_splits`
        # to warn the user.
        self.layers_that_will_not_be_shown = set([])

        for view in self.views:
            # here we will populate runinfo['views'] with json objects.
            view_dict = self.views[view]['dict']
            view_headers = self.views[view]['header']

            json_object = []

            # (1) set the header line with the first entry:
            json_header = ['contigs']

            # (2) then add length and GC content IF we have sequences available
            if self.splits_basic_info:
                basic_info_headers = ['length', 'gc_content']
                json_header.extend(basic_info_headers)

            # (3) then add the view!
            json_header.extend(view_headers)

            # (4) then add 'additional' headers as the outer ring:
            if self.items_additional_data_keys:
                self.items_additional_data_keys = sorted(self.get_layer_names_with_at_least_one_hit_for_splits(self.items_additional_data_keys, self.items_additional_data_dict))
                json_header.extend(self.items_additional_data_keys)

            # (5) finally add hmm search results
            if self.hmm_searches_dict:
                if not self.split_hmm_layers:
                    layer_names_for_hmm_hits = self.get_layer_names_with_at_least_one_hit_for_splits([tpl[1] for tpl in self.hmm_searches_header], self.hmm_searches_dict)

                    # because keys for HMMs are weird to make sure we can easily switch between
                    # the regular form of displaying HMMs and the form comes with the user flag
                    # `self.split_hmm_layers`, we have to do something more here to clean up
                    # the global `self.hmm_searches_header` list that is used again down
                    # below:
                    self.hmm_searches_header = [tpl for tpl in self.hmm_searches_header if tpl[1] in layer_names_for_hmm_hits]
                else:
                    # if the user asked for HMM layers to be split, we don't want to do anything.
                    # it gets too prickly otherwise.
                    pass

                # let's sort these layers alphabetically
                self.hmm_searches_header = sorted(self.hmm_searches_header)

                json_header.extend([tpl[0] for tpl in self.hmm_searches_header])

            # (6) add taxonomy, if exitsts:
            if len(self.splits_taxonomy_dict):
                json_header.extend(['taxonomy'])

            # (7) and finalize it (yay):
            json_object.append(json_header)

            for split_name in view_dict:
                # (1)
                json_entry = [split_name]

                # (2)
                if self.splits_basic_info:
                    json_entry.extend([self.splits_basic_info[split_name][header] for header in basic_info_headers])

                # (3) adding essential data for the view
                json_entry.extend([view_dict[split_name][header] if header in view_dict[split_name] else 0 for header in view_headers])

                # (4) adding additional layers
                if self.items_additional_data_keys:
                    json_entry.extend([self.items_additional_data_dict[split_name][header] if split_name in self.items_additional_data_dict else None for header in self.items_additional_data_keys])

                # (5) adding hmm stuff
                if self.hmm_searches_dict:
                    if self.split_hmm_layers:
                        json_entry.extend([self.hmm_searches_dict[split_name][header] if split_name in self.hmm_searches_dict else None for header in [tpl[0] for tpl in self.hmm_searches_header]])
                    else:
                        json_entry.extend([len(self.hmm_searches_dict[split_name][header]) if split_name in self.hmm_searches_dict else 0 for header in [tpl[1] for tpl in self.hmm_searches_header]])

                # (6)
                if self.splits_taxonomy_dict:
                    if split_name in self.splits_taxonomy_dict:
                        json_entry.extend([self.splits_taxonomy_dict[split_name]])
                    else:
                        json_entry.extend([None])

                # (7) send it along!
                json_object.append(json_entry)

            self.views[view] = json_object

        if len(self.layers_that_will_not_be_shown):
            self.run.warning(f"Even though the following layer names were associated with your data, they "
                             f"will not be shown in your interactive display since none of the splits you "
                             f"are in terested at this point had any hits in those layers: "
                             f"\"{', '.join(sorted(list(self.layers_that_will_not_be_shown)))}\". If you "
                             f"would like every layer to be shown in the interface even when there are no"
                             f"hits, please include the flag `--show-all-layers` in your command.")


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
            raise RefineError('%s of your bin names %s NOT unique, and already exist%s in the database. You must rename '
                               '%s to something else: %s' % (
                                                              'One' if len(bad_bin_names) == 1 else len(bad_bin_names),
                                                              'is' if len(bad_bin_names) == 1 else 'are',
                                                              's' if len(bad_bin_names) == 1 else '',
                                                              'this one' if len(bad_bin_names) == 1 else 'these',
                                                              ', '.join(bad_bin_names)
                                                             ))

        collections = TablesForCollections(self.profile_db_path)

        # remove bins that should be updated in the database:
        for bin_id in self.ids_for_already_refined_bins:
            collections.delete_bin(self.collection_name, bin_id)

        # zero it out
        self.ids_for_already_refined_bins = set([])

        if anvio.DEBUG:
            run.info('collection from db', f"{collection_dict}")
            run.info('bins info from db', f"{bins_info_dict}")
            run.info_single('')

            run.info('incoming collection data', f"{refined_bin_data}")
            run.info('incoming bins info', f"{refined_bins_info_dict}")
            run.info_single('')

        for bin_id in refined_bin_data:
            self.ids_for_already_refined_bins.add(bin_id)

        if anvio.DEBUG:
            run.info('resulting collection', f"{collection_dict}")
            run.info('resulting bins info', f"{bins_info_dict}")
            run.info_single('')

        collections.append(self.collection_name, refined_bin_data, refined_bins_info_dict, drop_collection=False)

        collections.refresh_collections_info_table(self.collection_name)

        run.info_single('"%s" collection is updated!' % self.collection_name)


    def end(self):
        # FIXME: remove temp files and stuff
        pass


    def get_anvio_news(self):
        if self.skip_news:
            return
        else:
            try:
                self.anvio_news = utils.get_anvio_news()
            except:
                pass


class StructureInteractive(VariabilitySuper, ContigsSuperclass):
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.args = args
        self.mode = 'structure'

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        # database
        self.profile_db_path = A('profile_db', null)
        self.contigs_db_path = A('contigs_db', null)
        self.structure_db_path = A('structure_db', null)
        # genes
        self.available_gene_caller_ids = A('gene_caller_ids', null)
        self.available_genes_path = A('genes_of_interest', null)
        self.available_genes = A('available_genes', list) or []
        # samples
        self.available_samples_path = A('samples_of_interest', null)
        self.available_samples = A('available_samples', set) or set([])
        # others
        self.saavs_only = A('SAAVs_only', bool)
        self.scvs_only = A('SCVs_only', bool)
        self.variability_table_path = A('variability_profile', null)
        self.min_departure_from_consensus = A('min_departure_from_consensus', float) or 0

        # states
        self.states_table = TablesForStates(self.structure_db_path)

        # For now, only true if self.variability_table_path. Otherwise variability is computed on the fly
        self.store_full_variability_in_memory = True if self.variability_table_path else False
        self.full_variability = None
        self.variability_storage = {}

        self.num_reported_frequencies = 5

        self.sanity_check()

        ContigsSuperclass.__init__(self, self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        # Init the amino acid additional data table object
        args = argparse.Namespace(contigs_db=self.contigs_db_path)
        self.amino_acid_additional_data = TableForAminoAcidAdditionalData(args)

        if self.store_full_variability_in_memory:
            self.profile_full_variability_data()

        self.available_genes = self.get_available_genes()
        self.available_engines = self.get_available_engines()

        # can save significant memory if available genes is a fraction of genes in full variability
        if self.full_variability:
            self.filter_full_variability()
            self.available_genes = self.get_available_genes() # genes can be lost during self.filter_full_variability()
                                                              # so we re-calculate available genes :\
            self.process_full_variability()

        # default gene is the first gene of interest
        self.profile_gene_variability_data(self.available_genes[0])

        self.available_samples = self.get_available_samples()
        self.sample_groups = self.create_sample_groups_dict()


    def get_gene_function_info(self, gene_callers_id):
        """Returns gene function info from the gene_function_calls_dict of ContigsSuperclass

        If gene functions have not been initialized, ContigsSuperclass.init_functions is called

        Notes
        =====
        - This initializes genes for the entire contigs database. When this becomes intolerably
          inefficient, ContigsSuperclass.init_functions should be modified to take a genes of
          interest flag
        """

        if not self.gene_function_calls_initiated:
            self.init_functions()

        return {'functions': self.gene_function_calls_dict.get(gene_callers_id)}


    def get_model_info(self, gene_callers_id):
        """Returns information about the protein model, e.g. template IDs, DOPE_score, etc"""

        structure_db = structureops.StructureDatabase(self.structure_db_path, 'none', ignore_hash=True)

        models = structure_db.db.get_table_as_dataframe(
            'models',
            columns_of_interest=['GA341_score', 'DOPE_score', 'molpdf', 'picked_as_best'],
            where_clause='corresponding_gene_call = %d' % gene_callers_id,
        ).rename(columns={
            'DOPE_score': 'DOPE',
            'GA341_score': 'GA341',
        })
        models['Models tried'] = models.shape[0]
        models = models.loc[models['picked_as_best'] == 1, ['DOPE', 'GA341', 'molpdf', 'Models tried']].reset_index(drop=True)

        templates = structure_db.db.get_table_as_dataframe(
            'templates',
            columns_of_interest=['pdb_id', 'chain_id'],
            where_clause='corresponding_gene_call = %d' % gene_callers_id,
        ).rename(columns={
            'pdb_id': 'PDB',
            'chain_id': 'Chain',
        })

        # Check if the columns exist before accessing them
        if 'percent_similarity' in templates.columns and 'align_fraction' in templates.columns:
            templates['%Identity'] = templates['percent_similarity']
            templates['Align fraction'] = templates['align_fraction']
            templates = templates[['PDB', 'Chain', '%Identity', 'Align fraction']]
        else:
            # Handle the case when columns are not present
            # You can set default values or handle it in a way that makes sense for your application
            templates['%Identity'] = 0
            templates['Align fraction'] = 0
            templates = templates[['PDB', 'Chain', '%Identity', 'Align fraction']]


        structure_db.disconnect()

        return {
            'models': models.to_json(orient='index'),
            'templates': templates.to_json(orient='index'),
        }


    def get_search_results_for_gene_functions(self, search_terms):
        """FIXME Currently unused"""
        items, full_report = ContigsSuperclass.search_for_gene_functions(self, search_terms, verbose=True)


    def filter_full_variability(self):
        try:
            self.full_variability.filter_data(criterion="corresponding_gene_call",
                                              subset_filter=self.available_genes)
        except self.EndProcess:
            raise ConfigError("This is really sad. There is no overlap between the gene IDs in your "
                              "structure database and the gene IDs in your variability table.")

        try:
            self.full_variability.filter_data(criterion="departure_from_consensus",
                                              min_filter=self.min_departure_from_consensus)
        except self.EndProcess:
            raise ConfigError("This is really sad. There are no entries in your variability table "
                              "with a departure_from_consensus less than {}. Try setting "
                              "--min-departure-from-consensus to 0.".format(self.min_departure_from_consensus))


    def process_full_variability(self):
        self.full_variability.convert_counts_to_frequencies()


    def create_sample_groups_dict(self):
        """This method converts the native data structure of additional_layer_dict:

           {sample1: {"temperature": "hot", "diet": "high-fat", ...},
            sample2: {"temperature": "hot", "diet": "low-fat", ...}}

           to a data structure designed with the structure interactive interface in mind:

           {"temperature": {"hot":      [sample1, sample2], "cold":    []},
            "diet":        {"high-fat": [sample1],          "low-fat": [sample2]}}
        """
        if not self.profile_db_path:
            return {"merged"  : {"all": sorted(list(self.available_samples))},
                    "samples" : {s:[s] for s in self.available_samples}}

        layer_names, additional_layer_dict = self.load_additional_layer_data()

        d = {}
        # important merged group (all samples in one group)
        d["merged"] = {"all": sorted(list(self.available_samples))}

        # filter additional_layer_dict to only include available samples
        additional_layer_dict = {k: v for k, v in additional_layer_dict.items() if k in self.available_samples}

        # Now begins the horrible process of converting between the two data structures. get ready
        # for triple-nested `for` loops and an unreadable mess. Or... a cheeky little pandas dataframe
        # whose only purpose is to enable a dictionary comprehension. Imagine how many lines this
        # would be using native python.
        df = pd.DataFrame(additional_layer_dict).T.fillna("NA").reset_index().rename(columns={"index":"samples"})
        d.update({col: {val: list(df.loc[df[col]==val, "samples"]) for val in numpy.sort(df[col].unique())} for col in df.columns}) # V/\

        return d


    def load_additional_layer_data(self, profile_db_path=None):
        """Get additional layer data for grouping

        Notes
        =====
        - Only includes layers if they are string-like
        """

        if not profile_db_path:
            profile_db_path = self.profile_db_path

        x = argparse.Namespace(pan_or_profile_db=profile_db_path, target_data_table="layers")
        additional_layers_table = TableForLayerAdditionalData(args=x)

        layer_names, additional_layer_dict = additional_layers_table.get()

        # The remainder of this method filters layer_names and additional_layer_dict to only include
        # layers that have string-like data. This is because we group based on these values, which
        # would not make sense with integer data

        samples_in_layer_data = additional_layer_dict.keys()
        layers_to_remove = []

        for layer_name in layer_names:
            # Assume we are removing this layer until we find evidence to keep it
            remove_column = True

            # If there is a "!" it is of stacked-bar type
            if not "!" in layer_name:
                for sample_name in samples_in_layer_data:
                    # loops through samples until it finds evidence the column is string-type
                    if isinstance(additional_layer_dict[sample_name][layer_name], str):
                        remove_column = False
                        continue

            if remove_column:
                layers_to_remove.append(layer_name)

        for sample_name in additional_layer_dict:
            for bad_layer in layers_to_remove:
                del additional_layer_dict[sample_name][bad_layer]

        return layer_names, additional_layer_dict


    def get_column_info(self, gene_callers_id, engine):
        var = self.variability_storage[gene_callers_id][engine]['var_object']
        if var.data.empty:
            return []

        FIND_MIN = lambda c, buff=0.01: var.data[c].min() - buff if c in var.data.columns else 0
        FIND_MAX = lambda c, buff=0.01: var.data[c].max() + buff if c in var.data.columns else 1

        info = [
            {
                'name': 'departure_from_consensus',
                'title': 'Departure from consensus',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 0.01,
                'min': 0,
                'max': 1,
            },
            {
                'name': 'departure_from_reference',
                'title': 'Departure from reference',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 0.01,
                'min': 0,
                'max': 1,
            },
            {
                'name': 'n2n1ratio',
                'title': 'Ratio of 2nd to 1st',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 0.01,
                'min': float(FIND_MIN('n2n1ratio')),
                'max': float(FIND_MAX('n2n1ratio'))
            },
            {
                'name': 'prevalence',
                'title': 'Prevalence',
                'as_view': True,
                'as_filter': False,
                'merged_only': True,
                'data_type': 'float',
                'min': 0,
                'max': 1
            },
            {
                'name': 'occurrence',
                'title': 'occurrence',
                'as_view': False,
                'as_filter': False,
                'merged_only': True,
                'data_type': 'integer',
            },
            {
                'name': 'contact_numbers',
                'as_view': False,
                'as_filter': False,
                'data_type': 'text',
            },
            {
                'name': 'mean_normalized_coverage',
                'title': 'Site coverage normalized by gene average',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 0.01,
                'min': float(FIND_MIN('mean_normalized_coverage')),
                'max': float(FIND_MAX('mean_normalized_coverage'))
            },
            {
                'name': 'coverage',
                'title': 'Site coverage',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1,
                'min': int(FIND_MIN('coverage', buff=1)),
                'max': int(FIND_MAX('coverage', buff=1))
            },
            {
                'name': 'log_pN_popular_consensus',
                'title': 'log10(pN) [popular consensus]',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1e-3,
                'min': float(FIND_MIN('log_pN_popular_consensus', buff=1e-3)),
                'max': float(FIND_MAX('log_pN_popular_consensus', buff=1e-3))
            },
            {
                'name': 'log_pN_consensus',
                'title': 'log10(pN) [consensus]',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1e-3,
                'min': float(FIND_MIN('log_pN_consensus', buff=1e-3)),
                'max': float(FIND_MAX('log_pN_consensus', buff=1e-3))
            },
            {
                'name': 'log_pN_reference',
                'title': 'log10(pN) [reference]',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1e-3,
                'min': float(FIND_MIN('log_pN_reference', buff=1e-3)),
                'max': float(FIND_MAX('log_pN_reference', buff=1e-3))
            },
            {
                'name': 'log_pS_popular_consensus',
                'title': 'log10(pS) [popular consensus]',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1e-3,
                'min': float(FIND_MIN('log_pS_popular_consensus', buff=1e-3)),
                'max': float(FIND_MAX('log_pS_popular_consensus', buff=1e-3))
            },
            {
                'name': 'log_pS_consensus',
                'title': 'log10(pS) [consensus]',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1e-3,
                'min': float(FIND_MIN('log_pS_consensus', buff=1e-3)),
                'max': float(FIND_MAX('log_pS_consensus', buff=1e-3))
            },
            {
                'name': 'log_pS_reference',
                'title': 'log10(pS) [reference]',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1e-3,
                'min': float(FIND_MIN('log_pS_reference', buff=1e-3)),
                'max': float(FIND_MAX('log_pS_reference', buff=1e-3))
            },
            {
                'name': 'entropy',
                'title': 'Entropy',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 0.01,
                'min': 0,
                'max': float(FIND_MAX('entropy'))
            },
            {
                'name': 'rel_solvent_acc',
                'title': 'Relative solvent accessibility',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 0.01,
                'min': 0,
                'max': 1,
            },
            {
                'name': 'sec_struct',
                'title': 'Secondary structure',
                'as_view': True,
                'as_filter': 'checkbox',
                'data_type': 'text',
                'choices': ['C', 'S', 'G', 'H', 'T', 'I', 'E', 'B']
            },
            {
                'name': 'phi',
                'title': 'Phi',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1,
                'min': -180,
                'max': 180,
            },
            {
                'name': 'psi',
                'title': 'Psi',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'float',
                'step': 1,
                'min': -180,
                'max': 180,
            },
            {
                'name': 'BLOSUM62',
                'title': 'BLOSUM62',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'integer',
                'step': 1,
                'min': -4,
                'max': 11,
            },
            {
                'name': 'BLOSUM90',
                'title': 'BLOSUM90',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'integer',
                'step': 1,
                'min': -6,
                'max': 11,
            },
            {
                'name': 'codon_order_in_gene',
                'title': 'Codon index',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'integer',
                'step': 1,
                'min': 0,
                'max': var.data["gene_length"].iloc[0]/3 - 1,
            },
            {
                'name': 'codon_number',
                'title': 'Codon number',
                'as_view': True,
                'as_filter': 'slider',
                'data_type': 'integer',
                'step': 1,
                'min': 1,
                'max': var.data["gene_length"].iloc[0]/3,
            },
            {
                'name': var.competing_items,
                'title': 'Competing Amino Acids' if engine == "AA" else 'Competing Codons',
                'as_view': True,
                'as_filter': 'checkbox',
                'data_type': 'text',
                'choices': list(var.data[var.competing_items].value_counts().sort_values(ascending=False).index)
            },
            {
                'name': 'reference',
                'title': 'Reference',
                'as_view': True,
                'as_filter': 'checkbox',
                'data_type': 'text',
                'choices': list(var.data['reference'].value_counts().sort_values(ascending=False).index)
            },
            {
                'name': 'consensus',
                'title': 'Consensus',
                'as_view': True,
                'as_filter': 'checkbox',
                'data_type': 'text',
                'choices': list(var.data['consensus'].value_counts().sort_values(ascending=False).index)
            },
        ]

        # coverage values
        for item in var.items:
            info.append(
                {
                    'name': item,
                    'as_view': False,
                    'as_filter': False,
                    'data_type': 'integer',
                },
            )

        # top n freqs
        for x in range(self.num_reported_frequencies):
            info.extend([
                {
                    'name': str(x) + '_item',
                    'as_view': False,
                    'as_filter': False,
                    'data_type': 'text',
                    'merged_only': True,
                },
                {
                    'name': str(x) + '_item_AA',
                    'as_view': False,
                    'as_filter': False,
                    'data_type': 'text',
                    'merged_only': True,
                },
                {
                    'name': str(x) + '_freq',
                    'as_view': False,
                    'as_filter': False,
                    'data_type': 'float',
                    'merged_only': True,
                },
            ])

        # keep those that the variability data table has. also keep those with no controller. these
        # entries may not exist in table, but will when data is merged
        info = [v for v in info if v['name'] in var.data.columns or v.get('merged_only', False)]
        return info


    def get_genes_of_interest_from_bin_id(self):
        pass


    def get_available_engines(self):
        if self.variability_table_path:
            engine = self.full_variability.engine
            if engine == 'NT':
                self.run.warning("You passed a variabilty table of SNVs (`--engine NT`). We did not "
                                 "develop this program with the intention of displaying NT variants "
                                 "yet it very randomly works and we're pretty happy about it. So "
                                 "have fun, yet understand some features won't work.")
            return [engine]
        elif self.saavs_only:
            return ['AA']
        elif self.scvs_only:
            return ['CDN']
        else:
            return ['AA', 'CDN']


    def get_available_genes(self, available_genes_path="", available_gene_caller_ids=""):
        """It is essential to note that available_genes_path and available_gene_caller_ids are "" by
           default, not None. The programmer can pass None to avoid the argument defaulting to a
           class-wide attribute (first code block of this method), which may not exist for classes
           inheriting this method.
        """
        # use class-wide attributes if no parameters are passed
        if available_gene_caller_ids == "" and available_genes_path == "":
            available_gene_caller_ids = self.available_gene_caller_ids
            available_genes_path = self.available_genes_path

        if available_gene_caller_ids:
            available_gene_caller_ids = [int(gene) for gene in available_gene_caller_ids.split(',')]

        # method inherited from VariabilitySuper
        requested_available_genes = self.get_genes_of_interest(available_genes_path, available_gene_caller_ids)

        # load in structure to compare genes of interest with those in requested_available_genes
        structure_db = structureops.StructureDatabase(self.structure_db_path, 'none', ignore_hash=True)

        if requested_available_genes:

            # check for genes that do not appear in the structure database
            unrecognized_genes = [g for g in requested_available_genes if g not in structure_db.genes_queried]
            if unrecognized_genes:
                some_to_report = unrecognized_genes[:5] if len(unrecognized_genes) <= 5 else unrecognized_genes
                raise ConfigError("{} of the gene caller ids you provided {} not known to the "
                                  "structure database. {}: {}. They might exist in the contigs "
                                  "database used to generate the structure database, but not all "
                                  "genes in the contigs database exist in the corresponding "
                                  "structure database :(. You can try running "
                                  "anvi-gen-structure-database again, this time with these missing "
                                  "genes included".format(len(unrecognized_genes),
                                                 "is" if len(unrecognized_genes) == 1 else "are",
                                                 "Here are a few of those ids" if len(some_to_report) > 1 else "Its id is",
                                                 ", ".join([str(x) for x in some_to_report])))

            # check for genes that structures were attempted for, but failed
            unavailable_genes = [g for g in requested_available_genes if g in structure_db.genes_without_structure]
            available_genes = [g for g in requested_available_genes if g not in unavailable_genes]
            if unavailable_genes:
                some_to_report = unavailable_genes[:5] if len(unavailable_genes) <= 5 else unavailable_genes
                run.warning("When your structure database was first generated "
                            "(anvi-gen-structure-database), anvi'o attempted--but failed--to predict "
                            "the structures for genes with the following ids: {}. Yet, these genes "
                            "are all included in your genes of interest. This is just a heads up so "
                            "that you're not surprised when these genes don't show up in your "
                            "display.".format(", ".join([str(x) for x in some_to_report])))

        else:
            available_genes = list(structure_db.genes_with_structure)

        # if full variability table is loaded, we further demand variability exists for the genes
        if self.full_variability:
            genes_in_variability = self.full_variability.data["corresponding_gene_call"].unique()
            available_genes = [x for x in available_genes if x in genes_in_variability]

        structure_db.disconnect()
        # We are done with self.args.gene_caller_ids and self.args.genes_of_interest, so we set them
        # to None. Now downstream class instances initialized with self.args will not process our
        # already processed gene caller ids or genes of interest path.
        self.args.gene_caller_ids = None
        self.args.genes_of_interest = None

        return sorted(available_genes)


    def get_available_samples(self, available_samples_path=""):
        """There may be confusion within this method regarding the difference between available
           samples and samples of interest. Available samples refer to those samples which are
           capable of being displayed, whereas samples of interest refers to a subset of the
           available samples that have been explicitly selected from within the interactive
           interface. This is despite the fact that users modify available samples using the flag
           --samples-of-interest.

           It is essential to note that available_samples_path is "" by default, not None. The
           programmer can pass None to avoid the argument defaulting to a class-wide attribute
           (first code block of this method), which may not exist for classes inheriting this
           method.
        """
        # use class-wide attributes if no parameters are passed
        if available_samples_path == "":
            available_samples_path = self.available_samples_path

        # method inherited from VariabilitySuper
        available_samples = self.get_sample_ids_of_interest(available_samples_path)

        if self.full_variability:
            all_sample_ids = self.full_variability.data["sample_id"].unique()
        else:
            profile_db = dbops.ProfileDatabase(self.profile_db_path)
            all_sample_ids = sorted(list(profile_db.samples))
            profile_db.disconnect()

        if available_samples:

            # check for samples that do not appear in the profile database
            unrecognized_samples = [g for g in available_samples if g not in all_sample_ids]
            if unrecognized_samples:
                raise ConfigError("{} of the sample ids you provided are not in the variability "
                                  "table. Here they are: {}. They may exist in the profile database "
                                  "that generated the variability table "
                                  "(anvi-gen-variability-profile), but were filtered out at some "
                                  "point. Or you made a mistake, but what are the chances of "
                                  "that?".format(len(unrecognized_samples),
                                                 ", ".join([str(x) for x in unrecognized_samples])))

        else:
            available_samples = set(all_sample_ids)

        # We are done with self.args.samples_of_interest, so we set it to None. Now downstream class
        # instances initialized with self.args will not process our already processed samples of
        # interest path.
        self.args.samples_of_interest = None

        return available_samples


    def sanity_check(self):
        if not self.structure_db_path:
            raise ConfigError("Must provide a structure database.")
        utils.is_structure_db(self.structure_db_path)

        if not self.contigs_db_path:
            raise ConfigError("Must provide a contigs database.")
        utils.is_contigs_db(self.contigs_db_path)

        if not self.profile_db_path and not self.variability_table_path:
            raise ConfigError("You have to provide either a variability table generated from "
                              "anvi-gen-variability-profile, or a profile database from "
                              "which sequence variability will be computed.")

        if self.variability_table_path:
            run.warning("You opted to work with a variability table previously generated from "
                        "anvi-gen-variability-profile. As a word of caution, keep in mind that any "
                        "filters applied when the table was generated now persist in the "
                        "following visualizations.")
            if not self.profile_db_path:
                run.warning("Anvi'o offers a nice way to visualize sequence variability for "
                            "user-defined groups of samples, as opposed to individual samples. "
                            "However, these groups are defined as additional layers in the profile "
                            "database used to generate your variability table "
                            "(anvi-gen-variability-profile) which you did not provide. If you "
                            "already have these tables in your profile database, include your "
                            "profile database with the flag `-p`. If you want to create groupings, "
                            "you can read about how to create them here: "
                            "http://merenlab.org/2017/12/11/additional-data-tables/#layers-additional-data-table.")

        if self.saavs_only and self.scvs_only:
            raise ConfigError("--SAAVs-only and --SCVs-only are not compatible with one another. Pick one.")


    def profile_full_variability_data(self):
        """Creates self.full_variability, which houses the full variability... well, the full
           variability of all genes with structures in the structure database
        """
        verbose, stealth = (True, False) if anvio.DEBUG else (False, True)

        self.progress.new("Loading full variability table"); self.progress.update("...")
        self.full_variability = variabilityops.VariabilityData(self.args, p=terminal.Progress(verbose=verbose), r=terminal.Run(verbose=verbose))
        self.full_variability.stealth_filtering = stealth
        self.progress.end()


    def profile_gene_variability_data(self, gene_callers_id):
        """Variability data is computed on the fly if a profile and contigs database is provided.
           If the variability table is provided, the full table is stored in memory and a gene
           subset is created.
        """
        if anvio.DEBUG:
            verbose, stealth = True, False
        else:
            verbose, stealth = False, True

        if gene_callers_id in self.variability_storage:
            # already profiled.
            return

        self.progress.new("Profiling variability for gene %s" % str(gene_callers_id))
        gene_var = {}
        for engine in self.available_engines:
            self.progress.update("Fetching %s..." % ("SAAVs" if engine == "AA" else "SCVs"))
            gene_var[engine] = {}

            if self.store_full_variability_in_memory:
                # if the full variability is in memory, make a deep copy, then filter
                var = copy.deepcopy(self.full_variability)
                var.filter_data(criterion="corresponding_gene_call", subset_filter=set([gene_callers_id]))
            else:
                # if not, we profile from scratch, passing as an argument our gene of interest
                self.args.engine = engine
                self.args.genes_of_interest_set = set([gene_callers_id])
                self.args.compute_gene_coverage_stats = True
                var = variability_engines[engine](self.args, p=terminal.Progress(verbose=verbose), r=terminal.Run(verbose=verbose))
                var.stealth_filtering = stealth

                if engine == 'CDN':
                    var.process_functions.append((var.calc_pN_pS, dict(grouping='site', comparison='reference', add_potentials=False, log_transform=True)))
                    var.process_functions.append((var.calc_pN_pS, dict(grouping='site', comparison='consensus', add_potentials=False, log_transform=True)))
                    var.process_functions.append((var.calc_pN_pS, dict(grouping='site', comparison='popular_consensus', add_potentials=False, log_transform=True)))

                # we convert counts to frequencies so high-covered samples do not skew averaging
                # across samples
                var.process_functions.append((var.convert_counts_to_frequencies, {}))
                var.process(exit_if_data_empty=False)

            gene_var[engine]['var_object'] = var
            self.run.info_single('%s for gene %s are loaded' % ('SAAVs' if engine == 'AA' else 'SCVs', gene_callers_id),
                                 progress = self.progress)

        self.variability_storage[gene_callers_id] = gene_var
        for engine in self.available_engines:
            self.variability_storage[gene_callers_id][engine]['column_info'] = self.get_column_info(gene_callers_id, engine)

        self.progress.end()


    def get_initial_data(self):
        output = {}
        output['available_gene_callers_ids'] = list(self.available_genes)
        output['available_engines'] = self.available_engines
        output['sample_groups'] = self.sample_groups
        return output


    def get_structure(self, gene_callers_id):
        """Calls self.profile_gene_variability_data, then returns `summary`

        Returns
        =======
        summary : dict
            summary has the following keys: {'pdb_content', 'residue_info', 'histograms'}
        """

        # Put the gene in self.variability_storage if it isn't already there
        self.profile_gene_variability_data(gene_callers_id)

        # Build summary
        summary = {}

        structure_db = structureops.StructureDatabase(self.structure_db_path, 'none', ignore_hash=True)
        summary['pdb_content'] = structure_db.get_pdb_content(gene_callers_id)
        structure_db.disconnect()

        residue_info, residue_info_types = self.get_residue_info_for_gene(gene_callers_id)
        summary['residue_info'] = residue_info.to_json(orient='index')
        summary['residue_info_types'] = residue_info_types.to_json(orient='index')

        summary['histograms'] = {}
        for engine in self.available_engines:
            summary['histograms'][engine] = self.get_histograms(self.variability_storage[gene_callers_id][engine]['var_object'],
                                                                self.variability_storage[gene_callers_id][engine]['column_info'])

        return summary


    def get_residue_info_for_gene(self, gene_callers_id):
        """Grab the residue info from both the structure database and the contigs database

        This grabs residue info from both `residue_info` in the structure database as well as
        any potential data present in `amino_acid_additional_data` in the contigs database.
        """

        # Get residue info from `residue_info`
        structure_db = structureops.StructureDatabase(self.structure_db_path, 'none', ignore_hash=True)
        from_structure_db = structure_db.get_residue_info_for_gene(gene_callers_id)
        structure_db.disconnect()

        # Get residue info from `amino_acid_additional_data`
        from_contigs_db = self.amino_acid_additional_data.get_gene_dataframe(gene_callers_id)

        residue_info = from_structure_db.merge(from_contigs_db, 'left', on='codon_order_in_gene')
        residue_info.drop(['entry_id', 'corresponding_gene_call', 'codon_order_in_gene'], axis=1, inplace=True)

        residue_info_types = residue_info.drop('codon_number', axis=1).aggregate([numpy.min, numpy.max, numpy.dtype]).T
        residue_info_types['dtype'] = residue_info_types['dtype'].astype(str)

        return residue_info, residue_info_types


    def get_variability(self, options):
        selected_engine = options['engine']
        gene_callers_id = int(options['gene_callers_id'])
        column_info = self.variability_storage[gene_callers_id][selected_engine]['column_info']

        self.progress.new('Filtering %s' % ('SCVs' if selected_engine == 'CDN' else 'SAAVs'), discard_previous_if_exists=True)
        output = {}

        for group in options['groups']:
            self.progress.update('Group `%s`...' % group)
            samples_in_group = options['groups'][group]

            # var becomes a filtered subset of variability_storage. it is a deepcopy so that filtering is not irreversible
            var = copy.deepcopy(self.variability_storage[gene_callers_id][selected_engine]['var_object'])
            self.compute_merged_variability(var, column_info, samples_in_group)
            self.wrangle_merged_variability(var)

            var.filter_batch_parameters(options['filter_params'], name='merged')

            output[group] = {
                'data': var.merged.to_json(orient='index'),
                'entries_after_filtering': var.merged.shape[0]
            }

        self.progress.end()
        return output


    def store_variability(self, options):
        """Store the variability currently displayed on the user's screen in a user-specified path"""

        try:
            filesnpaths.is_output_file_writable(options['path'], ok_if_exists=False)

            self.progress.new('Reporting', discard_previous_if_exists=True)
            self.progress.update('Currently viewed variants')

            selected_engine = options['engine']
            gene_callers_id = int(options['gene_callers_id'])

            # prior to filtering, var starts as a copy from variability_storage
            var = copy.deepcopy(self.variability_storage[gene_callers_id][selected_engine]['var_object'])
            var.convert_frequencies_to_counts()

            # add selected samples into the filter parameters dictionary
            selected_samples = []
            for group in options['groups']:
                selected_samples.extend(options['groups'][group])
            options["filter_params"]['sample_id'] = {'sample_ids_of_interest': set(selected_samples)}

            var.filter_batch_parameters(options['filter_params'])

            var.output_file_path = options['path']
            var.report()

            self.progress.end()

            return {'success': 'Success! Saved to %s' % options['path']}

        except Exception as e:
            return {'failure': 'Error: %s' % e}


    def store_structure_as_pdb(self, options):
        """Store the structure currently displayed on the user's screen in a user-specified path"""

        try:
            filesnpaths.is_output_file_writable(options['path'], ok_if_exists=False)

            self.progress.new('Storing', discard_previous_if_exists=True)
            self.progress.update('Currently viewed structure')

            gene_callers_id = int(options['gene_callers_id'])
            structure_db = structureops.StructureDatabase(self.structure_db_path, 'none', ignore_hash=True)
            structure_db.export_pdb_content(gene_callers_id, options['path'], ok_if_exists=False)

            self.progress.end()

            return {'success': 'Success! Saved to %s' % options['path']}

        except Exception as e:
            return {'failure': 'Error: %s' % e}


    def get_histograms(self, var_object, column_info_list):
        numpy.seterr(invalid='ignore')
        self.progress.new('Generating %s histograms' % ("SAAV" if var_object.engine else "SCV"))

        # subset info list to columns that occur in var_object.data
        column_info_list = [info for info in column_info_list if info["name"] in var_object.data.columns]

        histograms = {}
        for column_info in column_info_list:
            column = column_info["name"]
            self.progress.update("`%s`..." % column_info["name"])

            histograms[column] = {}

            if column_info["as_filter"] in ["slider"]:
                # make a number histogram
                histogram_args = {
                    'range': (column_info["min"], column_info["max"]),
                    'bins': 30,
                }
                values, bins = var_object.get_histogram(column, fix_offset=False, **histogram_args)

            elif column_info["as_filter"] in ["checkbox"]:
                # make a bar chart (categorical)
                category_counts_df = var_object.data[column].value_counts().reset_index()
                values = category_counts_df[column]
                bins = category_counts_df["index"]

            elif not column_info["as_filter"]:
                continue

            else:
                self.progress.end()
                raise ConfigError("StructureInteractive :: %s is not a recognizable controller type" % (column_info["as_filter"]))

            histograms[column]['counts'] = values.tolist()
            histograms[column]['bins'] = bins.tolist()

        self.progress.end()
        return histograms


    def compute_merged_variability(self, var, column_info, sample_group_to_merge):
        """
           var : class
               instance that inherits VariabilitySuper
           sample_names : list-like
               sample_ids to merge

           This function merges rows in self.data that share the same unique position identifier.
           For example if you have samples s01, s02, and s03, each with an entry with unique
           position identifier = 1234, this function will return dataframe containing only one entry
           with unique position identifier = 1234, with entries in this row being the mean of the
           three previous rows for columns that have numerical values and the most common value for
           columns that have categorical values. For example,

           unique_pos  sec_struc  dfc
           1234        H          1.0       becomes       unique_pos  sec_struc  dfc
           1234        H          0.8      ========>      1234        H          0.9
           1234        B          0.9
        """
        var.merged = copy.deepcopy(var.data)
        if var.merged.empty:
            return

        var.filter_data(name = 'merged', criterion = 'sample_id', subset_filter = sample_group_to_merge)

        columns_and_types_dict = {e['name']: e['data_type'] for e in column_info if not e.get('merged_only')}

        # all statistical measures
        generic_operation_dictionary = {'text': [
                                            ('',                  lambda x: x.mode().iloc[0] if not x.mode().empty else x.iloc[0]), # most common gets no suffix
                                            ('_most_common_freq', lambda x: x.value_counts().iloc[0] / x.count() if not x.value_counts().empty else 0),
                                             ],
                                        'float': [
                                            ('',                  lambda x: x.mean()), # mean gets no suffix
                                            ('_std',              lambda x: x.std() if len(x) > 1 else 0.0),
                                            #('_mini',             lambda x: x.min()),
                                            #('_maxi',             lambda x: x.max()),
                                            #('_median',           lambda x: x.median()),
                                            #('_percentile_25',    lambda x: x.quantile(0.25)),
                                            #('_percentile_50',    lambda x: x.quantile(0.50)),
                                            #('_percentile_75',    lambda x: x.quantile(0.75)),
                                             ],
                                        'integer': [
                                            ('',                  lambda x: x.mean()),
                                            ('_std',              lambda x: x.std()),
                                            #('_mini',             lambda x: x.min()),
                                            #('_maxi',             lambda x: x.max()),
                                            #('_median',           lambda x: x.median()),
                                            #('_percentile_25',    lambda x: x.quantile(0.25)),
                                            #('_percentile_50',    lambda x: x.quantile(0.50)),
                                            #('_percentile_75',    lambda x: x.quantile(0.75)),
                                             ],
                                       }

        # e.g. {'dfc':[('dfc_min', mini), ...], 'sec_struc':[('sec_struct_most_common', most_common), ...]}
        column_operations = {k: [(k + x, y) for x, y in generic_operation_dictionary[v]] for k, v in columns_and_types_dict.items()}

        # ad-hoc operations
        # occurrence: the number of samples that contained a given SAAV in the sample group
        # prevalence: the frequency of samples that contained a given SAAV in the sample group
        var.merged['occurrence'] = 0 # initialize column
        var.merged['prevalence'] = 0 # initialize column
        column_operations.update({'occurrence': [('occurrence', lambda x: x.count())],
                                  'prevalence': [('prevalence', lambda x: x.count() / len(sample_group_to_merge))],
                                  'sample_id':  [('sample_ids', lambda x: ", ".join(list(x.unique())))]})

        var.merged = var.merged.groupby('unique_pos_identifier').agg(column_operations)
        var.merged.columns = var.merged.columns.droplevel()


    def wrangle_merged_variability(self, var):
        if var.merged.empty:
            return
        # determine top n items and frequencies per variant
        self.num_reported_frequencies = 5
        for i, s in var.merged[var.items].iterrows():
            top_freqs = s.sort_values(ascending=False).iloc[:self.num_reported_frequencies]
            for x in range(self.num_reported_frequencies):
                var.merged.loc[i, str(x)+'_item'] = top_freqs.index[x]
                var.merged.loc[i, str(x)+'_freq'] = top_freqs.iloc[x]
                if var.engine == 'CDN':
                    var.merged.loc[i, str(x)+'_item_AA'] = constants.codon_to_AA[top_freqs.index[x]]

        # delete all over freq data
        columns_to_drop = []
        columns = var.merged.columns
        for item in var.items:
            columns_to_drop.extend([x for x in columns if x.startswith(item)])
        var.merged.drop(labels=columns_to_drop, axis=1, inplace=True)


class MetabolismInteractive():
    def __init__(self, args, run=run, progress=progress):
        self.mode = "metabolism"

        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db') #TODO delete if we don't need this

        self.estimator = kegg.KeggMetabolismEstimator(args)


    def get_metabolism_data(self):
        return self.estimator.get_metabolism_data_for_visualization()


class PangraphInteractive(PanSuperclass):
    def __init__(self, args, run=run, progress=progress):
        self.mode = "pangraph"

        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.pan_graph_json_path = A('pan_graph_json')
        self.pan_graph_summary_path = A('pan_graph_summary')

        if self.pan_graph_json_path:
            filesnpaths.is_file_json_formatted(self.pan_graph_json_path)

        if self.pan_graph_summary_path:
            filesnpaths.is_file_json_formatted(self.pan_graph_summary_path)

        if not self.pan_graph_json_path:
            raise ConfigError("Unfortunately you can only use this program with the `--pan-graph-json` parameter.")

        if A('pan_db_path'):
            PanSuperclass.__init__(self, self.args)
        else:
            self.run.warning("Since you did not provide anvi'o with a pan-db file, this program will initiate "
                             "the pangenome graph interface with the JSON file you have provided",
                             header="NO PAN DB, BUT ALL GOOD 👍", lc="yellow")

        self.pan_graph_json = self.get_pangraph_json()
        self.pan_graph_summary = self.get_pangraph_summary()


    def get_pangraph_json(self):
        """A function to 'get' pangraph JSON data from wherever appropriate.

        Currently this function returns the user-provided JSON data, but it
        will take into consideration other parameters and pan-db to determine
        which resource is the most appropriate to retrieve the pangraph data.
        """

        if self.pan_graph_json_path:
            json_data = json.load(open(self.pan_graph_json_path))
        else:
            # FIXME: this is where we will get things from the pan-db, but it
            #        is not yet implemented
            raise ConfigError("Not implemented.")

        # make sure the JSON data we have access is compatible with the codebase
        # we are running.
        version_we_have = str(json_data['infos']['meta']['version']) if 'version' in json_data['infos']['meta'] else None
        version_we_want = str(t.pangraph_json_version)

        if version_we_have == None:
            raise ConfigError("Bad news: the pan-graph data you have is outdated :/ You will have to re-run "
                              "the program `anvi-pan-graph` on these data to generate an up-to-date "
                              "pan-graph artifact that your current codebase can work with.")

        if version_we_have != version_we_want:
            raise ConfigError(f"We have a problem here. The data for your pan-graph is at version {version_we_have}. "
                              f"Yet the anvi'o codebase you have here can only work with pan-graph version {version_we_want}. "
                              f"The best solution to address this mismatch is to re-run the program `anvi-pan-graph` "
                              f"on your dataset.")

        # all good with the version stuff:
        return json_data


    def get_pangraph_summary(self):
        if self.pan_graph_summary_path:
            summary_data = json.load(open(self.pan_graph_summary_path))

            return summary_data


class ContigsInteractive():
    def __init__(self, args, run=run, progress=progress):
        self.mode = 'contigs'

        self.args = args
        self.run = run
        self.progress = progress

        self.contigs_stats = {}

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.input_contig_db_paths = A('input')

        if not len(self.input_contig_db_paths):
            raise ConfigError("ContigsInteractive should be inherited with an args object with a valid `input` "
                              "member. Not like the way you tried it with no input paths whatsoever :/")

        for contig_db_path in self.args.input:
            self.contigs_stats[contig_db_path] = summarizer.ContigSummarizer(contig_db_path).get_summary_dict_for_assembly()

        self.tables = {}
        self.generate_tables()


    def generate_tables(self):
        # let's keep track of all keys we will need to access later from the interface. if
        # we don't do this, non-standard keys (such as 'Gene caller (prodigal)' becomes very
        # inaccessable when we need to access to it the way we access to 'N50' or 'Contig
        # Lengths'):
        self.human_readable_keys = []

        self.tables['header'] = [c['project_name'] for c in self.contigs_stats.values()]

        basic_stats = []

        ##
        ##  Table for basic stats
        ##
        self.progress.new('Generating stats tables')

        self.progress.update('Basic stats ...')
        basic_stats.append(['Total Length'] + [c['total_length'] for c in self.contigs_stats.values()])
        basic_stats.append(['Num Contigs'] + [c['num_contigs'] for c in self.contigs_stats.values()])

        self.progress.update('Number of contigs ...')
        contig_lengths_for_all = [c['contig_lengths'] for c in self.contigs_stats.values()]
        X = lambda n: [len([count for count in contig_lengths_for_one if count >= n]) for contig_lengths_for_one in contig_lengths_for_all]
        basic_stats.append(['Num Contigs > 100 kb'] + X(100000))
        basic_stats.append(['Num Contigs > 50 kb'] + X(50000))
        basic_stats.append(['Num Contigs > 20 kb'] + X(20000))
        basic_stats.append(['Num Contigs > 10 kb'] + X(10000))
        basic_stats.append(['Num Contigs > 5 kb'] + X(5000))
        basic_stats.append(['Num Contigs > 2.5 kb'] + X(2500))

        self.progress.update('Contig lengths ...')
        contig_lengths_for_all = [c['contig_lengths'] for c in self.contigs_stats.values()]
        MAX_L = lambda: [max(lengths) for lengths in contig_lengths_for_all]
        MIN_L = lambda: [min(lengths) for lengths in contig_lengths_for_all]
        basic_stats.append(['Longest Contig'] + MAX_L())
        basic_stats.append(['Shortest Contig'] + MIN_L())

        self.progress.update('Number of genes ...')
        contig_lengths_for_all = [c['contig_lengths'] for c in self.contigs_stats.values()]
        basic_stats.append(['Num Genes (' + constants.default_gene_caller + ')'] + [c['num_genes'] for c in self.contigs_stats.values()])

        self.progress.update('N/L values ...')
        n_values = [c['n_values'] for c in self.contigs_stats.values()]
        N = lambda n: [n_value[n]['num_contigs'] for n_value in n_values]
        L = lambda n: [n_value[n]['length'] for n_value in n_values]
        basic_stats.append(['L50'] + N(49))
        basic_stats.append(['L75'] + N(74))
        basic_stats.append(['L90'] + N(89))
        basic_stats.append(['N50'] + L(49))
        basic_stats.append(['N75'] + L(74))
        basic_stats.append(['N90'] + L(89))

        self.tables['basic_stats'] = basic_stats
        self.human_readable_keys.extend([e[0] for e in basic_stats])

        ##
        ##  Table for hmm hits
        ##
        self.progress.update('HMMs summary ...')
        all_hmm_sources = set()
        for c in self.contigs_stats.values():
            for source in c['gene_hit_counts_per_hmm_source'].keys():
                all_hmm_sources.add(source)

        hmm_table = []
        for source in sorted(all_hmm_sources):
            line = [source]
            for c in self.contigs_stats.values():
                if source in c['gene_hit_counts_per_hmm_source']:
                    line.append(sum(c['gene_hit_counts_per_hmm_source'][source].values()))
                else:
                    line.append('n/a')

            self.human_readable_keys.append(line[0])
            hmm_table.append(line)

        self.tables['hmm'] = hmm_table

        ##
        ##  Table for SCG genome prediction
        ##
        self.progress.update('Num genome prediction ...')
        source_to_domain = {}
        all_scg_sources = set()
        for c in self.contigs_stats.values():
            for source in c['num_genomes_per_SCG_source_dict'].keys():
                all_scg_sources.add(source)

                if source not in source_to_domain:
                    source_to_domain[source] = c['num_genomes_per_SCG_source_dict'][source]['domain']

        scg_table = []
        for source in all_scg_sources:
            line = [source_to_domain[source] + ' (' + source + ')']
            for c in self.contigs_stats.values():
                if source in c['num_genomes_per_SCG_source_dict']:
                    line.append(c['num_genomes_per_SCG_source_dict'][source]['num_genomes'])
                else:
                    line.append('n/a')

            self.human_readable_keys.append(line[0])
            scg_table.append(line)

        self.tables['scg'] = scg_table

        self.progress.end()
