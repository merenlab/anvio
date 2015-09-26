# -*- coding: utf-8
"""The module that curates data for the interactive interface"""

import os
import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.dbops import ProfileSuperclass, ContigsSuperclass, SamplesInformationDatabase, TablesForStates, ProfileDatabase
from anvio.dbops import is_profile_db_and_contigs_db_compatible, is_profile_db_and_samples_db_compatible
from anvio.errors import ConfigError

with terminal.SuppressAllOutput():
    from ete2 import Tree


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


class InputHandler(ProfileSuperclass, ContigsSuperclass):
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args, external_clustering = None):
        self.args = args
        self.views = {}
        self.states_table = None
        self.p_meta = {}
        self.title = 'Unknown Project'

        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.manual_mode = A('manual_mode')
        self.state = A('state')
        self.split_hmm_layers = A('split_hmm_layers')
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
        self.skip_check_names = A('skip_check_names')

        self.split_names_ordered = None
        self.additional_layers = None

        self.samples_information_dict = {}
        self.samples_order_dict = {}


        self.external_clustering = external_clustering

        self.collections = ccollections.Collections()

        ContigsSuperclass.__init__(self, self.args)

        if self.samples_information_db_path:
            samples_information_db = SamplesInformationDatabase(self.samples_information_db_path)
            self.samples_information_dict, self.samples_order_dict = samples_information_db.get_samples_information_and_order_dicts()

        if self.contigs_db_path:
            self.completeness = completeness.Completeness(self.contigs_db_path)
            self.collections.populate_sources_dict(self.contigs_db_path, anvio.__contigs__version__)
        else:
            self.completeness = None

        # make sure we are not dealing with apples and oranges here.
        if self.contigs_db_path and self.profile_db_path:
            is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        # make sure the samples information database, if there is one, is in fact compatible with the profile database
        if self.profile_db_path and self.samples_information_db_path:
            is_profile_db_and_samples_db_compatible(self.profile_db_path, self.samples_information_db_path)

        self.P = lambda x: os.path.join(self.p_meta['output_dir'], x)
        self.cwd = os.getcwd()

        # here is where the big deal stuff takes place:
        if self.manual_mode:
            self.load_from_user_files(args)
        else:
            self.load_from_anvio_files(args)

        if self.external_clustering:
            self.p_meta['clusterings'] = self.clusterings = self.external_clustering['clusterings']
            self.p_meta['available_clusterings'] = self.clusterings.keys()
            self.p_meta['default_clustering'] = self.external_clustering['default_clustering']

        if not self.p_meta['clusterings']:
            if self.p_meta['merged']:
                raise ConfigError, "This merged profile database does not seem to have any hierarchical clustering\
                                    that is required by the interactive interface. It may have been generated\
                                    by anvi-merge with `--skip-hierarchical-clustering` flag, or hierarchical\
                                    clustering step may have been skipped automatically by the platform. Please\
                                    read the help menu for anvi-merge, and/or refer to the tutorial: \
                                    http://merenlab.org/2015/05/01/anvio-tutorial/#clustering-during-merging"
            else:
                raise ConfigError, "This single profile database does not seem to have any hierarchical clustering\
                                    that is required by the interactive interface. You must use `--cluster-contigs`\
                                    flag for single profiles to access to this functionality. Please read the help\
                                    menu for anvi-profile, and/or refer to the tutorial."

        tree = Tree(self.p_meta['clusterings'][self.p_meta['default_clustering']]['newick'], format = 1)

        # self.split_names_ordered is going to be the 'master' names list. everything else is going to
        # need to match these names:
        self.split_names_ordered = [n.name for n in tree.get_leaves()]

        # now we knot what splits we are interested in (self.split_names_ordered), we can get rid of all the
        # unnecessary splits stored in views dicts.
        self.prune_view_dicts()

        # if there are any HMM search results in the contigs database other than 'singlecopy' sources,
        # we would like to visualize them as additional layers. following function is inherited from
        # Contigs DB superclass and will fill self.hmm_searches_dict if appropriate data is found in
        # search tables:
        self.init_non_singlecopy_gene_hmm_sources(self.split_names_ordered, return_each_gene_as_a_layer = self.split_hmm_layers)

        if self.additional_layers_path:
            filesnpaths.is_file_tab_delimited(self.additional_layers_path)
            self.additional_layers = self.additional_layers_path

        self.check_names_consistency()
        self.convert_view_data_into_json()


    def load_from_user_files(self, args):
        if self.contigs_db_path:
            raise ConfigError, "When you want to use the interactive interface in an ad hoc manner, you must\
                                not use a contigs database."

        if not self.profile_db_path:
            raise ConfigError, "Even when you want to use the interactive interface in an ad hoc manner by\
                                using the '--manual-mode' flag, you still need to declare a profile database.\
                                The profile database in this mode only used to read or store the 'state' of\
                                the display for visualization purposes. You DO NOT need to point to an already\
                                existing database, as anvi'o will generate an empty one for your if there is no\
                                profile database."

        if (not self.fasta_file) or (not self.view_data_path) or (not self.tree):
            raise ConfigError, "When you are running the interactive interface in manual mode, you must declare\
                                each of '-f', '-d', and '-t' parameters. Please see the help menu for more info."

        if self.view:
            raise ConfigError, "You can't use '--view' parameter when you are running the interactive interface\
                                in manual mode"

        if self.show_views:
            raise ConfigError, "Sorry, there are no views to show in manual mode :/"


        # create a new, empty profile database for ad hoc operations
        if not os.path.exists(self.profile_db_path):
            profile_db = ProfileDatabase(self.profile_db_path)
            profile_db.create({'db_type': 'profile', 'contigs_db_hash': None})

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path, anvio.__profile__version__)

        # also populate collections, if there are any
        self.collections.populate_sources_dict(self.profile_db_path, anvio.__profile__version__)

        view_data_path = os.path.abspath(self.view_data_path)
        self.p_meta['splits_fasta'] = os.path.abspath(self.fasta_file)
        self.p_meta['output_dir'] = None
        self.p_meta['views'] = {}
        self.p_meta['default_view'] = 'single'
        self.p_meta['default_clustering'] = 'default'
        self.p_meta['available_clusterings'] = ['default']
        self.p_meta['clusterings'] = {'default': {'newick': open(os.path.abspath(self.tree)).read()}}

        self.default_view = self.p_meta['default_view']

        # sanity of the view data
        filesnpaths.is_file_tab_delimited(view_data_path)
        view_data_columns = utils.get_columns_of_TAB_delim_file(view_data_path, include_first_column=True)
        if not view_data_columns[0] == "contig":
            raise ConfigError, "The first row of the first column of the view data file must\
                                      say 'contig', which is not the case for your view data file\
                                      ('%s'). Please make sure this is a properly formatted view data\
                                      file." % (view_data_path)

        # store view data as view:
        self.views[self.default_view] = {'header': view_data_columns[1:],
                                         'dict': utils.get_TAB_delimited_file_as_dictionary(view_data_path)}
        self.split_names_ordered = self.views[self.default_view]['dict'].keys()

        filesnpaths.is_file_fasta_formatted(self.p_meta['splits_fasta'])
        self.split_sequences = utils.get_FASTA_file_as_dictionary(self.p_meta['splits_fasta'])

        # setup a mock splits_basic_info dict
        self.splits_basic_info = {}
        for split_id in self.split_names_ordered:
            self.splits_basic_info[split_id] = {'length': len(self.split_sequences[split_id]),
                                                'gc_content': utils.get_GC_content_for_sequence(self.split_sequences[split_id])}

        if self.title:
            self.title = self.title


    def load_from_anvio_files(self, args):
        if not self.contigs_db_path:
            raise ConfigError, "Anvi'o needs the contigs database to make sense of this run."

        ProfileSuperclass.__init__(self, args)

        # this is a weird place to do it, but we are going to ask ContigsSuperclass function to load
        # all the split sequences since only now we know the mun_contig_length that was used to profile
        # this stuff
        self.init_split_sequences(self.p_meta['min_contig_length'])

        self.collections.populate_sources_dict(self.profile_db_path, anvio.__profile__version__)

        self.p_meta['self_path'] = self.profile_db_path
        self.p_meta['output_dir'] = os.path.join(os.getcwd(), os.path.dirname(self.profile_db_path))

        # create an instance of states table
        self.states_table = TablesForStates(self.profile_db_path, anvio.__profile__version__)

        # load views from the profile database
        self.load_views()
        self.default_view = self.p_meta['default_view']

        # if the user wants to see available views, show them and exit.
        if self.show_views:
            run.warning('', header = 'Available views (%d)' % len(self.views), lc = 'green')
            for view in self.views:
                run.info(view,
                         'Via "%s" table' % self.views[view]['table_name'],
                         lc='crimson',
                         mc='green' if view == self.default_view else 'crimson')
            print
            sys.exit()

        # if the user has an additional view data, load it up into the self.views dict.
        if self.additional_view_path:
            filesnpaths.is_file_tab_delimited(self.additional_view_path)
            additional_view_columns = utils.get_columns_of_TAB_delim_file(self.additional_view_path)

            if not additional_view_columns[-1] == '__parent__':
                raise ConfigError, "The last column of the additional view must be '__parent__' with the proper\
                                    parent information for each split."

            column_mapping = [str] + [float] * (len(additional_view_columns) - 1) + [str]

            self.views['user_view'] = {'table_name': 'NA',
                                       'header': additional_view_columns,
                                       'dict': utils.get_TAB_delimited_file_as_dictionary(self.additional_view_path, column_mapping = column_mapping)}

        # if the user specifies a view, set it as default:
        if self.view:
            if not self.view in self.views:
                raise ConfigError, "The requested view ('%s') is not available for this run. Please see\
                                          available views by running this program with --show-views flag." % self.view

            self.default_view = self.view

        self.p_meta['clusterings'] = self.clusterings 

        if self.tree:
            entry_id = os.path.basename(self.tree).split('.')[0]
            if not self.p_meta['clusterings']:
                self.p_meta['default_clustering'] = entry_id
                self.p_meta['available_clusterings'] = [entry_id]
                self.p_meta['clusterings'] = {entry_id: {'newick': open(os.path.abspath(self.tree)).read()}}
                run.info('Additional Tree', "Splits will be organized based on '%s'." % entry_id)
            else:
                self.p_meta['clusterings'][entry_id] = {'newick': open(os.path.abspath(self.tree)).read()}
                run.info('Additional Tree', "'%s' has been added to available trees." % entry_id)

        # set title
        if self.title:
            self.title = self.title
        else:
            self.title = self.p_meta['sample_id'].replace('-', ' ').replace('_', ' ')

        # do we have auxiliary data available?
        if not self.auxiliary_data_available:
            summary_cp_available = os.path.exists(os.path.join(os.path.dirname(self.profile_db_path), 'SUMMARY.cp'))
            self.run.warning("Auxiliary data is not available; which means you will not be able to perform\
                              certain operations (i.e., the inspect menu in the interactive interface will\
                              not work, etc). %s" % ('' if not summary_cp_available else "Although, you have\
                              a SUMMARY.cp file in your work directory, which means you are working with an\
                              outdated anvi'o run. You can convert your SUMMARY.cp into an auxiliary data file\
                              by using `anvi-script-generate-auxiliary-data-from-summary-cp` script."))

    def check_names_consistency(self):
        if self.skip_check_names:
            return

        splits_in_tree = set(self.split_names_ordered)
        splits_in_view_data = set(self.views[self.default_view]['dict'].keys())
        splits_in_database = set(self.split_sequences)

        splits_in_tree_but_not_in_view_data = splits_in_tree - splits_in_view_data
        splits_in_tree_but_not_in_database = splits_in_tree - splits_in_database

        if splits_in_tree_but_not_in_view_data:
            raise ConfigError, 'Some split names found in your tree are missing in your view data. Hard to\
                                know what cuased this, but here is a couple of them that: %s'\
                                    % ', '.join(splits_in_tree_but_not_in_view_data[0:5])

        if splits_in_tree_but_not_in_database:
            raise ConfigError, 'Some split names found in your tree are missing from your database. Hard to\
                                know why is this the case, but here is a couple of them: %s'\
                                    % ', '.join(list(splits_in_tree_but_not_in_database)[0:5])

        if self.additional_layers_path:
            splits_in_additional_layers = set(sorted([l.split('\t')[0] for l in open(self.additional_layers_path).readlines()[1:]]))
            splits_only_in_additional_layers = []
            for split_name in splits_in_additional_layers:
                if split_name not in splits_in_tree:
                    splits_only_in_additional_layers.append(split_name)
            if len(splits_only_in_additional_layers):
                one_example = splits_only_in_additional_layers[-1]
                num_all = len(splits_only_in_additional_layers)
                run.warning("Some of the contigs in your addtional view data file does not\
                            appear to be in anywhere else. Additional view data file is not\
                            required to list all contigs (which means, there may be contigs\
                            in the database that are not in the additional view data file),\
                            however, finding contigs that are only in the additional view data\
                            file usually means trouble. anvio will continue, but please\
                            go back and check your files if you think there may be something\
                            wrong. Here is a random contig name that was only in the\
                            view data file: '%s'. And there were %d of them in total. You\
                            are warned!" % (one_example, num_all))


    def prune_view_dicts(self):
        self.progress.new('Pruning view dicts')
        self.progress.update('...')
        splits_in_views = set(self.views.values()[0]['dict'].keys())
        splits_to_remove = splits_in_views - set(self.split_names_ordered)

        for view in self.views:
            self.progress.update('processing view "%s"' % view)
            for split_name in splits_to_remove:
                self.views[view]['dict'].pop(split_name)

        self.progress.end()


    def convert_view_data_into_json(self):
        '''This function's name must change to something more meaningful.'''

        additional_layers_dict, additional_layer_headers = None, []
        if self.additional_layers_path:
            additional_layers_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_layers_path)
            additional_layer_headers = utils.get_columns_of_TAB_delim_file(self.additional_layers_path)

        for view in self.views:
            # here we will populate runinfo['views'] with json objects.
            view_dict = self.views[view]['dict']
            view_headers = self.views[view]['header']

            json_object = []

            # (1) set the header line with the first entry:
            json_header = ['contigs']

            # (2) then add contigs db stuff, if exists
            if len(self.genes_in_splits_summary_dict):
                json_header.extend(self.genes_in_splits_summary_headers[1:])

            # (3) then add length and GC content
            basic_info_headers = ['length', 'gc_content']
            json_header.extend(basic_info_headers)

            # (4) then add the view!
            json_header.extend(view_headers)

            # (5) then add 'additional' headers as the outer ring:
            if additional_layer_headers:
                json_header.extend(additional_layer_headers)

            # (6) finally add hmm search results
            if self.hmm_searches_header:
                json_header.extend([tpl[0] for tpl in self.hmm_searches_header])

            # (7) and finalize it (yay):
            json_object.append(json_header)

            for split_name in view_dict:
                # (1)
                json_entry = [split_name]

                # (2)
                if self.genes_in_splits_summary_dict:
                    json_entry.extend([self.genes_in_splits_summary_dict[split_name][header] for header in self.genes_in_splits_summary_headers[1:]])

                # (3)
                json_entry.extend([self.splits_basic_info[split_name][header] for header in basic_info_headers])

                # (4) adding essential data for the view
                json_entry.extend([view_dict[split_name][header] for header in view_headers])

                # (5) adding additional layers
                json_entry.extend([additional_layers_dict[split_name][header] if additional_layers_dict.has_key(split_name) else None for header in additional_layer_headers])

                # (6) adding hmm stuff
                if self.hmm_searches_dict:
                    if self.split_hmm_layers:
                        json_entry.extend([self.hmm_searches_dict[split_name][header] if self.hmm_searches_dict.has_key(split_name) else None for header in [tpl[0] for tpl in self.hmm_searches_header]])
                    else:
                        json_entry.extend([len(self.hmm_searches_dict[split_name][header]) if self.hmm_searches_dict.has_key(split_name) else 0 for header in [tpl[1] for tpl in self.hmm_searches_header]])

                # (7) send it along!
                json_object.append(json_entry)

            self.views[view] = json_object


    def end(self):
        # FIXME: remove temp files and stuff
        pass
