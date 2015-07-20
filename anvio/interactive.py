# -*- coding: utf-8
"""The module that curates data for the interactive interface"""

import os
import sys

import anvio
import anvio.utils as utils
import anvio.dictio as dictio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.dbops import ProfileSuperclass, AnnotationSuperclass, TablesForStates, is_annotation_and_profile_dbs_compatible
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


class InputHandler(ProfileSuperclass, AnnotationSuperclass):
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args, external_clustering = None):
        self.args = args
        self.views = {}
        self.states_table = None
        self.p_meta = {}
        self.title = 'Unknown Project'

        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        self.state = A('state')
        self.split_hmm_layers = A('split_hmm_layers')
        self.additional_metadata_path = A('additional_metadata')
        self.additional_view_path = A('additional_view')
        self.profile_db_path = A('profile_db')
        self.annotation_db_path = A('annotation_db')
        self.view = A('view')
        self.fasta_file = A('fasta_file')
        self.metadata = A('metadata')
        self.tree = A('tree')
        self.title = A('title')
        self.summary_index = A('summary_index')
        self.output_dir = A('output_dir')
        self.show_views = A('show_views')
        self.skip_check_names = A('skip_check_names')

        self.split_names_ordered = None
        self.splits_summary_index = {}
        self.additional_metadata = None
        self.external_clustering = external_clustering

        self.collections = ccollections.Collections()

        AnnotationSuperclass.__init__(self, self.args)

        if self.annotation_db_path:
            self.completeness = completeness.Completeness(self.annotation_db_path)
            self.collections.populate_sources_dict(self.annotation_db_path, anvio.__annotation__version__)
        else:
            self.completeness = None

        if self.annotation_db_path and self.profile_db_path:
            # make sure we are not dealing with apples and oranges here.
            is_annotation_and_profile_dbs_compatible(self.annotation_db_path, self.profile_db_path)

        self.P = lambda x: os.path.join(self.p_meta['output_dir'], x)
        self.cwd = os.getcwd()

        # here is where the big deal stuff takes place:
        if self.profile_db_path:
            if not self.annotation_db_path:
                raise ConfigError, "Anvi'o needs the annotation database to make sense of this run."

            ProfileSuperclass.__init__(self, args)

            # this is a weird place to do it, but we are going to ask AnnotationSuperclass function to load
            # all the split sequences since only now we know the mun_contig_length that was used to profile
            # this stuff
            self.init_split_sequences(self.p_meta['min_contig_length'])

            self.collections.populate_sources_dict(self.profile_db_path, anvio.__profile__version__)

            self.load_from_profile_database(args)
        else:
            self.load_from_files(args)

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

        # if there are any HMM search results in the annotation database other than 'singlecopy' sources,
        # we would like to visualize them as additional layers. following function is inherited from
        # Annotation DB superclass and will fill self.hmm_searches_dict if appropriate data is found in
        # search tables:
        self.init_non_singlecopy_gene_hmm_sources(self.split_names_ordered, return_each_gene_as_a_layer = self.split_hmm_layers)

        if self.additional_metadata_path:
            filesnpaths.is_file_tab_delimited(self.additional_metadata_path)
            self.additional_metadata = self.additional_metadata_path

        self.check_names_consistency()
        self.convert_metadata_into_json()


    def load_from_files(self, args):
        if (not self.fasta_file) or (not self.metadata) or (not self.tree) or (not self.output_dir):
            raise ConfigError, "If you do not have a RUNINFO dict, you must declare each of\
                                           '-f', '-m', '-t' and '-o' parameters. Please see '--help' for\
                                           more detailed information on them."

        if self.view:
            raise ConfigError, "You can't use '-v' parameter when this program is not called with a RUNINFO.cp"

        if self.show_views:
            raise ConfigError, "Sorry, there are no views to show when there is no RUNINFO.cp :/"

        metadata_path = os.path.abspath(self.metadata)
        self.p_meta['splits_fasta'] = os.path.abspath(self.fasta_file)
        self.p_meta['output_dir'] = os.path.abspath(self.output_dir)
        self.p_meta['views'] = {}
        self.p_meta['default_view'] = 'single'
        self.p_meta['default_clustering'] = 'default'
        self.p_meta['available_clusterings'] = ['default']
        self.p_meta['clusterings'] = {'default': {'newick': open(os.path.abspath(self.tree)).read()}}

        self.default_view = self.p_meta['default_view']

        if self.summary_index:
            self.p_meta['profile_summary_index'] = os.path.abspath(self.summary_index)
            self.splits_summary_index = dictio.read_serialized_object(self.p_meta['profile_summary_index'])

        # sanity of the metadata
        filesnpaths.is_file_tab_delimited(metadata_path)
        metadata_columns = utils.get_columns_of_TAB_delim_file(metadata_path, include_first_column=True)
        if not metadata_columns[0] == "contig":
            raise ConfigError, "The first row of the first column of the metadata file must\
                                      say 'contig', which is not the case for your metadata file\
                                      ('%s'). Please make sure this is a properly formatted metadata\
                                      file." % (metadata_path)

        # store metadata as view:
        self.views[self.default_view] = {'header': metadata_columns[1:],
                                         'dict': utils.get_TAB_delimited_file_as_dictionary(metadata_path)}
        self.split_names_ordered = self.views[self.default_view]['dict'].keys()

        filesnpaths.is_file_fasta_formatted(self.p_meta['splits_fasta'])
        self.split_sequences = utils.get_FASTA_file_as_dictionary(self.p_meta['splits_fasta'])

        # setup a mock splits_basic_info dict
        self.splits_basic_info = {}
        for split_id in self.split_names_ordered:
            self.splits_basic_info[split_id] = {'length': len(self.split_sequences[split_id]),
                                                'gc_content': utils.get_GC_content_for_sequence(self.split_sequences[split_id])}

        # reminder: this is being stored in the output dir provided as a commandline parameter:
        self.p_meta['self_path'] = os.path.join(self.p_meta['output_dir'], 'RUNINFO.cp')

        if self.title:
            self.title = self.title

        filesnpaths.gen_output_directory(self.p_meta['output_dir'])


    def load_from_profile_database(self, args):
        if self.p_meta['version'] != anvio.__profile__version__:
            raise ConfigError, "The profile database has a version number that differs from the version that is valid\
                                for this codebase (the profile database is at '%s', and the codebase is at '%s'). Very\
                                unfortunately, you need to re-profile and re-merge this project using the current anvi'o :("

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

        # is summary being overwritten?
        if self.summary_index:
            run.info('Warning', "The default summary index in RUNINFO is being overriden by '%s'." % self.summary_index)
            self.p_meta['profile_summary_index'] = os.path.abspath(self.summary_index)

        if os.path.exists(self.P('SUMMARY.cp')):
            self.splits_summary_index = dictio.read_serialized_object(self.P('SUMMARY.cp'))
        else:
            self.splits_summary_index = None
            run.warning("SUMMARY.cp is missing for your run. Anvi'o will continue working (well, at least\
                         it will attempt to do it), but things may behave badly with the absence of\
                         SUMMARY.cp (first and foremost, you will not be able to inspect individual\
                         contigs through any of the interactive interfaces). Please investigate it\
                         if you were not expecting this.")

        # set title
        if self.title:
            self.title = self.title + ' (%s)' % self.default_view
        else:
            self.title = self.p_meta['sample_id'] + ' (%s)' % self.default_view


    def check_names_consistency(self):
        if self.skip_check_names:
            return

        splits_in_tree = set(self.split_names_ordered)
        splits_in_metadata = set(self.views[self.default_view]['dict'].keys())
        splits_in_database = set(self.split_sequences)

        splits_in_tree_but_not_in_metadata = splits_in_tree - splits_in_metadata
        splits_in_tree_but_not_in_database = splits_in_tree - splits_in_database

        if splits_in_tree_but_not_in_metadata:
            raise ConfigError, 'Some split names found in your tree are missing in your metadata. Hard to\
                                know what cuased this, but here is a couple of them that: %s'\
                                    % ', '.join(splits_in_tree_but_not_in_metadata[0:5])

        if splits_in_tree_but_not_in_database:
            raise ConfigError, 'Some split names found in your tree are missing from your database. Hard to\
                                know why is this the case, but here is a couple of them: %s'\
                                    % ', '.join(splits_in_tree_but_not_in_database[0:5])

        if self.additional_metadata_path:
            splits_in_additional_metadata = set(sorted([l.split('\t')[0] for l in open(self.additional_metadata_path).readlines()[1:]]))
            splits_only_in_additional_metadata = []
            for split_name in splits_in_additional_metadata:
                if split_name not in splits_in_tree:
                    splits_only_in_additional_metadata.append(split_name)
            if len(splits_only_in_additional_metadata):
                one_example = splits_only_in_additional_metadata[-1]
                num_all = len(splits_only_in_additional_metadata)
                run.warning("Some of the contigs in your addtional metadata file does not\
                            appear to be in anywhere else. Additional metadata file is not\
                            required to list all contigs (which means, there may be contigs\
                            in the database that are not in the additional metadata file),\
                            however, finding contigs that are only in the additional metadata\
                            file usually means trouble. anvio will continue, but please\
                            go back and check your files if you think there may be something\
                            wrong. Here is a random contig name that was only in the\
                            metadata file: '%s'. And there were %d of them in total. You\
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


    def convert_metadata_into_json(self):
        '''This function's name must change to something more meaningful.'''

        additional_dict, additional_headers = None, []
        if self.additional_metadata_path:
            additional_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_metadata_path)
            additional_headers = utils.get_columns_of_TAB_delim_file(self.additional_metadata_path)

        for view in self.views:
            # here we will populate runinfo['views'] with json objects.
            view_dict = self.views[view]['dict']
            view_headers = self.views[view]['header']

            json_object = []

            # (1) set the header line with the first entry:
            json_header = ['contigs']

            # (2) then add annotation db stuff, if exists
            if len(self.genes_in_splits_summary_dict):
                json_header.extend(self.genes_in_splits_summary_headers[1:])

            # (3) then add length and GC content
            basic_info_headers = ['length', 'gc_content']
            json_header.extend(basic_info_headers)

            # (4) then add the view!
            json_header.extend(view_headers)

            # (5) then add 'additional' headers as the outer ring:
            if additional_headers:
                json_header.extend(additional_headers)

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

                # (5) adding additional data
                json_entry.extend([additional_dict[split_name][header] if additional_dict.has_key(split_name) else None for header in additional_headers])

                # (6) adding hmm stuff
                if self.hmm_searches_dict:
                    if self.split_hmm_layers:
                        json_entry.extend([self.hmm_searches_dict[split_name][header] if self.hmm_searches_dict.has_key(split_name) else None for header in [tpl[0] for tpl in self.hmm_searches_header]])
                    else:
                        json_entry.extend([len(self.hmm_searches_dict[split_name][header]) if self.hmm_searches_dict.has_key(split_name) else 0 for header in [tpl[0] for tpl in self.hmm_searches_header]])

                # (7) send it along!
                json_object.append(json_entry)

            self.views[view] = json_object


    def end(self):
        # FIXME: remove temp files and stuff
        pass
