# -*- coding: utf-8
"""The module that curates data for the interactive interface"""

import os
import sys

import anvio.tables as t
import anvio.utils as utils
import anvio.dictio as dictio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.dbops import ProfileSuperclass, AnnotationSuperclass
from anvio.errors import ConfigError

with terminal.SuppressAllOutput():
    from ete2 import Tree


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


progress = terminal.Progress()
run = terminal.Run()


class InputHandler(ProfileSuperclass, AnnotationSuperclass):
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args):
        self.args = args
        self.views = {}
        self.runinfo = {}
        self.title = 'Unknown Project'

        self.collections = ccollections.Collections()

        AnnotationSuperclass.__init__(self, self.args)
        if self.annotation_db_path:
            self.completeness = completeness.Completeness(self.annotation_db_path)
            self.collections.populate_sources_dict(self.annotation_db_path, t.annotation_db_version)
        else:
            self.completeness = None

        self.profile_db_path = None

        self.split_names_ordered = None
        self.splits_summary_index = {}
        self.additional_metadata_path = None

        self.P = lambda x: os.path.join(self.runinfo['output_dir'], x)
        self.cwd = os.getcwd()

        self.state = args.state

        # here is where the big deal stuff takes place:
        if args.runinfo:
            if not self.annotation_db_path:
                raise ConfigError, "anvio needs the annotation database to make sense of this run."

            self.runinfo = self.read_runinfo_dict(args)

            # this is a weird place to do it, but we are going to ask AnnotationSuperclass function to load
            # all the split sequences since only now we know the mun_contig_length that was used to profile
            # this stuff
            self.init_split_sequences(self.runinfo['min_contig_length'])

            args.profile_db = self.P(self.runinfo['profile_db'])
            ProfileSuperclass.__init__(self, args)
            self.collections.populate_sources_dict(self.profile_db_path, t.profile_db_version)

            self.load_from_runinfo_dict(args)
        else:
            self.load_from_files(args)

        tree = Tree(self.runinfo['clusterings'][self.runinfo['default_clustering']]['newick'])
        self.split_names_ordered = [n.name for n in tree.get_leaves()]

        # if there are any HMM search results in the annotation database other than 'singlecopy' sources,
        # we would like to visualize them as additional layers. following function is inherited from
        # Annotation DB superclass and will fill self.hmm_searches_dict if appropriate data is found in
        # search tables:
        self.init_non_singlecopy_gene_hmm_sources(self.split_names_ordered, return_each_gene_as_a_layer = args.split_hmm_layers)

        if args.additional_metadata:
            filesnpaths.is_file_tab_delimited(args.additional_metadata)
            self.additional_metadata_path = args.additional_metadata

        self.check_names_consistency()
        self.convert_metadata_into_json()


    def load_from_files(self, args):
        if (not args.fasta_file) or (not args.metadata) or (not args.tree) or (not args.output_dir):
            raise ConfigError, "If you do not have a RUNINFO dict, you must declare each of\
                                           '-f', '-m', '-t' and '-o' parameters. Please see '--help' for\
                                           more detailed information on them."

        if args.view:
            raise ConfigError, "You can't use '-v' parameter when this program is not called with a RUNINFO.cp"

        if args.show_views:
            raise ConfigError, "Sorry, there are no views to show when there is no RUNINFO.cp :/"

        metadata_path = os.path.abspath(args.metadata)
        self.runinfo['splits_fasta'] = os.path.abspath(args.fasta_file)
        self.runinfo['output_dir'] = os.path.abspath(args.output_dir)
        self.runinfo['views'] = {}
        self.runinfo['default_view'] = 'single'
        self.runinfo['default_clustering'] = 'default'
        self.runinfo['available_clusterings'] = ['default']
        self.runinfo['clusterings'] = {'default': {'newick': open(os.path.abspath(args.tree)).read()}}

        self.default_view = self.runinfo['default_view']

        if args.summary_index:
            self.runinfo['profile_summary_index'] = os.path.abspath(args.summary_index)
            self.splits_summary_index = dictio.read_serialized_object(self.runinfo['profile_summary_index'])

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

        filesnpaths.is_file_fasta_formatted(self.runinfo['splits_fasta'])
        self.split_sequences = utils.get_FASTA_file_as_dictionary(self.runinfo['splits_fasta'])

        # setup a mock splits_basic_info dict
        self.splits_basic_info = {}
        for split_id in self.split_names_ordered:
            self.splits_basic_info[split_id] = {'length': len(self.split_sequences[split_id]),
                                                'gc_content': utils.get_GC_content_for_sequence(self.split_sequences[split_id])}

        # reminder: this is being stored in the output dir provided as a commandline parameter:
        self.runinfo['self_path'] = os.path.join(self.runinfo['output_dir'], 'RUNINFO.cp')

        if args.title:
            self.title = args.title

        filesnpaths.gen_output_directory(self.runinfo['output_dir'])


    def read_runinfo_dict(self, args):
        if args.fasta_file or args.metadata:
            raise ConfigError, "You declared a RUNINFO dict with '-r'. You are not allowed to\
                                      declare any of '-f', '-m', or '-t' parameters if you have a\
                                      RUNINFO dict. Please refer to the documentation."
 
        if not os.path.exists(args.runinfo):
            raise ConfigError, "'%s'? No such file." % (args.runinfo)

        r = dictio.read_serialized_object(args.runinfo)

        if not r.has_key('runinfo'):
            raise ConfigError, "'%s' does not seem to be a anvio RUNINFO.cp." % (args.runinfo)

        r['self_path'] = args.runinfo
        r['output_dir'] = os.path.join(os.getcwd(), os.path.dirname(args.runinfo))

        return r


    def load_from_runinfo_dict(self, args):
        if not self.runinfo.has_key('profiler_version') or self.runinfo['profiler_version'] != t.profile_db_version:
            raise ConfigError, "RUNINFO.cp seems to be generated from an older version of anvio\
                                           profiler that is not compatible with the current interactive interface\
                                           anymore. You need to re-run anvio profiler on these projects."

        # load views from the profile database
        self.load_views()
        self.default_view = self.p_meta['default_view']

        # if the user wants to see available views, show them and exit.
        if args.show_views:
            run.warning('', header = 'Available views (%d)' % len(self.views), lc = 'green')
            for view in self.views:
                run.info(view,
                         'Via "%s" table' % self.views[view]['table_name'],
                         lc='crimson',
                         mc='green' if view == self.default_view else 'crimson')
            print
            sys.exit()


        # if the user specifies a view, set it as default:
        if args.view:
            if not args.view in self.views:
                raise ConfigError, "The requested view ('%s') is not available for this run. Please see\
                                          available views by running this program with --show-views flag." % args.view

            self.p_meta = args.view

        # set clusterig
        self.runinfo['clusterings'] = self.clusterings 
        if args.tree:
            entry_id = os.path.basename(args.tree).split('.')[0]
            run.info('Additional Tree', "'%s' has been added to available trees." % entry_id)
            self.runinfo['clusterings'][entry_id] = {'newick': open(os.path.abspath(args.tree)).read()}

        # is summary being overwritten?
        if args.summary_index:
            run.info('Warning', "The default summary index in RUNINFO is being overriden by '%s'." % args.summary_index)
            self.runinfo['profile_summary_index'] = os.path.abspath(args.summary_index)
        self.splits_summary_index = dictio.read_serialized_object(self.P(self.runinfo['profile_summary_index']))

        # set title
        if args.title:
            self.title = args.title + ' (%s)' % self.default_view
        else:
            self.title = self.runinfo['sample_id'] + ' (%s)' % self.default_view


    def check_names_consistency(self):
        if self.args.skip_check_names:
            return

        splits_in_tree = sorted(self.split_names_ordered)
        splits_in_metadata = sorted(self.views[self.default_view]['dict'].keys())
        splits_in_database = sorted(self.split_sequences)

        try:
            assert(splits_in_database == splits_in_tree == splits_in_metadata)
        except:
            S = lambda x, y: "agrees" if x == y else "does not agree"
            raise ConfigError, "Entries found in the annotation database, the tree file and the\
                                      metadata need to match perfectly. It seems it is not the\
                                      case for the input you provided (the metadata %s with the tree,\
                                      the tree %s with the database, the database %s with the metadata;\
                                      HTH!)." % (S(splits_in_metadata, splits_in_tree),
                                                 S(splits_in_database, splits_in_tree),
                                                 S(splits_in_metadata, splits_in_database))

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


    def update_runinfo_on_disk(self):
        path = self.runinfo.pop('self_path')
        dictio.write_serialized_object(self.runinfo, path)


    def convert_metadata_into_json(self):
        '''This function's name must change to something more meaningful.'''

        if self.annotation_db_path:
            # FIXME: Gotta think about more carefully;
            self.args.simplify_taxonomy = False
            if self.args.simplify_taxonomy:
                for split_name in self.genes_in_splits_summary_dict:
                    s = self.genes_in_splits_summary_dict[split_name]
                    if s['taxonomy']:
                        s['taxonomy'] = s['taxonomy'].split()[0]

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
                json_header.extend(self.hmm_searches_header)

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

                # (6) hmm stuff
                if self.hmm_searches_dict:
                    json_entry.extend([self.hmm_searches_dict[split_name][header] if self.hmm_searches_dict.has_key(split_name) else None for header in self.hmm_searches_header])

                # (7) send it along!
                json_object.append(json_entry)

            self.views[view] = json_object


    def end(self):
        # FIXME: remove temp files and stuff
        pass
