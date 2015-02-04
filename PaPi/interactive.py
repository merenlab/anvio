#!/usr/bin/env python
# -*- coding: utf-8

"""
Copyright (C) 2015, PaPi Authors

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.

Please read the COPYING file.
"""

import os
import sys
from ete2 import Tree

import PaPi.db
import PaPi.profiler
import PaPi.annotation
import PaPi.fastalib as u
import PaPi.utils as utils
import PaPi.dictio as dictio
import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths


progress = terminal.Progress()
run = terminal.Run()


class InputHandler:
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args):
        self.args = args
        self.views = {}
        self.runinfo = {}
        self.title = 'Unknown Project'
        self.annotation = None
        self.contig_names_ordered = None
        self.contig_rep_seqs = {}
        self.contigs_summary_index = {}
        self.contig_lengths = {}
        self.additional_metadata_path = None

        # from annotation db if exists (these will be populated in self.init_annotation_db()):
        self.annotation_dict = {}
        self.genes_in_splits = {}
        self.split_to_genes_in_splits_ids = {} # for fast access to all self.genes_in_splits entries for a given split

        self.P = lambda x: os.path.join(self.runinfo['output_dir'], x)

        self.cwd = os.getcwd()

        self.state = args.state

        # here is where the big deal stuff takes place:
        if args.runinfo:
            self.load_from_runinfo_dict(args)
        else:
            self.load_from_files(args.fasta_file, args.metadata, args.tree, args.output_dir)

        # if we have an annotation_db, lets load it up. the split length used when the annotation db
        # created must match with the split length used to profile these merged runs. the problem is,
        # if interactive binning is being called without a runinfo, we are going to have to ask the
        # the user to run interactive interface without an annotation db.
        if args.annotation_db:
            self.init_annotation_db(args.annotation_db)


        if args.additional_metadata:
            filesnpaths.is_file_tab_delimited(args.additional_metadata)
            self.additional_metadata_path = args.additional_metadata


        tree = Tree(self.P(self.runinfo['tree_file']))
        self.contig_names_ordered = [n.name for n in tree.get_leaves()]

        self.check_names_consistency()
        self.convert_metadata_into_json()

        fasta = u.SequenceSource(self.P(self.runinfo['splits_fasta']))

        while fasta.next():
            self.contig_rep_seqs[fasta.id] = fasta.seq
            self.contig_lengths[fasta.id] = len(fasta.seq)


    def init_annotation_db(self, db_path):
        self.annotation = PaPi.annotation.Annotation(db_path)
        self.annotation.init_database()

        profiling_split_length = int(self.runinfo['split_length'])
        annotation_split_length = int(self.annotation.db.get_meta_value('split_length'))
        if profiling_split_length != annotation_split_length:
            raise utils.ConfigError, "The split length (-L) used to profile these merged runs (which is '%s') seem\
                                      to differ from the split length used to generate %s (which is\
                                      '%s'). Probably the best option is to re-create the annotation database with\
                                      an identical split length parameter." % (terminal.pretty_print(profiling_split_length),
                                                                               os.path.basename(db_path),
                                                                               terminal.pretty_print(annotation_split_length))

        self.annotation_dict = self.annotation.db.get_table_as_dict('annotation')
        self.genes_in_splits = self.annotation.db.get_table_as_dict('genes_in_splits')
        for entry_id in self.genes_in_splits:
            split_name = self.genes_in_splits[entry_id]['split']
            if split_name in self.split_to_genes_in_splits_ids:
                self.split_to_genes_in_splits_ids[split_name].add(entry_id)
            else:
                self.split_to_genes_in_splits_ids[split_name] = set([entry_id])

        annotation_source = self.annotation.db.get_meta_value('annotation_source')
        run.info('annotation_db initialized', '%s (v. %s) (via "%s")' % (db_path,
                                                                         self.annotation.db.version,
                                                                         annotation_source))


    def load_from_files(self, fasta_file, metadata, tree, output_dir):
        if (not fasta_file) or (not metadata) or (not tree) or (not output_dir):
            raise utils.ConfigError, "If you do not have a RUNINFO dict, you must declare each of\
                                           '-f', '-m', '-t' and '-o' parameters. Please see '--help' for\
                                           more detailed information on them."

        if args.view:
            raise utils.ConfigError, "You can't use '-v' parameter when this program is not called with a RUNINFO.cp"

        if args.show_views:
            raise utils.ConfigError, "Sorry, there are no views to show when there is no RUNINFO.cp :/"

        metadata_path = os.path.abspath(metadata)
        self.runinfo['splits_fasta'] = os.path.abspath(args.fasta_file)
        self.runinfo['output_dir'] = os.path.abspath(args.output_dir)
        self.runinfo['default_view'] = 'single'
        self.runinfo['default_clustering'] = 'default'
        self.runinfo['clusterings'] = {'default': os.path.abspath(args.tree)}
        self.runinfo['tree_file'] = os.path.abspath(args.tree)

        if args.summary_index:
            self.runinfo['profile_summary_index'] = os.path.abspath(args.summary_index)
            self.contigs_summary_index = dictio.read_serialized_object(self.runinfo['profile_summary_index'])

        # sanity of the metadata
        filesnpaths.is_file_tab_delimited(metadata_path)
        metadata_columns = utils.get_columns_of_TAB_delim_file(metadata_path, include_first_column=True)
        if not metadata_headers[0] == "contig":
            raise utils.ConfigError, "The first row of the first column of the metadata file must\
                                      say 'contig', which is not the case for your metadata file\
                                      ('%s'). Please make sure this is a properly formatted metadata\
                                      file." % (metadata_path)

        # store metadata as view:
        self.views['single'] = {'header': metadata_headers,
                                'dict': get_TAB_delimited_file_as_dictionary(metadata_path)}

        filesnpaths.is_file_fasta_formatted(self.runinfo['splits_fasta'])

        # reminder: this is being stored in the output dir provided as a commandline parameter:
        self.runinfo['self_path'] = os.path.join(self.runinfo['output_dir'], 'RUNINFO.cp')

        if args.title:
            self.title = args.title

        filesnpaths.gen_output_directory(self.runinfo['output_dir'])


    def load_from_runinfo_dict(self, args):
        if args.fasta_file or args.metadata:
            raise utils.ConfigError, "You declared a RUNINFO dict with '-r'. You are not allowed to\
                                      declare any of '-f', '-m', or '-t' parameters if you have a\
                                      RUNINFO dict. Please refer to the documentation."
 
        if not os.path.exists(args.runinfo):
            raise utils.ConfigError, "'%s'? No such file." % (args.runinfo)

        self.runinfo = dictio.read_serialized_object(args.runinfo)
        self.views = self.runinfo['views']
        self.runinfo['views'] = {}

        # if the user wants to see available views, show them and exit.
        if args.show_views:
            num_views = len(self.views)
            print "* %d view%s available for this run is listed below." % (num_views,
                                                                           's' if num_views > 1 else '')
            run.info('Available views', None, header = True)
            for view in self.views:
                run.info(view, 'Via "%s" table' % self.views[view])
            print
            sys.exit()

        # if the user specifies a view, set it as default:
        if args.view:
            if not args.view in self.views:
                raise utils.ConfigError, "The requested view ('%s') is not available for this run. Please see\
                                          available views by running this program with --show-views flag." % args.view

            self.runinfo['default_view'] = args.view


        base_dir = os.path.dirname(args.runinfo)
        self.runinfo['output_dir'] = os.path.join(os.getcwd(), base_dir)

        if not self.runinfo.has_key('runinfo'):
            raise utils.ConfigError, "'%s' does not seem to be a PaPi RUNINFO.cp." % (args.runinfo)

        # connect to the PROFILE.db
        self.profile_db = PaPi.db.DB(self.P(self.runinfo['profile_db']), PaPi.profiler.__version__)


        # if there is a tree supplied through the command line, let that tree override what is in the
        # runinfo. kinda sketchy, but it is necessary.
        if args.tree:
            run.info('Warning', "The default tree in RUNINFO is being overriden by '%s'." % args.tree)
            self.runinfo['tree_file'] = os.path.abspath(args.tree)
        else:
            self.runinfo['tree_file'] = self.runinfo['clusterings'][self.runinfo['default_clustering']]

        if args.summary_index:
            run.info('Warning', "The default summary index in RUNINFO is being overriden by '%s'." % args.summary_index)
            self.runinfo['profile_summary_index'] = os.path.abspath(args.summary_index)

        if not self.runinfo.has_key('profiler_version') or self.runinfo['profiler_version'] != PaPi.profiler.__version__:
            raise utils.ConfigError, "RUNINFO.cp seems to be generated from an older version of PaPi\
                                           profiler that is not compatible with the current interactive interface\
                                           anymore. You need to re-run PaPi profiler on these projects."

        if args.title:
            self.title = args.title + ' (%s)' % self.runinfo['default_view']
        else:
            self.title = self.runinfo['sample_id'] + ' (%s)' % self.runinfo['default_view']

        # read available views from the profile database:
        for view in self.views:
            table = self.views[view]
            self.views[view] = {'header': self.profile_db.get_table_structure(table)[1:],
                                'dict': self.profile_db.get_table_as_dict(table)}

        self.contigs_summary_index = dictio.read_serialized_object(self.P(self.runinfo['profile_summary_index']))

        self.runinfo['self_path'] = args.runinfo


    def check_names_consistency(self):
        contigs_in_tree = sorted(self.contig_names_ordered)
        contigs_in_metadata = sorted(self.views[self.runinfo['default_view']]['dict'].keys())
        contigs_in_fasta = sorted(utils.get_all_ids_from_fasta(self.P(self.runinfo['splits_fasta'])))

        try:
            assert(contigs_in_fasta == contigs_in_tree == contigs_in_metadata)
        except:
            S = lambda x, y: "agrees" if x == y else "does not agree"
            raise utils.ConfigError, "Contigs name found in the FASTA file, the tree file and the\
                                           metadata needs to match perfectly. It seems it is not the\
                                           case for the input you provided (the metadata %s with the tree,\
                                           the tree %s with the fasta, the fasta %s with the metadata;\
                                           HTH!)." % (S(contigs_in_metadata, contigs_in_tree),
                                                      S(contigs_in_fasta, contigs_in_tree),
                                                      S(contigs_in_metadata, contigs_in_fasta))

        if self.additional_metadata_path:
            contigs_in_additional_metadata = set(sorted([l.split('\t')[0] for l in open(self.additional_metadata_path).readlines()[1:]]))
            for contig_name in contigs_in_additional_metadata:
                if contig_name not in contigs_in_tree:
                    raise utils.ConfigError, "Some contig names in the additional metadata file is not found in contigs\
                                              found in other files (such as this one: '%s'). Additional metadata file\
                                              does not have to list all contigs, but whenever there is a contig, it must\
                                              be present in other files. Bad news :/" % (contig_name)


    def update_runinfo_on_disk(self):
        path = self.runinfo.pop('self_path')
        dictio.write_serialized_object(self.runinfo, path)


    def convert_metadata_into_json(self):
        '''This function's name must change to something more meaningful.'''

        annotation_dict, annotation_headers = None, []
        splits_header = None
        if self.annotation:
            annotation_dict = self.annotation.db.get_table_as_dict('splits')
            annotation_headers = self.annotation.db.get_table_structure('splits')[1:]

        additional_dict, additional_headers = None, []
        if self.additional_metadata_path:
            additional_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_metadata_path)
            additional_headers = utils.get_columns_of_TAB_delim_file(self.additional_metadata_path)

        for view in self.views:
            # here we will populate runinfo['views'] with json objects.
            view_dict = self.views[view]['dict']
            view_headers = self.views[view]['header']

            json_object = []

            # set the header line:
            json_header = ['contigs']
            # first annotation, if exists
            if annotation_headers:
                json_header.extend(annotation_headers)
            # then, the view!
            json_header.extend(view_headers)
            # additional headers as the outer ring:
            if additional_headers:
                json_header.extend(additional_headers)
            # finalize it:
            json_object.append(json_header)

            for split_name in view_dict:
                json_entry = [split_name]

                json_entry.extend([annotation_dict[split_name][header] for header in annotation_headers])
                json_entry.extend([view_dict[split_name][header] for header in view_headers])
                json_entry.extend([additional_dict[split_name][header] if additional_dict.has_key(split_name) else None for header in additional_headers])

                json_object.append(json_entry)
                 
            self.runinfo['views'][view] = json_object


    def end(self):
        # FIXME: remove temp files and stuff
        pass
