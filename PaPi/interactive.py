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
from ete2 import Tree

import PaPi.db
import PaPi.profiler
import PaPi.annotation
import PaPi.fastalib as u
import PaPi.utils as utils
import PaPi.dictio as dictio
import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths
from PaPi.constants import levels_of_taxonomy


progress = terminal.Progress()
run = terminal.Run()


class InputHandler:
    """The class that loads everything for the interactive interface. Wow. Such glory."""
    def __init__(self, args):
        self.args = args
        self.runinfo = {}
        self.title = 'Unknown Project'
        self.annotation = None
        self.contig_names_ordered = None
        self.contig_rep_seqs = {}
        self.contigs_summary_index = {}
        self.contig_lengths = {}
        self.taxonomy = None
        self.taxonomic_level = None
        self.additional_metadata_path = None

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

        # take care of taxonomy if a file is provided. otherwise self.taxonomy remain None
        if args.taxonomy:
            if not args.taxonomic_level:
                raise utils.ConfigError, "When a taxonomy file is declared, a taxonomic level to visualize\
                                               has to be defined."

            self.taxonomic_level = args.taxonomic_level.lower()

            if self.taxonomic_level not in levels_of_taxonomy:
                raise utils.ConfigError, "The taxonomic level '%s' is not a valid one. Taxonomic level has\
                                               to match one of these: %s" % (args.taxonomic_level,
                                                                             levels_of_taxonomy)

            filesnpaths.is_file_tab_delimited(args.taxonomy)
            self.taxonomy = utils.get_TAB_delimited_file_as_dictionary(args.taxonomy, expected_fields = levels_of_taxonomy)

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

        self.runinfo['splits_fasta'] = os.path.abspath(args.fasta_file)
        self.runinfo['output_dir'] = os.path.abspath(args.output_dir)
        self.runinfo['default_view'] = 'single'
        self.runinfo['views'] = {'single': os.path.abspath(args.metadata)}
        self.runinfo['default_clustering'] = 'default'
        self.runinfo['clusterings'] = {'default': os.path.abspath(args.tree)}
        self.runinfo['tree_file'] = os.path.abspath(args.tree)

        if args.summary_index:
            self.runinfo['profile_summary_index'] = os.path.abspath(args.summary_index)
            self.contigs_summary_index = dictio.read_serialized_object(self.runinfo['profile_summary_index'])

        # sanity of the metadata
        filesnpaths.is_file_tab_delimited(self.runinfo['views']['single'])
        if not open(self.runinfo['views']['single']).readline().split('\t')[0] == "contig":
            raise utils.ConfigError, "The first row of the first column of the metadata file must\
                                      say 'contig', which is not the case for your metadata file\
                                      ('%s'). Please make sure this is a properly formatted metadata\
                                      file." % (self.runinfo['views']['single'])
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

        # if the user wants to see available views, show them and exit.
        if args.show_views:
            num_views = len(self.runinfo['views'])
            print "* %d view%s available for this run is listed below." % (num_views,
                                                                           's' if num_views > 1 else '')
            run.info('Available views', None, header = True)
            for view in self.runinfo['views']:
                run.info(view, 'Via "%s" table' % self.runinfo['views'][view])
            print
            sys.exit()

        # if the user specifies a view, set it as default:
        if args.view:
            if not args.view in self.runinfo['views']:
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

        # FIXME: how it came to this point is a long story. but here we are trying to comply with the design
        # that was in place by replacing 'table names' for each self.runinfo['views'] item with TAB delimited
        # file paths on this ... by .. literally .. rading them from the database, and converting them into
        # TAB delimited files, JUST TO BE READ later on as JSON formatted objects. This is ridiculous, because
        # what is being done here is 2 steps with a lot of IO when it could be one very fast step. This
        # has to be fixed at some point (tl;dr, for each view, self.runinfo['views'][view] must contain a JSON
        # formatted output of the 'view' table):
        for view in self.runinfo['views']:
            table = self.runinfo['views'][view]
            table_rows = self.profile_db.get_all_rows_from_table(table)
            tmp_file_path = filesnpaths.get_temp_file_path()
            table_structure = self.profile_db.get_table_structure(table)
            utils.store_array_as_TAB_delimited_file(table_rows, tmp_file_path, table_structure)
            self.runinfo['views'][view] = tmp_file_path

        self.contigs_summary_index = dictio.read_serialized_object(self.P(self.runinfo['profile_summary_index']))

        self.runinfo['self_path'] = args.runinfo


    def check_names_consistency(self):
        contigs_in_tree = sorted(self.contig_names_ordered)
        contigs_in_metadata = sorted([l.split('\t')[0] for l in open(self.runinfo['views'][self.runinfo['default_view']]).readlines()[1:]])
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

        if self.taxonomy:
            contigs_in_taxonomy = set(self.taxonomy.keys())
            for contig in contigs_in_metadata:
                if not contig in contigs_in_taxonomy:
                    raise utils.ConfigError, "Contig names in taxonomy do not include all contig names present\
                                                   in other files. Bad news :/"

        if self.additional_metadata_path:
            contigs_in_additional_metadata = set(sorted([l.split('\t')[0] for l in open(self.additional_metadata_path).readlines()[1:]]))
            for contig in contigs_in_tree:
                if contig not in contigs_in_additional_metadata:
                    raise utils.ConfigError, "Contig names in additional metadata file do not include all contig names present\
                                                   in other files. Bad news :/"


    def update_runinfo_on_disk(self):
        path = self.runinfo.pop('self_path')
        dictio.write_serialized_object(self.runinfo, path)


    def convert_metadata_into_json(self):
        metadata_file_path = self.runinfo['views'][self.runinfo['default_view']]

        if self.taxonomy:
            metadata_headers = utils.get_columns_of_TAB_delim_file(metadata_file_path)
            metadata_dict = utils.get_TAB_delimited_file_as_dictionary(metadata_file_path)
            for contig in metadata_dict.keys():
                metadata_dict[contig][self.taxonomic_level] = self.taxonomy[contig][self.taxonomic_level]

            new_metatada_file_path = filesnpaths.get_temp_file_path()
            headers = ['contigs', self.taxonomic_level] + metadata_headers
            utils.store_dict_as_TAB_delimited_file(metadata_dict, new_metatada_file_path, headers)
            metadata_file_path = new_metatada_file_path


        if self.annotation:
            splits_additional_info = {}
            progress.new('Generating splits additional info')
            num_contigs = len(self.contig_names_ordered)

            metadata_headers = utils.get_columns_of_TAB_delim_file(metadata_file_path)
            metadata_dict = utils.get_TAB_delimited_file_as_dictionary(metadata_file_path)

            splits = self.annotation.db.get_table_as_dict('splits', PaPi.annotation.splits_table_structure)
            for split in metadata_dict:
                for key in PaPi.annotation.splits_table_structure[1:]:
                    metadata_dict[split][key] = splits[split][key]

            new_metatada_file_path = filesnpaths.get_temp_file_path()
            headers = ['contigs'] + PaPi.annotation.splits_table_structure[1:] + metadata_headers
            utils.store_dict_as_TAB_delimited_file(metadata_dict, new_metatada_file_path, headers)
            metadata_file_path = new_metatada_file_path


        if self.additional_metadata_path:
            metadata_headers = utils.get_columns_of_TAB_delim_file(metadata_file_path)
            additional_headers = utils.get_columns_of_TAB_delim_file(self.additional_metadata_path)
            metadata_dict = utils.get_TAB_delimited_file_as_dictionary(metadata_file_path)
            metadata_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_metadata_path, dict_to_append = metadata_dict)
            
            new_metatada_file_path = filesnpaths.get_temp_file_path()
            headers = ['contigs'] + metadata_headers + additional_headers
            utils.store_dict_as_TAB_delimited_file(metadata_dict, new_metatada_file_path, headers)
            metadata_file_path = new_metatada_file_path


        json_obj = utils.get_json_obj_from_TAB_delim_metadata(metadata_file_path)

        temp_file_path = filesnpaths.get_temp_file_path()
        f = open(temp_file_path, 'w')
        f.write(json_obj)
        f.close()
        self.runinfo['metadata_json'] = temp_file_path


    def end(self):
        # FIXME: remove temp files and stuff
        pass