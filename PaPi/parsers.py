#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

"""
    This is the file that keeps all the parser classes for different sources of annotation.
"""

import os

import PaPi.filesnpaths as filesnpaths
import PaPi.annotation as annotation
from PaPi.utils import ConfigError
from PaPi.utils import get_TAB_delimited_file_as_dictionary as get_dict
from PaPi.utils import store_dict_as_TAB_delimited_file as store_dict


class Parser(object):
    def __init__(self, name, input_file_paths, files_expected = {}, files_structure = {}, output_file_prefix = "ANNOTATION"):
        self.name = name
        self.output_file_prefix = output_file_prefix
        self.files_structure = files_structure
        self.input_file_paths = input_file_paths
        self.input_file_names = [os.path.basename(p) for p in input_file_paths]
        self.files_expected = files_expected
        self.paths = {}
        self.dicts = {}

        if sorted(files_expected.keys()) != sorted(files_structure.keys()):
            raise ConfigError, "Items in files_expected and files_structure must match."

        missing_files = []
        for f in self.files_expected.values():
            if f not in self.input_file_names:
                missing_files.append(f)
        if missing_files:
            if sorted(missing_files) == sorted(self.files_expected.values()):
                raise ConfigError, "%s parser requires these files: %s"\
                                     % (name,
                                        ', '.join(self.files_expected))

            raise ConfigError, "%s parser requires %d files (%s). %s missing from your input: %s"\
                                     % (name,
                                        len(self.files_expected),
                                        ', '.join(self.files_expected.values()),
                                        "These files were" if len(missing_files) > 1 else "This file was",
                                        ", ".join(missing_files))

        for alias in self.files_expected:
            for i in range(0, len(self.input_file_names)):
                file_name = self.input_file_names[i]
                if self.files_expected[alias] == file_name:
                    self.paths[alias] = self.input_file_paths[i]

        for alias in self.files_expected:
            self.dicts[alias] = get_dict(self.paths[alias],
                                         column_names = self.files_structure[alias]['col_names'],
                                         column_mapping = self.files_structure[alias]['col_mapping'])

    def store_annotations(self):
        if self.output_file_prefix.lower().endswith('.txt'):
            self.output_file_prefix = self.output_file_prefix[:-4]

        A = annotation.Annotation()
        A.init_from_matrix(source = self.annotations_dict)

        A.store_annotation_matrix(self.output_file_prefix + '.txt')
        A.store_annotation_dict(self.output_file_prefix + '.cp')


class MyRast(Parser):
    def __init__(self, input_file_paths, output_file_prefix):
        files_expected = {'functions': 'functions.tbl', 'gene_otus': 'gene_otus.tbl', 'peg': 'peg.tbl'}

        files_structure = {'functions': 
                                {'col_names': ['prot', 'figfam', 'field3', 'field4', 'field5', 'function'],
                                 'col_mapping': [str, str, int, int, int, str]},
                           'gene_otus': 
                                {'col_names': ['prot', 'taxonomy'],
                                 'col_mapping': None},
                           'peg':
                                {'col_names': ['prot', 'contig', 'start', 'end'],
                                 'col_mapping': [str, str, int, int]},}

        Parser.__init__(self, 'MyRAST', input_file_paths, files_expected, files_structure, output_file_prefix)

        self.annotations_dict = self.get_annotations_dict()
        self.store_annotations()


    def get_annotations_dict(self):
        proteins = set([])
        for alias in self.dicts:
            for prot in self.dicts[alias].keys():
                proteins.add(prot)

        annotations_dict = {}

        for prot in proteins:
            entry = {}
            d = self.dicts['peg'][prot]
            if d['start'] < d['end']:
                entry['start'] = d['start']
                entry['end'] = d['end']
                entry['direction'] = 'f'
            else:
                entry['start'] = d['end']
                entry['end'] = d['start']
                entry['direction'] = 'r'
            entry['contig'] = d['contig']

            if self.dicts['gene_otus'].has_key(prot):
                d = self.dicts['gene_otus'][prot]
                taxonomy_str = d['taxonomy']
                if len(taxonomy_str.split()) > 2:
                    taxonomy_str = ' '.join(taxonomy_str.split()[0:2])
                entry['taxonomy'] = taxonomy_str
            else:
                entry['taxonomy'] = None

            if self.dicts['functions'].has_key(prot):
                d = self.dicts['functions'][prot]
                entry['figfam'] = d['figfam']
                entry['function'] = d['function']
            else:
                entry['figfam'] = None
                entry['function'] = None

            annotations_dict[prot] = entry

        return annotations_dict


parsers = {"myrast": MyRast}