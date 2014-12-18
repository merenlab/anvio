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
from PaPi.utils import get_FASTA_file_as_dictionary as get_dict_f
from PaPi.utils import store_dict_as_TAB_delimited_file as store_dict


class Parser(object):
    def __init__(self, annotation_source, input_file_paths, files_expected = {}, files_structure = {}):
        self.annotation_source = annotation_source
        self.input_file_paths = input_file_paths
        self.files_expected = files_expected
        self.files_structure = files_structure
        self.input_file_names = [os.path.basename(p) for p in input_file_paths]
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
                                     % (self.annotation_source,
                                        ', '.join(self.files_expected))

            raise ConfigError, "%s parser requires %d files (%s). %s missing from your input: %s"\
                                     % (self.annotation_source,
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
            f = self.files_structure[alias]
            if f.has_key('type'):
                if f['type'] == 'fasta':
                    self.dicts[alias] = get_dict_f(self.paths[alias])
                else:
                    raise ConfigError, "Parser class does not know about file type '%s' :/" % f['type']
            else:
                # then it is tab-delimited
                indexing_field = f['indexing_field'] if f.has_key('indexing_field') else 0
                self.dicts[alias] = get_dict(self.paths[alias],
                                             column_names = self.files_structure[alias]['col_names'],
                                             column_mapping = self.files_structure[alias]['col_mapping'],
                                             indexing_field = indexing_field)


    def store_annotations(self, contigs_fasta, annotations_dict, split_length, output_file_prefix = "ANNOTATION"):
        if output_file_prefix.lower().endswith('.txt'):
            output_file_prefix = output_file_prefix[:-4]

        A = annotation.Annotation(output_file_prefix + '.db')
        A.create_new_database(contigs_fasta, annotations_dict, split_length, parser=self.annotation_source)


class MyRastCMDLine(Parser):
    def __init__(self, contigs_fasta, input_file_paths, output_file_prefix, split_length = 20000):
        files_expected = {'functions': 'functions.tbl', 'genes': 'genes.peg'}

        files_structure = {'functions': 
                                {'col_names': ['t_species', 'field2', 'prot', 'function'],
                                 'col_mapping': [str, int, str, str],
                                 'indexing_field': 2},
                           'genes': 
                                {'type': 'fasta'},}

        Parser.__init__(self, 'MyRastCMDLine', input_file_paths, files_expected, files_structure)

        annotations_dict = self.get_annotations_dict()
        self.store_annotations(contigs_fasta, annotations_dict, split_length, output_file_prefix)


    def get_annotations_dict(self):
        # start with the fasta dict to identify start and stop. fasta file looked like this, as a reminder:
        #
        #    >prot_00001 D23-1_contig_1_127_1215
        #    MTDCSXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #    >prot_00002 D23-1_contig_1_1281_1499
        #    MKDNAERKAKRRIFLXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #    >prot_00003 D23-1_contig_1_2630_1626
        #    MAKQKIRIRLKAYDHRVIDQSAEKIVETAKRSGADVSGPIPLPTE
        #    >prot_00004 D23-1_contig_1_2987_3376
        #    MKNGPKLSLALIGIFLILCEFFYGIPFLGATFILSFGWQPLIFNA
        #    >prot_00005 D23-1_contig_1_4901_3567
        #    MKNYFQFDKYGTNFKREILGGITTFLSMAYILAVNPQVLSLAGVK
        #    >prot_00006 D23-1_contig_1_7149_5014
        #    MKSLILAEKPSVARDIADALQINQKRNGYFENNQYIVTWALGHLV
        #    >prot_00007 D23-1_contig_1_7266_8231
        #    MLISLLTFISVEILYNKSNKKYGGNDMSIVQLYDITQIKSFIEHS

        annotations_dict = {}

        for key in self.dicts['genes']:
            entry = {}
            for field in annotation.annotation_table_structure:
                entry[field] = None

            prot, remainder = key.split()
            contig = '_'.join(remainder.split('_')[:-2])
            start, stop = [t for t in remainder.split('_')[-2:]]
            entry['start'] = int(start)
            entry['stop'] = int(stop)
            entry['contig'] = contig

            annotations_dict[prot] = entry

        for prot in self.dicts['functions']:
            d = self.dicts['functions'][prot]

            t_species_str = d['t_species']
            if len(t_species_str.split()) > 2:
                t_species_str = ' '.join(t_species_str.split()[0:2])
            annotations_dict[prot]['t_species'] = t_species_str
            annotations_dict[prot]['function'] = d['function']

        return annotations_dict


class MyRastGUI(Parser):
    def __init__(self, contigs_fasta, input_file_paths, output_file_prefix, split_length = 20000):
        files_expected = {'functions': 'functions.tbl', 'gene_otus': 'gene_otus.tbl', 'peg': 'peg.tbl'}

        files_structure = {'functions': 
                                {'col_names': ['prot', 'figfam', 'field3', 'field4', 'field5', 'function'],
                                 'col_mapping': [str, str, int, int, int, str]},
                           'gene_otus': 
                                {'col_names': ['prot', 't_species'],
                                 'col_mapping': None},
                           'peg':
                                {'col_names': ['prot', 'contig', 'start', 'stop'],
                                 'col_mapping': [str, str, int, int]},}

        Parser.__init__(self, 'MyRastGUI', input_file_paths, files_expected, files_structure)

        annotations_dict = self.get_annotations_dict()
        self.store_annotations(contigs_fasta, annotations_dict, split_length, output_file_prefix)


    def get_annotations_dict(self):
        proteins = set([])
        for alias in self.dicts:
            for prot in self.dicts[alias].keys():
                proteins.add(prot)

        annotations_dict = {}

        for prot in proteins:
            entry = {}
            for key in annotation.annotation_table_structure:
                entry[key] = None

            d = self.dicts['peg'][prot]
            if d['start'] < d['stop']:
                entry['start'] = d['start']
                entry['stop'] = d['stop']
                entry['direction'] = 'f'
            else:
                entry['start'] = d['stop']
                entry['stop'] = d['start']
                entry['direction'] = 'r'
            entry['contig'] = d['contig']

            if self.dicts['gene_otus'].has_key(prot):
                d = self.dicts['gene_otus'][prot]
                t_species_str = d['t_species']
                if len(t_species_str.split()) > 2:
                    t_species_str = ' '.join(t_species_str.split()[0:2])
                entry['t_species'] = t_species_str

            if self.dicts['functions'].has_key(prot):
                d = self.dicts['functions'][prot]
                entry['figfam'] = d['figfam']
                entry['function'] = d['function']
            else:
                entry['figfam'] = None
                entry['function'] = None

            annotations_dict[prot] = entry

        return annotations_dict


parsers = {"myrast_gui": MyRastGUI,
           "myrast_cmdline": MyRastCMDLine}