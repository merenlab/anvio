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


from PaPi.parsers.base import Parser


class DefaultMatrix(Parser):
    def __init__(self, input_file_paths, annotation_table_structure):
        matrix_txt = input_file_paths[0]
        files_expected = {'matrix': matrix_txt}

        files_structure = {'matrix': 
                                {'col_names': ['prot', 'contig', 'start', 'stop', 'direction', 'figfam', 'function', 't_phylum', 't_class', 't_order', 't_family', 't_genus', 't_species'],
                                 'col_mapping': [str, str, int, int, str, str, str, str, str, str, str, str, str],
                                 }
                          }

        self.annotation_table_structure = annotation_table_structure
        Parser.__init__(self, 'DefaultMatrix', [matrix_txt], files_expected, files_structure)


    def get_annotations_dict(self):
        return self.dicts['matrix']
