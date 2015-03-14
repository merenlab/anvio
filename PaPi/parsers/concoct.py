#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2015, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

"""
    This is the module to parse CONCOCT clusters.
"""

from PaPi.parsers.base import Parser


class CONCOCT(Parser):
    def __init__(self, input_files, contigs = 'False'):
        if type(input_files) != type(list()):
            input_files = [input_files]

        files_expected = {'clusters': input_files[0]}

        files_structure = {'clusters': 
                                {'col_names': ['split', 'cluster_id'],
                                 'col_mapping': [str, str],
                                 'separator': ',',
                                 'indexing_field': -1,
                                 'no_header': True
                                 },
                           }
        Parser.__init__(self, 'CONCOCT', input_files, files_expected, files_structure)


    def get_clusters_dict(self):
        return self.dicts['clusters']

