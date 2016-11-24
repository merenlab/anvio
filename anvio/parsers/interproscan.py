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


from anvio.parsers.base import Parser


class InterProScan(Parser):
    def __init__(self, input_file_paths):
        input_file_path = input_file_paths[0]
        files_expected = {'matrix': input_file_path}

        files_structure = {'matrix':
                                {'col_names': ['gene_callers_id', 'hash', 'length', 'source', 'accession', 'function', 'start', 'stop', 'e_value', 'status', 'date'],
                                 'col_mapping': [int, str, int, str, str, str, int, int, str, str, str],
                                 'indexing_field': -1,
                                 'no_header': True},
                            }

        Parser.__init__(self, 'InterProScan', input_file_paths, files_expected, files_structure)


    def get_dict(self):
        d = self.dicts['matrix']
        for entry in d:
            try:
                d[entry]['e_value'] = float(d[entry]['e_value'])
            except:
                d[entry]['e_value'] = 0.0

            for key in d[entry]:
                if d[entry][key] == 'None':
                    d[entry][key] = None

            if not d[entry]['function'] and d[entry]['accession']:
                d[entry]['function'] = d[entry]['accession']

        return d
