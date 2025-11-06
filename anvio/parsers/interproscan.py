#!/usr/bin/env python
# -*- coding: utf-8

import pandas as pd

import anvio
from anvio.parsers.base import Parser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


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

        # Drop duplicate function annotations.
        df = pd.DataFrame.from_dict(d, orient='index')
        df = df.drop_duplicates(subset=['gene_callers_id', 'hash', 'source', 'accession', 'function'])
        d = df.to_dict(orient='index')

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
