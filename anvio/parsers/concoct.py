#!/usr/bin/env python

"""
    Module to parse CONCOCT clusters.
"""

from anvio.parsers.base import Parser


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class CONCOCT(Parser):
    def __init__(self, input_files, contigs='False'):
        if type(input_files) != type(list()):
            input_files = [input_files]

        files_expected = {'clusters': input_files[0]}

        files_structure = {'clusters':
                                {'col_names': ['split', 'bin_name'],
                                 'col_mapping': [str, str],
                                 'separator': ',',
                                 'indexing_field': -1,
                                 'no_header': True
                                 },
                           }
        Parser.__init__(self, 'CONCOCT', input_files, files_expected, files_structure)


    def get_clusters_dict(self):
        c = self.dicts['clusters']

        clusters_dict = {}
        for bin_name in set([e['bin_name'] for e in list(c.values())]):
            clusters_dict[bin_name] = []

        for entry in list(c.values()):
            clusters_dict[entry['bin_name']].append(entry['split'])

        return clusters_dict

