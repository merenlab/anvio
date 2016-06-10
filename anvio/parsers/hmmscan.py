# -*- coding: utf-8

"""Parser for HMMer's hmmscan output"""

from anvio.parsers.base import Parser


class HMMScan(Parser):
    def __init__(self, hmm_scan_hits_txt):
        files_expected = {'hits': hmm_scan_hits_txt}

        files_structure = {'hits':
                                {'col_names': ['gene_name', 'gene_hmm_id', 'gene_callers_id', 'f', 'e_value', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f'],
                                 'col_mapping': [str, str, int, str, float, str, str, str, str, str, str, str, str, str, str, str, str, str],
                                 'indexing_field': -1,
                                 'no_header': True
                                 },
                        }

        Parser.__init__(self, 'HMMScan', [hmm_scan_hits_txt], files_expected, files_structure)


    def get_search_results(self):
        annotations_dict = {}

        # this is the stuff we are going to try to fill with this:
        # search_table_structure = ['entry_id', 'source', 'search_type', 'contig', 'gene_callers_id' 'gene_name', 'gene_hmm_id', 'e_value']

        entry_id = 0
        for hit in self.dicts['hits'].values():
            entry = {'entry_id': entry_id,
                     'gene_name': hit['gene_name'],
                     'gene_hmm_id': hit['gene_hmm_id'],
                     'gene_callers_id': hit['gene_callers_id'],
                     'e_value': hit['e_value']}

            entry_id += 1
            annotations_dict[entry_id] = entry

        return annotations_dict

