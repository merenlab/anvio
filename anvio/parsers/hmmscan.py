# -*- coding: utf-8

"""Parser for HMMer's hmmscan output"""

from anvio.parsers.base import Parser


class HMMScan(Parser):
    def __init__(self, proteins_in_contigs_fasta, hmm_scan_hits_txt):
        files_expected = {'proteins': proteins_in_contigs_fasta, 'hits': hmm_scan_hits_txt}

        files_structure = {'hits': 
                                {'col_names': ['gene_name', 'gene_id', 'query_name', 'f', 'e_value', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f','f', 'f','f'],
                                 'col_mapping': [str, str, str, str, float, str, str, str, str, str, str, str, str, str, str, str, str, str],
                                 'indexing_field': -1,
                                 'no_header': True
                                 },
                           'proteins': 
                                {'type': 'fasta'},}

        Parser.__init__(self, 'HMMScan', [proteins_in_contigs_fasta, hmm_scan_hits_txt], files_expected, files_structure)


    def get_search_results(self):
        annotations_dict = {}

        # this is the stuff we are going to try to fill with this:
        # search_table_structure = ['entry_id', 'source', 'search_type', 'contig',
        #                           'start' , 'stop'  , 'gene_name', 'gene_id', 'e_value']

        orfs_dict = {}
        for defline in self.dicts['proteins'].keys():
            fields = [f.strip() for f in defline.split('#')]
            orfs_dict[fields[0]] = {'contig': '_'.join(fields[0].split('_')[:-1]),
                                    'start': int(fields[1]),
                                    'stop': int(fields[2])}

        entry_id = 0
        for hit in self.dicts['hits'].values():
            orf = orfs_dict[hit['query_name']]
            entry = {'entry_id': entry_id,
                     'contig': orf['contig'],
                     'start': orf['start'],
                     'stop': orf['stop'],
                     'gene_name': hit['gene_name'],
                     'gene_id': hit['gene_id'],
                     'e_value': hit['e_value']}

            entry_id += 1
            annotations_dict[entry_id] = entry

        return annotations_dict

