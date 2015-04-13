# -*- coding: utf-8

"""Parser for files generated through MyRAST graphical user interface"""

from anvio.parsers.base import Parser


class MyRastGUI(Parser):
    def __init__(self, input_file_paths, genes_table_structure):
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

        self.genes_table_structure = genes_table_structure
        Parser.__init__(self, 'MyRastGUI', input_file_paths, files_expected, files_structure)


    def get_annotations_dict(self):
        proteins = set([])
        for alias in self.dicts:
            for prot in self.dicts[alias].keys():
                proteins.add(prot)

        annotations_dict = {}

        for prot in proteins:
            entry = {}
            for key in self.genes_table_structure:
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
