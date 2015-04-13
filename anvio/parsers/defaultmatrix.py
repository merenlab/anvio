# -*- coding: utf-8

"""Parser for simple matrix (see col_names for specification)"""

from anvio.parsers.base import Parser


class DefaultMatrix(Parser):
    def __init__(self, input_file_paths, genes_table_structure):
        matrix_txt = input_file_paths[0]
        files_expected = {'matrix': matrix_txt}

        files_structure = {'matrix': 
                                {'col_names': ['prot', 'contig', 'start', 'stop', 'direction', 'figfam', 'function', 't_phylum', 't_class', 't_order', 't_family', 't_genus', 't_species'],
                                 'col_mapping': [str, str, int, int, str, str, str, str, str, str, str, str, str],
                                 }
                          }

        self.genes_table_structure = genes_table_structure
        Parser.__init__(self, 'DefaultMatrix', [matrix_txt], files_expected, files_structure)


    def get_annotations_dict(self):
        return self.dicts['matrix']
