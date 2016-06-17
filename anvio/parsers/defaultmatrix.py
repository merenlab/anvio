#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2016, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import anvio.terminal as terminal

from anvio.parsers.base import Parser
from anvio.parsers.base import TaxonomyHelper


class DefaultMatrix(Parser):
    def __init__(self, input_file_paths, taxonomy_table_structure, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        matrix_txt = input_file_paths[0]
        files_expected = {'matrix': matrix_txt}

        files_structure = {'matrix':
                                {'col_names': ['gene_callers_id', 't_phylum', 't_class', 't_order', 't_family', 't_genus', 't_species'],
                                 'col_mapping': [int, str, str, str, str, str, str],
                                 'only_expected_fields': True,
                                 }
                          }

        self.taxonomy_table_structure = taxonomy_table_structure
        Parser.__init__(self, 'DefaultMatrix', [matrix_txt], files_expected, files_structure)


    def process(self):
        """Parses the simple matrix file, returns two dicts: genes_taxonomy, taxon_names."""

        self.run.info('Total num hits found', len(self.dicts['matrix']))

        genes_taxonomy, taxon_names = TaxonomyHelper(self.dicts['matrix']).get_genes_taxonomy_and_taxon_names_dicts()

        return (genes_taxonomy, taxon_names)
