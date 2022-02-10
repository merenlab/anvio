#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import csv
import pandas as pd

import anvio
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.parsers.base import Parser


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


class signalp6(Parser):
    def __init__(self, input_file_paths, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress
        self.just_do_it = False

        input_file_path = self.fix_input_file(input_file_paths[0])

        files_expected = {'signalp6_output': input_file_path}

        files_structure = {'signalp6_output':
                                {'col_names': ['gene_callers_id', 'Prediction', 'OTHER', 'SP(Sec/SPI)', 'LIPO(Sec/SPII)', 'TAT(Tat/SPI)', 'TATLIPO(Sec/SPII)', 'PILIN(Sec/SPIII)', 'CS Position'],
                                 'col_mapping': [int, str, float, float, float, float, float, float, str],
                                 'indexing_field': -1,
                                 'separator': '\t'},
                            }

        self.progress.new('Initializing the parser')
        self.progress.update('...')
        Parser.__init__(self, 'signalp6', [input_file_path], files_expected, files_structure)
        self.progress.end()

        # This is where I would write specific sanity checks for signalp6


    def fix_input_file(self, input_file_path):
        """Select columns for anvio and remove duplicate rows"""
        self.progress.new('Making signalp6 output anvio friendly')
        self.progress.update('...')

        temp_file_path = filesnpaths.get_temp_file_path()

        file = open(input_file_path, "r")
        tsv_file = csv.reader(file, delimiter="\t")

        lists_from_csv = []
        for row in tsv_file:
            lists_from_csv.append(row)

        # Remove first line
        del lists_from_csv[0]

        # Get column names
        col_names = lists_from_csv[0]
        col_names[0] = col_names[0][2:]

        # create df
        df = pd.DataFrame(lists_from_csv[1:], columns = col_names)
        df = df.rename(columns={'ID': 'gene_callers_id'})
        df = df.drop_duplicates(subset=["gene_callers_id"])

        # remove "Other" annotations which means no signal peptide  
        df = df[df['Prediction'] != "OTHER"]

        df.to_csv(temp_file_path, sep = '\t', index = False, na_rep = 'NA')

        self.progress.end()

        return temp_file_path


    def get_dict(self):
        """Convert angostos output into functions dict"""
        d = self.dicts['signalp6_output']

        # Parse signalp6 output to make functions_dict
        df = pd.DataFrame.from_dict(d, orient='index')
        df['source'] = "signalp6"
        df['e_value'] = 0
        df['accession'] = df['Prediction']
        df['function'] = df['Prediction']
        df = df.drop_duplicates(subset=['gene_callers_id', 'source', 'accession', 'function'])
        d = df.to_dict(orient='index')

        return d
